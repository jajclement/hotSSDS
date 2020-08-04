#!/usr/bin/env nextflow
/*
========================================================================================
                        SSDS Pipeline version 1.8_NF_pa
			Author : Kevin Brick (version 1.8_NF)
========================================================================================
 SSDS nextflow pipeline
 #### Homepage / Documentation
 https://github.com/kevbrick/SSDSnextflowPipeline
 Adapted from version 1.8_NF (Pauline Auffret, 2020)
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Single-Stranded-DNA-Sequencing (SSDS) Pipeline : Align & Parse ssDNA
Pipeline overview:
1. trimming    Trimmomatic for quality trimming, adapters removal and hard trimming
2. fastqc    FastQC for sequencing reads quality control
3. bwaAlign     Custom BWA alignment against reference genome
4. filterBam
5. parseITRs
6. gatherOutputs
7. makeDeeptoolsBigWig
8. samStats
9. toFRBigWig
10. makeSSreport
11. multiqc
----------------------------------------------------------------------------------------
*/
def helpMessage() { 
    log.info"""
=========================================
  SSDS Pipeline version 1.8_NF
=========================================
    Usage:

    nextflow run main.nf -c conf.nf

=============================================================================

Input data Arguments:                  
                              
    --fqdir     DIR     PATH TO PAIRED-END FASTQ(.GZ) DIRECTORY (i.e. /path/to/fq/*{1,2}.fq.gz)
OR  --sra_ids   STRING  SRA SAMPLE ID(s) (Comma separated list of SRA IDS, i.e ['ERR908507', 'ERR908506']
OR  --bamdir    DIR     PATH TO BAM DIRECTORY (i.e. /path/to/bam/*.bam)
    --genomedir	DIR     PATH TO GENOME DIRECTORY
    --genome	STRING	REFERENCE GENOME NAME (i.e "mm10", must correspond to a folder in --genomedir)
    --genome_fasta  FILE	PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA
    --genome2screen STRING	GENOMES TO SCREEN FOR FASTQC SCREENING (Comma separated list of genomes to screen reads for contamination, names must correspond to folders in --genomedir)
    --adapters  FILE	PATH TO ADAPTERS FILE FOR TRIMMOMATIC (special formatting see http://www.usadellab.org/cms/?page=trimmomatic)
                                                                  
Output and Tempory directory Arguments:                            
                                                                  
    --name      STRING    RUN NAME       
    --outdir    DIR       PATH TO OUTPUT DIRECTORY            
    --scratch   DIR       PATH TO TEMPORARY DIRECTORY

Pipeline dependencies:

    --src	        DIR	 PATH TO SOURCE DIRECTORY (containing perl scripts)
    --custom_bwa        EXE	PATH TO CUSTOM BWA EXEC
    --custom_bwa_ra	EXE	PATH TO CUSTOM BWA_SRA EXEC
    --custom_multiqc	EXE	PATH TO CUSTOM MULTIQC EXEC
    --hotspots	        DIR	PATH TO HOTSPOTS FILES DIRECTORY

Trimming arguments

    --trim_minlen	INT	trimmomatic : minimum length of reads after trimming (default 25)
    --trim_crop         INT	trimmomatic : Cut the read to that specified length (default 50, set to initial length of reads if you want a different crop length for R1 and R2)
    --trim_cropR1	INT	fastx : Cut the R1 read to that specified length
    --trim_cropR2	INT	fastx : Cut the R2 read to that specified length
    --trim_slidingwin	STRING	trimmomatic : perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default "4:15")
    --trim_illumina_clip	STRING	trimmomatic :  Cut adapter and other illumina-specific sequences from the read (default "2:20:10")

Bam processing arguments

    --bamPGline	STRING	bam header (default '@PG\\tID:ssDNAPipeline1.8_nxf_KBRICK')
                                                                   
=================================================================================================
""".stripIndent()
}

//Show help message
params.help = false 
if (params.help){
    helpMessage()
    exit 0
} 

// Define global variables
// The channel initialization depends on the input data type (fastq, bam or sra IDs)
def inputType
if (params.sra_ids){inputType = 'sra'}
if (params.bamdir){inputType = 'bam'}
if (params.fqdir){inputType = 'fastq'}

// Create scratch directory
scrdir = file("${params.scratch}")
scrdir.mkdirs()

// Custom variables
def outNameStem = "${params.name}.SSDS.${params.genome}"
def tmpNameStem = "${params.name}.tmpFile.${params.genome}"

// In-house scripts (@Kevin Brick) used in the pipeline
def generateFastQCScreenConfig_script = "${params.src}/generateFastQCScreenConfig.pl"
def ITR_id_v2c_NextFlow2_script = "${params.src}/ITR_id_v2c_NextFlow2.pl"
def ssDNA_to_bigwigs_FASTER_LOMEM_script = "${params.src}/ssDNA_to_bigwigs_FASTER_LOMEM.pl"
def makeSSMultiQCReport_nextFlow_script= "${params.src}/makeSSMultiQCReport_nextFlow.pl"

// ******************* //
// BEGINNING PIPELINE  //
// ******************* //
// Get input files according to the input types : sra ; bam or fastq
switch (inputType) {
case 'sra':
    Channel
        .fromSRA(params.sra_ids, apiKey:params.ncbi_api_key)
        .set { fq_ch }
break
case 'bam':
    Channel
        .fromPath(params.bamdir, checkIfExists:true)
        .set { bam_ch }
    process bam2fastq {
        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
        input:
            file(bam) from bam_ch
        output:
            set val("${bam.baseName}"), file('*.fastq.gz') into fq_ch                                    
        script:
        """
        picard FixMateInformation I=${bam} O=fixmate.bam SORT_ORDER=queryname \
                    TMP_DIR=${params.scratch} VALIDATION_STRINGENCY=LENIENT
        samtools sort -n fixmate.bam -o sorted.bam
        bedtools bamtofastq -i sorted.bam -fq ${bam.baseName}_1.fastq -fq2 ${bam.baseName}_2.fastq
        gzip ${bam.baseName}_1.fastq ; gzip ${bam.baseName}_2.fastq
        """
    }
break
case 'fastq':
    Channel
        .fromFilePairs(params.fqdir, checkIfExists:true)
        .set { fq_ch }
break
}

// TRIMMING PROCESS : USE TRIMMOMATIC TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES
process trimming {
    tag "$sampleId" 
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*_report.txt"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.zip"
    input:
        set val(sampleId), file(reads) from fq_ch
    output:
        set val("${sampleId}"), file(reads) into fqc_ch 
        set val("${sampleId}"), file('*R1.fastq.gz'), file('*R2.fastq.gz') into trim_ch
        file '*_report.txt' into trimmomaticReport
        file '*.html' into trimmedFastqcReport
        file '*.zip' into trimmedFastqcData
    script:
    """
    trimmomatic PE -threads ${task.cpus} ${reads} \
                ${sampleId}_trim_R1.fastq.gz R1_unpaired.fastq.gz \
                ${sampleId}_trim_R2.fastq.gz R2_unpaired.fastq.gz \
                ILLUMINACLIP:${params.adapters}:${params.trim_illuminaclip} SLIDINGWINDOW:${params.trim_slidingwin} MINLEN:${params.trim_minlen} CROP:${params.trim_crop} \
                >& ${sampleId}_trim_${outNameStem}_trimmomatic_report.txt 2>&1
    fastqc -t ${task.cpus} ${sampleId}_trim_R1.fastq.gz ${sampleId}_trim_R2.fastq.gz 
    """
}

// QUALITY CONTROL PROCESS : RUN FASTQC AND FASTQSCREEN
process fastqc {
    tag "$sampleId"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.zip"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.png"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.txt"
    input:
        set val(sampleId), file(reads) from fqc_ch
    output:
        file '*zip'  into fqcZip
        file '*html' into repHTML
        file '*png'  into repPNG
        file '*txt'  into repTXT
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    perl ${generateFastQCScreenConfig_script} ${params.genomes2screen} ${task.cpus} ${params.genomedir} > fastq_screen.conf 
    fastq_screen --threads ${task.cpus} --force --aligner bwa --conf fastq_screen.conf ${reads}
    """
}


// MAPPING PROCESS : USE CUSTOM BWA TO ALIGN SSDS DATA
process bwaAlign {
    tag "$sampleId"
    input:
        set val(sampleId), file(fqR1), file(fqR2) from trim_ch
    output:
        file '*.sorted.bam' into sortedbam2filterbam, sortedbam2parseITR
    script:
  """
    if [ ${params.trim_cropR1} != ${params.trim_cropR2} ]
    then
        fastx_trimmer -f 1 -l ${params.trim_cropR1} -i ${fqR1} -o ${tmpNameStem}.R1.fastq
        fastx_trimmer -f 1 -l ${params.trim_cropR2} -i ${fqR2} -o ${tmpNameStem}.R2.fastq
    else
        mv ${fqR1} ${tmpNameStem}.R1.fastq
        mv ${fqR2} ${tmpNameStem}.R2.fastq
    fi
    ${params.custom_bwa} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R1.fastq \
            > ${tmpNameStem}.R1.sai
    ${params.custom_bwa_ra} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R2.fastq \
            > ${tmpNameStem}.R2.sai
    ${params.custom_bwa} sampe ${params.genome_fasta} ${tmpNameStem}.R1.sai ${tmpNameStem}.R2.sai ${tmpNameStem}.R1.fastq \
            ${tmpNameStem}.R2.fastq >${tmpNameStem}.unsorted.sam
    picard SamFormatConverter I=${tmpNameStem}.unsorted.sam O=${tmpNameStem}.unsorted.tmpbam VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${tmpNameStem}.unsorted.tmpbam O=${sampleId}.sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${sampleId}.sorted.bam 
    
    """
}

// PROCESS : FILTER BAM FILES
process filterBam {
    tag "$bam"
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        file(bam) from sortedbam2filterbam
    output:
        file '*.unparsed.bam' into bamAligned
        file '*.unparsed.bam.bai' into bamIDX 
        file '*.unparsed.suppAlignments.bam' into bamAlignedSupp
        file '*.unparsed.suppAlignments.bam.bai' into bamIDXSupp
        file '*MDmetrics.txt' into listz
    script:
    """
    samtools view -F 2048 -hb ${bam} > ${bam.baseName}.ok.bam
    picard MarkDuplicatesWithMateCigar I=${bam.baseName}.ok.bam O=${bam.baseName}.unparsed.bam \
        PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.MDmetrics.txt MINIMUM_DISTANCE=400 \
        CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${bam.baseName}.unparsed.bam
    samtools view -f 2048 -hb ${bam} > ${bam.baseName}.supp.bam
    picard MarkDuplicatesWithMateCigar I=${bam.baseName}.supp.bam O=${bam.baseName}.unparsed.suppAlignments.bam \
        PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.suppAlignments.MDmetrics.txt \
        MINIMUM_DISTANCE=400 CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${bam.baseName}.unparsed.suppAlignments.bam
    """
}

// PROCESS : PARSE ITRS
process parseITRs {
    tag "$bam"
    input:
        file(bam) from sortedbam2parseITR
    output:                                                                                               
        file '*.md.*bam'                       into ITRBAM mode flatten
        file '*.md.*bam.bai'                   into ITRBAMIDX mode flatten
        file '*bed'                            into ITRBED mode flatten
        file '*.md.*MDmetrics.txt'             into ITRMD mode flatten
        set file('*.md.*bam'),file('*.md.*bam.bai') into BAMwithIDXfr, BAMwithIDXss, BAMwithIDXdt mode flatten
    script:
    """
    perl ${ITR_id_v2c_NextFlow2_script} ${bam} ${params.genome_fasta}
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.ssDNA_type1.bed \
        -o ${bam}.ssDNA_type1.bed
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.ssDNA_type2.bed \
        -o ${bam}.ssDNA_type2.bed
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.dsDNA.bed \
        -o ${bam}.dsDNA.bed
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.dsDNA_strict.bed \
        -o ${bam}.dsDNA_strict.bed
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.unclassified.bed \
        -o ${bam}.unclassified.bed
    samtools view -H ${bam} > header.txt
    echo -e "${params.bamPGline}" >>header.txt
    cat header.txt ${bam}.ssDNA_type1.sam  >${bam}.ssDNA_type1.RH.sam
    cat header.txt ${bam}.ssDNA_type2.sam  >${bam}.ssDNA_type2.RH.sam
    cat header.txt ${bam}.dsDNA.sam        >${bam}.dsDNA.RH.sam
    cat header.txt ${bam}.dsDNA_strict.sam >${bam}.dsDNA_strict.RH.sam
    cat header.txt ${bam}.unclassified.sam >${bam}.unclassified.RH.sam
    samtools view -Shb ${bam}.ssDNA_type1.RH.sam  >${bam}.ssDNA_type1.US.bam
    samtools view -Shb ${bam}.ssDNA_type2.RH.sam  >${bam}.ssDNA_type2.US.bam
    samtools view -Shb ${bam}.dsDNA.RH.sam        >${bam}.dsDNA.US.bam
    samtools view -Shb ${bam}.dsDNA_strict.RH.sam >${bam}.dsDNA_strict.US.bam
    samtools view -Shb ${bam}.unclassified.RH.sam >${bam}.unclassified.US.bam
    picard SortSam I=${bam}.ssDNA_type1.US.bam  O=${bam}.ssDNA_type1.bam \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${bam}.ssDNA_type2.US.bam  O=${bam}.ssDNA_type2.bam \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${bam}.dsDNA.US.bam O=${bam}.dsDNA.bam \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${bam}.dsDNA_strict.US.bam O=${bam}.dsDNA_strict.bam \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${bam}.unclassified.US.bam O=${bam}.unclassified.bam \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard MarkDuplicatesWithMateCigar I=${bam}.ssDNA_type1.bam O=${bam}.md.ssDNA_type1.bam \
        PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.md.ssDNA_type1.MDmetrics.txt  CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard MarkDuplicatesWithMateCigar I=${bam}.ssDNA_type2.bam O=${bam}.md.ssDNA_type2.bam \
        PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.md.ssDNA_type2.MDmetrics.txt  CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard MarkDuplicatesWithMateCigar I=${bam}.dsDNA.bam O=${bam}.md.dsDNA.bam \
        PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.md.dsDNA.MDmetrics.txt  CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard MarkDuplicatesWithMateCigar I=${bam}.dsDNA_strict.bam O=${bam}.md.dsDNA_strict.bam \
        PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.md.dsDNA_strict.MDmetrics.txt  CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard MarkDuplicatesWithMateCigar I=${bam}.unclassified.bam O=${bam}.md.unclassified.bam \
        PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.md.unclassified.MDmetrics.txt  CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${bam}.md.ssDNA_type1.bam
    samtools index ${bam}.md.ssDNA_type2.bam
    samtools index ${bam}.md.dsDNA.bam
    samtools index ${bam}.md.dsDNA_strict.bam
    samtools index ${bam}.md.unclassified.bam
    """
}

/* MAKE CHANNEL WITH BAM TYPES
Channel
    .from(["ssDNA_type1","ssDNA_type2","dsDNA","dsDNA_strict","unclassified"])
    .set { parseType }

// GATHER BAM OUTPUTS PROCESS
process gatherOutputs {
    tag { pType }
    publishDir params.outdir,  mode: 'copy', overwrite: false 
    input:
        val pType  from parseType
        file T1bam from bamT1.collect()
        file T2bam from bamT2.collect()
        file Dbam  from bamD.collect()
        file DSbam from bamDS.collect()
        file UNbam from bamUN.collect()
                                       
        file T1bed from bedT1.collect()
        file T2bed from bedT2.collect()
        file Dbed  from bedD.collect()
        file DSbed from bedDS.collect()
        file UNbed from bedUN.collect()
    output:
        file '*bam'                       into mergeITRBAM
        file '*bam.bai'                   into mergeITRBAMIDX
        file '*bed'                       into mergeITRBED
        file '*MDmetrics.txt'             into mergeITRMD
        set file('*bam'),file('*bam.bai') into mergeBAMwithIDXfr, mergeBAMwithIDXss, mergeBAMwithIDXdt
    script:
    """
    picard MergeSamFiles O=${outNameStem}.${pType}.US.BAM `ls *${pType}.bam | sed 's/.*\$/I=& /'`
    picard SortSam I=${outNameStem}.${pType}.US.BAM   O=${outNameStem}.${pType}.S.BAM \
        SO=coordinate VALIDATION_STRINGENCY=LENIENT
    picard MarkDuplicatesWithMateCigar I=${outNameStem}.${pType}.S.BAM O=${outNameStem}.${pType}.bam \
    PG=Picard2.9.2_MarkDuplicates M=${outNameStem}.ssDNA_type1.MDmetrics.txt  CREATE_INDEX=false \
    VALIDATION_STRINGENCY=LENIENT
    samtools index ${outNameStem}.${pType}.bam 
    sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 -m *${pType}.bed  >${outNameStem}.${pType}.bed
    rm -f ${outNameStem}.${pType}.US.BAM
    rm -f ${outNameStem}.${pType}.S.BAM 
    """
}
*/

// BIGWIG PROCESS
process makeDeeptoolsBigWig { 
    tag "$bam"
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        set file(bam), file(bamidx) from BAMwithIDXdt
    output:
        file '*deeptools.*' into dtDeepTools
    shell:
    println(BAMwithIDXdt)
    """
    bamCoverage --bam ${bam} --normalizeUsing RPKM --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.RPKM.bigwig 
    plotCoverage --bamfiles ${bam} --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.coveragePlot.png
    plotFingerprint --bamfiles ${bam} --labels ${bam} --numberOfProcessors ${task.cpus} \
        --minMappingQuality 30 --skipZeros --plotFile ${bam.baseName}.deeptools.fingerprints.png \
        --outRawCounts ${bam.baseName}.deeptools.fingerprints.tab
    """
}

// BAM STATS PROCESS
process samStats {
    tag "$bam"
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        set file(bam),file(bamidx) from BAMwithIDXss
    output:
        file '*stats.tab' into dtSamStat
    script:
        """
        samtools idxstats ${bam} > ${bam.baseName}.idxstats.tab
        samtools stats ${bam} > ${bam.baseName}.samstats.tab
        """
}

// FWD/REV bigwig PROCESS
process toFRBigWig {
    tag "$bam"
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        set file(bam),file(bamidx) from BAMwithIDXfr
    output:
        file '*.bigwig' into frBW
    script:
        """
        perl ${ssDNA_to_bigwigs_FASTER_LOMEM_script} --bam ${bam} --g ${params.genome} --o ${bam.baseName}.out--s 100 --w 1000 --sc ${params.scratch} --gd ${params.genomedir}   -v
        """
}    

// FINAL PROCESS TO MAKE MULTIQC REPORT
process makeSSreport {
    tag "$bam"
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        file ssdsBEDs from ITRBED
        file bam      from bamAligned
    output:
        file '*SSDSreport*' into SSDSreport mode flatten
    script:
        """
        perl ${makeSSMultiQCReport_nextFlow_script} ${bam} $ssdsBEDs --g ${params.genome} --h ${params.hotspots}
        """
}

// MULTIQC PROCESS
process multiqc {
    publishDir params.outdir,  mode: 'copy', overwrite: false
    input:
        file SSDSreport
    output:
        file '*ultiQC*' into multiqcOut
    script:
    """
    python ${params.custom_multiqc} -m ssds -n ${outNameStem}.multiQC .
    """
}

// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    if (workflow.success){
        def newFile = new File("${params.outdir}/${outNameStem}.ssds_ok.done")
        newFile.createNewFile()
    }
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

