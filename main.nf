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
0. makeScreenConfigFile Generate configuration file for fastqscreen in accordance with params.genome2screen
1. trimming             Runs Trimmomatic or trim_galore for quality trimming, adapters removal and hard trimming
2. fastqc               Runs FastQC for sequencing reads quality control
3. bwaAlign             Runs Custom BWA algorithm (bwa-ra : bwa right align) against reference genome
4. filterBam            Uses PicardTools and Samtools for duplicates marking ; removing supplementary alignements ; bam sorting and indexing
5. parseITRs            Uses in-house perl script and samtools to parse BAM files into type1 ss dna, type 2 ss dna, ds dna, strict ds dna, unclassfied (5 types bam files)
6. makeDeeptoolsBigWig  Runs deeptool to make bigwig files for each 5 types bam files
7. samStats             Runs samtools to make stats on 5 types bam files
8. toFRBigWig           Uses in-house perl to make FWD bigwig files
9. makeSSreport         Uses in-house perl script to parse bed files
10. ssds_multiqc        Runs custom multiqc to edit report on SSDS alignement
11. general_multiqc     Runs multiqc to edit a general stat report
----------------------------------------------------------------------------------------
*/
def helpMessage() { 
    log.info"""
=========================================
  SSDS Pipeline version 1.8_NF_pa
=========================================
    Usage:

    nextflow run main.nf -c conf/igh.config --fqdir tests/fastq/*{R1,R2}.fastq --name "runtest" --trim_cropR1 36 --trim_cropR2 40 --with_trimgalore -profile conda -resume


=============================================================================

Input data parameters:                  
    --fqdir     		DIR     PATH TO PAIRED-END FASTQ(.GZ) DIRECTORY (e.g. /path/to/fq/*{R1,R2}.fq.gz)
OR  --sra_ids   		STRING  SRA SAMPLE ID(s) (Comma separated list of SRA IDS, e.g ['ERR908507', 'ERR908506']
OR  --bamdir    		DIR     PATH TO BAM DIRECTORY (e.g. /path/to/bam/*.bam)
    --genomebase		DIR	PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
    --genome    		STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
    --genomedir			DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
    --genome_name		STRING	REFERENCE GENOME NAME (e.g ".mm10", required if your reference genome is not present in your config file)
    --genome_fasta  		FILE	PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
    --fai			FILE	PATH TO GENOME FAI INDEX FILE (required if your reference genome is not present in your config file)
    --genome2screen 		STRING	GENOMES TO SCREEN FOR FASTQC SCREENING (default : ['mm10','hg19','dm3','dm6','hg38','sacCer2','sacCer3'], comma separated list of genomes to screen reads for contamination, names must correspond to existing genomes in your config file)
    --trimmomatic_adapters  	FILE	PATH TO ADAPTERS FILE FOR TRIMMOMATIC (default ${baseDir}/TruSeq2-PE.fa, special formatting see http://www.usadellab.org/cms/?page=trimmomatic)
                                                                  
Output and Tempory directory parameters:                            
    --name      		STRING    RUN NAME (default : "SSDS_pipeline")      
    --outdir    		DIR       PATH TO OUTPUT DIRECTORY (default : name.outdir)           
    --scratch   		DIR       PATH TO TEMPORARY DIRECTORY (default : scratch)

Pipeline dependencies:
    --src	        	DIR	PATH TO SOURCE DIRECTORY (default : accessoryFiles/SSDS/scripts ; contains perl scripts)
    --custom_bwa        	EXE	PATH TO CUSTOM BWA EXEC (default : accessoryFiles/SSDS/bwa_0.7.12)
    --custom_bwa_ra		EXE	PATH TO CUSTOM BWA_SRA EXEC (default : accessoryFiles/SSDS/bwa_ra_0.7.12)
    --custom_multiqc		EXE	PATH TO CUSTOM MULTIQC EXEC (default : accessoryFiles/SSDS/MultiQC_SSDS_Rev1/bin/multiqc)
    --hotspots	        	DIR	PATH TO HOTSPOTS FILES DIRECTORY (default : accessoryFiles/SSDS/hotspots)
    --blacklist                 FILE    PATH TO BLACKLIST BED FILE FOR PEAK CALLING (default : accessoryFiles/SSDS/blacklist/mm10/blackList.bed)
    --NCIS_dir                  DIR     PATH TO NCIS DIRECTORY (default : accessoryFiles/SSDS/NCIS)

QC parameters
    --with_ssds_multiqc		BOOL	RUN SSDS MULTIQC (need a functional conda environment, see multiqc-dev_conda-env parameter ; default : false)
    --multiqc_dev_conda_env     DIR	PATH TO MULTIQC-DEV CONDA ENVIRONMENT (used when --with_ssds-multiqc is true ; default : multiqc_dev)

Trimming parameters
    --with_trimgalore		BOOL	If you want to trim with trim-galore (default : false)
    --trimgalore_adapters	FILE	OPTIONAL : PATH TO ADAPTERS FILE FOR TRIMGALORE (default : none)
    --trimg_quality		INT	trim-galore : minimum quality (default 10)
    --trimg_stringency		INT	trim-galore : trimming stringency (default 6)
    --trim_minlen		INT	trimmomatic : minimum length of reads after trimming (default 25)
    --trim_crop         	INT	trimmomatic : Cut the read to that specified length (default 50, set to initial length of reads if you want a different crop length for R1 and R2)
    --trim_cropR1		INT	fastx : Cut the R1 read to that specified length (default 50)
    --trim_cropR2		INT	fastx : Cut the R2 read to that specified length (default 50)
    --trim_slidingwin		STRING	trimmomatic : perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default "4:15")
    --trim_illumina_clip	STRING	trimmomatic :  Cut adapter and other illumina-specific sequences from the read (default "2:20:10")

Bam processing parameters
    --bamPGline			STRING	bam header (default '@PG\\tID:ssDNAPipeline1.8_nxf_KBRICK')

Peak calling parameters
    --sctype                    STRING  Saturation curve type (either 'minimal', 'standard' or 'expanded' ; default : 'standard')
    --reps                      INT     Number of repetitions (default : 3)
    --macs_bw                   INT     MACS2 bandwidth paramter (default : 1000)
    --macs_slocal               INT     MACS2 slocal parameter (default : 5000)
    --shuffle_percent           INT     Default : 50                                                                   
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

// Scripts used in the pipeline
def ITR_id_v2c_NextFlow2_script = "${params.src}/ITR_id_v2c_NextFlow2.pl" //Author Kevin Brick
def ssDNA_to_bigwigs_FASTER_LOMEM_script = "${params.src}/ssDNA_to_bigwigs_FASTER_LOMEM.pl" //Author Kevin Brick
def makeSSMultiQCReport_nextFlow_script = "${params.src}/makeSSMultiQCReport_nextFlow.pl" //Author Kevin Brick
def check_design_script = "${params.src}/check_design.py" // From nf-core
def runNCIS_script = "${params.src}/runNCIS.R" //Author Kevin Brick
def pickNlines_script = "${params.src}/pickNlines.pl" //Author Kevin Brick
def satCurveHS_script = "${params.src}/satCurveHS.R" //Author Kevin Brick
def norm_script = "${params.src}/normalizeStrengthByAdjacentRegions.pl" //Author Kevin Brick
def reverse_script = "${params.src}/reverseStrandsForOriCalling.pl" // MISSING // Author Kevin Brick 

// Check if input file exists
if (params.inputcsv) { input_ch = file(params.inputcsv, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genome.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Genome variables
params.genome_fasta = params.genome ? params.genomes[ params.genome ].genome_fasta ?: false : false
params.genomedir = params.genome ? params.genomes[ params.genome ].genomedir ?: false : false
params.genome_name = params.genome ? params.genomes[ params.genome ].genome_name ?: false : false
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false

// Set saturation curve thresholds
if (params.satcurve){
  if (params.sctype == 'expanded'){
    satCurvePCs  = Channel.from(0.025,0.05,0.075,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00)
  }
  if (params.sctype == 'standard'){
    satCurvePCs  = Channel.from(0.20,0.40,0.60,0.80,1.00)
  }
  if (params.sctype == 'minimal'){
    satCurvePCs  = Channel.from(0.10,0.50,1.00)
  }
  useSatCurve  = true
  satCurveReps = params.reps-1
  }
else{
    satCurvePCs  = Channel.from(1.00)
    useSatCurve  = false
    satCurveReps = 0
}


// ******************* //
// BEGINNING PIPELINE  //
// ******************* //
// Get input files according to the input type : sra ; bam or fastq
/*
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
	label 'process_low'
        publishDir "${params.outdir}/preprocess", mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
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
*/
//CHECK INPUT DESIGN FILE
process CHECK_DESIGN {
    tag "${design}"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
        path design from input_ch 
    output:
        path 'design_reads.csv' into ch_design_reads_csv
        path 'design_controls.csv' into ch_design_controls_csv
    script:
    """
    ${check_design_script} $design design_reads.csv design_controls.csv
    """
}
    
//CREATE INPUT CHANNEL WITH FASTQ FILES
ch_design_reads_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
    .set { fq_ch }
    //.println()

ch_design_controls_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(),row.multipleGroups.toBoolean() ] }
    .set { ch_design_controls_csv }
    //.println()

// MAKE CONFIGURATION FILE FOR FASTQSCREEN
process makeScreenConfigFile {
    tag "${outNameStem}" 
    label 'process_basic'
    output:
        file "checkfile.ok" into fqscreen_conf_ok
    script:
        def glist=params.genome2screen
        File conf  = new File("${params.outdir}/conf.fqscreen")
        conf.write "This is a config file for FastQ Screen\n\n"
        conf << "THREADS ${task.cpus}\n\n"
        for (item in glist) {
	    if (params.genomes.containsKey(item)) {
            	fasta=params.genomes[ item ].genome_fasta
            	name=params.genomes[ item ].genome_name
            	conf << "DATABASE  ${name}    ${fasta}\n"
	    }
        }
	"""
	if [ -f "${params.outdir}/conf.fqscreen" ]; then
    		echo "${params.outdir}/conf.fqscreen exists." > checkfile.ok
                mkdir -p "${params.outdir}/fastqscreen"
                mv "${params.outdir}/conf.fqscreen" "${params.outdir}/fastqscreen/"
	else
		exit 1, "The configuration file for fastqscreen could not be generated. Please check your genome2screen parameter."
	fi
	"""
}

// TRIMMING : USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES
process trimming {
    tag "${sampleId}" 
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*_report.txt"
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*.zip"
    input:
        set val(sampleId), file(reads) from fq_ch
    output:
        set val("${sampleId}"), file(reads) into fqc_ch 
        set val("${sampleId}"), file('*R1.fastq.gz'), file('*R2.fastq.gz') into trim_ch
        val 'ok_multiqc' into trimming_ok
    script:
    	if (params.with_trimgalore && params.trimgalore_adapters)
	"""
	trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} --cores ${task.cpus} --adapter "file:${params.trimgalore_adapters}" --gzip --paired --basename ${sampleId} ${reads}
        mv ${sampleId}_val_1.fq.gz ${sampleId}_trim_R1.fastq.gz
        mv ${sampleId}_val_2.fq.gz ${sampleId}_trim_R2.fastq.gz
        fastqc -t ${task.cpus} ${sampleId}_trim_R1.fastq.gz ${sampleId}_trim_R2.fastq.gz
        """
        else if (params.with_trimgalore && !params.trimgalore_adapters)
        """
        trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} --cores ${task.cpus} --gzip --paired --basename ${sampleId} ${reads}
	mv ${sampleId}_val_1.fq.gz ${sampleId}_trim_R1.fastq.gz
	mv ${sampleId}_val_2.fq.gz ${sampleId}_trim_R2.fastq.gz
	fastqc -t ${task.cpus} ${sampleId}_trim_R1.fastq.gz ${sampleId}_trim_R2.fastq.gz
	"""
	else
	"""
    	trimmomatic PE -threads ${task.cpus} ${reads} \
                ${sampleId}_trim_R1.fastq.gz R1_unpaired.fastq.gz \
                ${sampleId}_trim_R2.fastq.gz R2_unpaired.fastq.gz \
                ILLUMINACLIP:${params.trimmomatic_adapters}:${params.trim_illuminaclip} SLIDINGWINDOW:${params.trim_slidingwin} MINLEN:${params.trim_minlen} CROP:${params.trim_crop} \
                >& ${sampleId}_trim_${outNameStem}_trimmomatic_report.txt 2>&1
    	fastqc -t ${task.cpus} ${sampleId}_trim_R1.fastq.gz ${sampleId}_trim_R2.fastq.gz 
    	"""
}

// QUALITY CONTROL : RUN FASTQC AND FASTQSCREEN
process fastqc {
    tag "${sampleId}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*.zip"
    publishDir "${params.outdir}/fastqscreen", mode: 'copy', overwrite: false, pattern: "*.png"
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: false, pattern: "*.txt"
    input:
        set val(sampleId), file(reads) from fqc_ch
        file(ok) from fqscreen_conf_ok
    output:
	file('*') into fastc_report
        val 'ok' into fastqc_ok
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    fastq_screen --threads ${task.cpus} --force --aligner bwa --conf ${params.outdir}/fastqscreen/conf.fqscreen ${reads}
    """
}

// MAPPING : USE CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA
process bwaAlign {
    tag "${sampleId}"
    label 'process_long'
    input:
        set val(sampleId), file(fqR1), file(fqR2) from trim_ch
    output:
        set val(sampleId), file('*.sorted.bam') into sortedbam2filterbam, sortedbam2parseITR
        val 'ok' into bwa_ok
    script:
  """
    if [ ${params.trim_cropR1} != ${params.trim_cropR2} ]
    then
        zcat ${fqR1} | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${tmpNameStem}.R1.fastq.gz
        zcat ${fqR2} | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${tmpNameStem}.R2.fastq.gz
    else
        mv ${fqR1} ${tmpNameStem}.R1.fastq.gz
        mv ${fqR2} ${tmpNameStem}.R2.fastq.gz
    fi
    ${params.custom_bwa} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R1.fastq.gz \
            > ${tmpNameStem}.R1.sai
    ${params.custom_bwa_ra} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R2.fastq.gz \
            > ${tmpNameStem}.R2.sai
    ${params.custom_bwa} sampe ${params.genome_fasta} ${tmpNameStem}.R1.sai ${tmpNameStem}.R2.sai ${tmpNameStem}.R1.fastq.gz \
            ${tmpNameStem}.R2.fastq.gz >${tmpNameStem}.unsorted.sam
    picard SamFormatConverter I=${tmpNameStem}.unsorted.sam O=${tmpNameStem}.unsorted.tmpbam VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${tmpNameStem}.unsorted.tmpbam O=${sampleId}.sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${sampleId}.sorted.bam 
    
    """
}

// BAM FILTERING
process filterBam {
    tag "${sampleId}"
    label 'process_medium'
    publishDir "${params.outdir}/filterbam",  mode: 'copy', overwrite: false, pattern: "*.txt"
    input:
        set val(sampleId), file(bam) from sortedbam2filterbam
    output:
        set val(sampleId), file('*.unparsed.bam') into bamAligned
        set val(sampleId), file('*.unparsed.bam.bai') into bamIDX 
        set val(sampleId), file('*.unparsed.suppAlignments.bam') into bamAlignedSupp
        set val(sampleId), file('*.unparsed.suppAlignments.bam.bai') into bamIDXSupp
        val 'ok' into filterbam_ok
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

// PARSE ITRS
process parseITRs {
    tag "${sampleId}"
    label 'process_medium'
    publishDir "${params.outdir}/parseitr",  mode: 'copy', overwrite: false, pattern: "*.txt"
    publishDir "${params.outdir}/parseitr",  mode: 'copy', overwrite: false, pattern: "*.bed"
    input:
        set val(sampleId), file(bam) from sortedbam2parseITR 
        //file(bam) from bamAligned
        //file(bai) from bamIDX
    output:
        set val(sampleId), file("${bam}"), file('*bed') into ITRBED
        set val(sampleId), file('*ssDNA_type1.bed') into T1BED
        set val(sampleId), file('*dsDNA.bed') into DSBED
        set val(sampleId), file('*.md.*bam'),file('*.md.*bam.bai') into BAMwithIDXfr, BAMwithIDXss, BAMwithIDXdt mode flatten
        val 'ok' into parseitr_ok
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


// MAKE DEEPTOOLS BIGWIG
process makeDeeptoolsBigWig { 
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: false, pattern: "*.png"
    publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: false, pattern: "*.bigwig"
    publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: false, pattern: "*.tab"
    input:
        set val(sampleId), file(bam), file(bamidx) from BAMwithIDXdt
    output:
        file('*') into bigwig
        val 'ok' into bigwig_ok
    shell:
    """
    bamCoverage --bam ${bam} --normalizeUsing RPKM --binSize ${params.binsize} --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.RPKM.bigwig 
    plotCoverage --bamfiles ${bam} --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.coveragePlot.png
    plotFingerprint --bamfiles ${bam} --labels ${bam} --numberOfProcessors ${task.cpus} \
        --minMappingQuality 30 --skipZeros --plotFile ${bam.baseName}.deeptools.fingerprints.png \
        --outRawCounts ${bam.baseName}.deeptools.fingerprints.tab
    """
}

// SAM STATS
process samStats {
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/samstats",  mode: 'copy', overwrite: false, pattern: "*.tab"
    input:
        set val(sampleId), file(bam), file(bamidx) from BAMwithIDXss
    output:
        file '*stats.tab' into dtSamStat
        val 'ok' into samstats_ok
    script:
        """
        samtools idxstats ${bam} > ${bam.baseName}.idxstats.tab
        samtools stats ${bam} > ${bam.baseName}.samstats.tab
        """
}

// FWD/REV bigwig 
process toFRBigWig {
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: false
    input:
        set val(sampleId), file(bam), file(bamidx) from BAMwithIDXfr
    output:
        file('*') into frbigwig
	val 'ok' into frbigwig_ok
    script:
        """
        perl ${ssDNA_to_bigwigs_FASTER_LOMEM_script} --bam ${bam} --g ${params.genome} --o ${bam.baseName}.out --s 100 --w 1000 --sc ${params.scratch} --gIdx ${params.fai} -v
        """
}    

// SSDS REPORT
process makeSSreport {
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/multiqc",  mode: 'copy', overwrite: false
    input:
        set val(sampleId), file(bam), file(ssdsBEDs) from ITRBED 
    output:
        set val(sampleId), file('*SSDSreport*') into SSDSreport2ssdsmultiqc
        val 'ok' into ssreport_ok
    script:
        """
        perl ${makeSSMultiQCReport_nextFlow_script} ${bam} ${ssdsBEDs} --g ${params.genome} --h ${params.hotspots}
        """
}

if (params.with_ssds_multiqc) {
    // SSDS MULTIQC
    process ssds_multiqc {
        tag "${sampleId}"
    	label 'process_basic'
    	conda "${params.multiqc_dev_conda_env}" 
        publishDir "${params.outdir}/multiqc",  mode: 'copy', overwrite: false
        input:
            set val(sampleId), file(report) from SSDSreport2ssdsmultiqc
        output:
            file('*') into ssdsmultiqc_report
        script:
        """
        multiqc -m ssds -n ${sampleId}.multiQC .
        """
    }
}

if (params.with_control) {
    // CREATE CHANNEL LINKING IP T1SSDNA BED WITH CONTROL DSDNA BED
    T1BED
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { T1BED }

    DSBED
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { DSBED }

    ch_design_controls_csv
        .combine(T1BED)
        .combine(DSBED)
        .filter { it[0] == it[5] && it[1] == it[7] }
        .map { it -> it[0,1,6,8].flatten() } 
        .set { T1BED_shuffle_ch }

    //BED SHUFFLING
    process shufBEDs_ct {
        tag "${id_ip}"
        label 'process_basic'
        publishDir "${params.outdir}/bed_shuffle",  mode: 'copy', overwrite: false
        input:
            tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed) from T1BED_shuffle_ch
        output:
            tuple val(id_ip), val(id_ct), path("*.IP.sq30.bed"), path("*.CT.sq30.bed") into SQ30BED_ch
            val 'ok' into shufbed_ok
        script:
            def ip_Q30_bed      = ip_bed.name.replaceFirst(".bed",".IP.q30.bed")
            def ip_Q30_shuf_bed = ip_bed.name.replaceFirst(".bed",".IP.sq30.bed")

            def ct_Q30_bed        = ct_bed.name.replaceFirst(".bed",".CT.q30.bed")
            def ct_Q30_shuf_bed   = ct_bed.name.replaceFirst(".bed",".CT.sq30.bed")
            """
            perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' ${ip_bed} >${ip_Q30_bed}
            perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' ${ct_bed} >${ct_Q30_bed}
            shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
            shuf ${ct_Q30_bed}   |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ct_Q30_shuf_bed}
            """
    }

    //PEAK CALLING WITH MACS2
    process callPeaks_ct {
        tag "${id_ip}"
        label 'process_basic'
        conda "${baseDir}/environment_callpeaks3.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', overwrite: true, pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', overwrite: true, pattern: "*peaks.be*"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', overwrite: true, pattern: "*peaks.xls"
        publishDir "${params.outdir}/model",                  mode: 'copy', overwrite: true, pattern: "*model*"
        input:
            tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed) from SQ30BED_ch
        output:
            path("*peaks_sc.bed") into allbed
            path("*peaks.xls") optional true into peaks_xls
            path("*peaks.bed") optional true into peaks_bed
            path("*peaks.bedgraph") optional true into peaks_bg
            path("*model.R") optional true into model_R
            path("*model.pdf") optional true into model_pdf
            val 'ok' into callPeaks_ok
        script:
        """
        nT=`cat ${ip_bed} |wc -l`
        nPC=`perl -e 'print int('\$nT'*${params.shuffle_percent})'`

        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## NCIS : Normalization for ChIp-Seq (just use chrom1: faster with analagous results)
        grep -w chr1 \$nPC.IP.bed >IP.cs1.bed
        grep -w chr1 \$ct_bed >CT.cs1.bed
        #Rscript ${runNCIS_script} IP.cs1.bed CT.cs1.bed ${params.NCIS_dir} NCIS.out
        #ratio=`cut -f1 NCIS.out`
 
        ## GET GENOME SIZE - BLACKLIST SIZE
        tot_sz=`cut -f3 ${params.fai} |tail -n1`
        bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
        genome_size=`expr \$tot_sz - \$bl_size`
    
        for i in {0..${satCurveReps}}; do
            thisName=${id_ip}'.N'\$nPC'_${params.shuffle_percent}pc.'\$i
            perl ${pickNlines_script} ${ip_bed} \$nPC >\$nPC.tmp
    
            sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed

            #macs2 callpeak  --ratio \$ratio \\ 
	    if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak \\
                    -g \$genome_size \\
                    -t \$nPC.IP.bed \\
                    -c ${ct_bed} \\
                    --bw ${params.macs_bw} \\
                    --keep-dup all \\
                    --slocal ${params.macs_slocal} \\
                    --name \$thisName \\
	            --nomodel \\
                    --extsize ${params.macs_extsize} \\
                    --pvalue ${params.macs_pv}
            else
                    macs2 callpeak \\
                    -g \$genome_size \\
                    -t \$nPC.IP.bed \\
                    -c ${ct_bed} \\
                    --bw ${params.macs_bw} \\
                    --keep-dup all \\
                    --slocal ${params.macs_slocal} \\
                    --name \$thisName \\
                    --nomodel \\
                    --extsize ${params.macs_extsize} \\
                    --qvalue ${params.macs_qv}
            fi

            intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
        
            cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done
  
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${params.shuffle_percent}.peaks_sc.bed
  
        if [ ${params.shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${params.shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
            cat *1.00*.r >${id_ip}.model.R
            R --vanilla <${id_ip}.model.R
            cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
            ## Calculate strength
            perl ${norm_script} --bed ${id_ip}.peaks.bed \
                 --in ${ip_bed} --out ${id_ip}.peaks.bedgraph --rc --rev_src ${reverse_script}
        fi
        """
    }

} else {

    //BED SHUFFLING
    process shufBEDs {
        tag "${id_ip}"
        label 'process_basic'
        publishDir "${params.outdir}/bed_shuffle",  mode: 'copy', overwrite: false
        input:
            tuple val(id_ip), path(ip_bed) from T1BED
        output:
            tuple val(id_ip), path("*.IP.sq30.bed") into SQ30BED_ch
            val 'ok' into shufbed_ok
        script:
            def ip_Q30_bed      = ip_bed.name.replaceFirst(".bed",".IP.q30.bed")
            def ip_Q30_shuf_bed = ip_bed.name.replaceFirst(".bed",".IP.sq30.bed")
            """
            perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' ${ip_bed} >${ip_Q30_bed}
            shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
            """
    }

    //PEAK CALLING WITH MACS2
    process callPeaks {
        tag "${id_ip}"
        label 'process_basic'
        conda "${baseDir}/environment_callpeaks3.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', overwrite: true, pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', overwrite: true, pattern: "*peaks.be*"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', overwrite: true, pattern: "*peaks.xls"
        publishDir "${params.outdir}/model",                  mode: 'copy', overwrite: true, pattern: "*model*"
        input:
            tuple val(id_ip), path(ip_bed) from SQ30BED_ch
        output:
            path("*peaks_sc.bed") into allbed
            path("*peaks.xls") optional true into peaks_xls
            path("*peaks.bed") optional true into peaks_bed
            path("*peaks.bedgraph") optional true into peaks_bg
            path("*model.R") optional true into model_R
            path("*model.pdf") optional true into model_pdf
            val 'ok' into callPeaks_ok
        script:
        """
        nT=`cat ${ip_bed} |wc -l`
        nPC=`perl -e 'print int('\$nT'*${params.shuffle_percent})'`

        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## GET GENOME SIZE - BLACKLIST SIZE
        tot_sz=`cut -f3 ${params.fai} |tail -n1`
        bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
        genome_size=`expr \$tot_sz - \$bl_size`
    
        for i in {0..${satCurveReps}}; do
            thisName=${id_ip}'.N'\$nPC'_${params.shuffle_percent}pc.'\$i
            perl ${pickNlines_script} ${ip_bed} \$nPC >\$nPC.tmp
    
            sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
            if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak \\
                    -g \$genome_size \\
                    -t \$nPC.IP.bed \\
                    --bw ${params.macs_bw} \\
                    --keep-dup all \\
                    --slocal ${params.macs_slocal} \\
                    --name \$thisName \\
	            --nomodel \\
                    --extsize ${params.macs_extsize} \\
                    --pvalue ${params.macs_pv}
            else
                macs2 callpeak \\
                -g \$genome_size \\
                -t \$nPC.IP.bed \\
                --bw ${params.macs_bw} \\
                --keep-dup all \\
                --slocal ${params.macs_slocal} \\
                --name \$thisName \\
                --nomodel \\
                --extsize ${params.macs_extsize} \\
                --qvalue ${params.macs_qv}
            fi


            intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
        
            cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done
  
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${params.shuffle_percent}.peaks_sc.bed
  
        if [ ${params.shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${params.shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
            cat *1.00*.r >${id_ip}.model.R
            R --vanilla <${id_ip}.model.R
            cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
            ## Calculate strength
            perl ${norm_script} --bed ${id_ip}.peaks.bed \
                 --in ${ip_bed} --out ${id_ip}.peaks.bedgraph --rc --rev_src ${reverse_script}
        fi
        """
    }
}

process makeSatCurve {
    tag "${id_ip}"
    label 'process_basic'
    conda "${baseDir}/environment_callpeaks3.yml"
    publishDir "${params.outdir}/saturation_curve",  mode: 'copy', overwrite: true
    input:
        file(saturation_curve_data) from allbed.collect()
    output:
        path("*satCurve.tab", emit: table) into satcurve_table
        path("*.png", emit: png) into curve
        val 'ok' into makeSatCurve_ok
    script:
    """
    #!/usr/bin/perl
    use strict;
    my \$tf = "satCurve.tab";
    open TMP, '>', \$tf;
    print TMP join("\\t","reads","pc","hs")."\\n";
    open my \$IN, '-|', 'cd ${params.outdir}/saturation_curve/peaks ; wc -l *peaks_sc.bed |grep -v total |sort -k1n,1n ';
    while (<\$IN>){
  	chomp;
  	next if (\$_ =~ /\\stotal\\s*\$/);
  	\$_ =~ /^\\s*(\\d+).+\\.N(\\d+)_([\\d\\.]+)pc.+\$/;
  	my (\$HS,\$N,\$pc) = (\$1,\$2,\$3*100);
  	print TMP join("\\t",\$N,\$pc,\$HS)."\\n";
    } 
    close TMP;
    close \$IN;
    my \$R = \$tf.'.R';
    makeRScript(\$R,"${params.name}.satCurve.tab","${params.name}");
    system('R --vanilla <'.\$R);
    sub makeRScript{
     	my (\$sName,\$data,\$sampleName) = @_;
  	open RS, '>', \$sName;
  	print RS 'source("${satCurveHS_script}")'."\\n";
  	print RS 'satCurveHS(fIN = "'.\$tf.'", sampleName = "'.\$sampleName.'")'."\\n";
  	close RS;
    }
    """
}

// GENERAL MULTIQC
process general_multiqc {
    tag "${outNameStem}"
    label 'process_basic'
    conda 'bioconda::multiqc=1.9'
    publishDir "${params.outdir}/multiqc",  mode: 'copy', overwrite: false
    input:
        val('trimming_ok') from trimming_ok.collect()
	val('fastqc_ok') from fastqc_ok.collect()
	val('bwa_ok') from bwa_ok.collect()
	val('filterbam') from filterbam_ok.collect()
	val('parseitr') from parseitr_ok.collect()
	val('bigwig') from bigwig_ok.collect()
	val('samstats') from samstats_ok.collect()
	val('frbigwig_ok') from frbigwig_ok.collect()
        val('ssreport_ok') from ssreport_ok.collect()
        val('shufbed_ok') from shufbed_ok.collect()
        val('callPeaks_ok') from callPeaks_ok.collect()
        val('makeSatCurve_ok') from makeSatCurve_ok.collect()
    output:
	file('*') into generalmultiqc_report
    script:
    """
    multiqc -n ${outNameStem}.multiQC ${params.outdir}/
    """
}

// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    scrdir.deleteDir()
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

