#!/usr/bin/env nextflow
/*
========================================================================================
                        SSDS Pipeline version 2.0
			Author : Kevin Brick (version 1.8_NF)
                        Pauline Auffret (version 2.0)
========================================================================================
 SSDS nextflow pipeline
 #### Homepage / Documentation
 https://github.com/kevbrick/SSDSnextflowPipeline
 https://github.com/kevbrick/callSSDSpeaks
 Adapted from version 1.8_NF (Pauline Auffret, 2020)
 https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Single-Stranded-DNA-Sequencing (SSDS) Pipeline : Align, Parse and call Peaks in ssDNA
Pipeline overview:
0. makeScreenConfigFile Generates configuration file for fastqscreen in accordance with params.genome2screen
1. trimming             Runs Trimmomatic or trim_galore for quality trimming, adapters removal and hard trimming
2. fastqc               Runs FastQC for sequencing reads quality control
3. bwaAlign             Runs BWA and Custom BWA algorithm (bwa-ra : bwa right align) to map reads to reference genome
4. filterBam            Uses PicardTools and Samtools for duplicates marking ; removing supplementary alignements ; bam sorting and indexing
5. parseITRs            Uses in-house perl script and samtools to parse BAM files into type1 ss dna, type 2 ss dna, ds dna, strict ds dna, unclassfied (5 types bam files)
6. makeDeeptoolsBigWig  Runs deeptool to make bigwig files for each 5 types bam files
7. samStats             Runs samtools to make stats on 5 types bam files
8. toFRBigWig           Uses in-house perl to make FWD bigwig files
9. makeSSreport         Uses in-house perl script to parse bed files
10. ssds_multiqc        Runs custom multiqc to edit report on SSDS alignement
11. shufBEDS            Perfoms random shuffling in bed files
12. callPeaks           Calls peaks in bed files
13. makeSatCurve        Makes saturation Curve
14. general_multiqc     Runs multiqc to edit a general stat report
----------------------------------------------------------------------------------------
*/
def helpMessage() { 
    log.info"""
=============================================================================
  SSDS Pipeline version 2.0
=============================================================================
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
=================================================================================================
""".stripIndent()
}

params.version = false
if (params.version){
    println("This is $workflow.manifest.name version $workflow.manifest.version")
    println("$workflow.manifest.description")
    exit 0
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
def pickNlines_script = "${params.src}/pickNlines.pl" //Author Kevin Brick
def satCurveHS_script = "${params.src}/satCurveHS.R" //Author Kevin Brick
def norm_script = "${params.src}/normalizeStrengthByAdjacentRegions.pl" //Author Kevin Brick
def reverse_script = "${params.src}/reverseStrandsForOriCalling.pl" // MISSING // Author Kevin Brick 
def getPeaksBedFiles_script = "${params.src}/getPeaksBedFiles.pl" //Author Kevin Brick, script adapted by Pauline Auffret
def runSatCurve_script = "${params.src}/runSatCurve.R" //Author Pauline Auffret
def encode_idr_script= "${params.src}/encode-dcc_chip-seq-pipeline2_src/encode_task_idr.py" //Author Jin Lee from https://github.com/ENCODE-DCC/chip-seq-pipeline2

// Check if input csv file exists
if (params.inputcsv) { input_ch = file(params.inputcsv, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genome.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Define genome variables
params.genome_fasta = params.genome ? params.genomes[ params.genome ].genome_fasta ?: false : false
params.genomedir = params.genome ? params.genomes[ params.genome ].genomedir ?: false : false
params.genome_name = params.genome ? params.genomes[ params.genome ].genome_name ?: false : false
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false

// Set saturation curve thresholds for processes callPeaks and makeSatCurve
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
// CHECK INPUT DESIGN FILE
// This process controls the input csv file (checks if there is the right columns, and valid file extension).
// It outputs 2 csv files for mapping fastq files to their sample ID and to map chIP files to corresponding control files if needed.
// This process uses the python script ${check_design_script} adapted from nf-core chipseq pipeline.
process CHECK_DESIGN {
    tag "${design}"
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
        path design from input_ch 
    output:
        path 'design_reads.csv' into ch_design_reads_csv
        path 'design_controls.csv' into ch_design_controls_csv
    script:
    """
    python ${check_design_script} $design design_reads.csv design_controls.csv
    # Check if the --with_control parameter value is consistent with csv input file
    nb_line_ctrl_file=`cat design_controls.csv | wc -l`
    if [[ \$nb_line_ctrl_file == "1" && ${params.with_control} == "true" ]]
    then
        echo "The option --with_control is set to true but there is no control files. Check the input csv file."
        exit 1
    fi
    """
}
    
//CREATE INPUT CHANNEL WITH SAMPLE ID AND FASTQ FILES
ch_design_reads_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
    .set { fq_ch }
    //.println()

//CREATE INPUT CHANNEL MAPPING chIP SAMPLE ID AND control SAMPLE ID
ch_design_controls_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(),row.multipleGroups.toBoolean() ] }
    .set { ch_design_controls_csv }
    //.println()

// MAKE CONFIGURATION FILE FOR FASTQSCREEN
// This process generates a configuration file for fastqscreen which will contain the list of genomes to be screened during general QC.
// The list of genomes is defined in the parameter --genome2screen
// The output is a text file named conf.fqscreen
process makeScreenConfigFile {
    tag "${outNameStem}" 
    label 'process_basic'
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
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
                echo "The configuration file for fastqscreen could not be generated. Please check your genome2screen parameter."
		exit 1
	fi
	"""
}

// TRIMMING : USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES
// This process runs trimmomatic (default) or Trim Galore (if option --with_trimagalore is set) for adapter & quality trimming.
// Several parameters can be set for the trimming, see help section.
// For trimmomatic an adapter file need to be set, and for trim galore the choice is yours.
// Finally hard trimming is done using fastx-trimmer on trimmed fastq files and QC is done with fastqc.
process trimming {
    tag "${sampleId}" 
    label 'process_low'
    publishDir "${params.outdir}/trim_fastqc", mode: 'copy', pattern: "*_report.txt"
    publishDir "${params.outdir}/trim_fastqc", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/trim_fastqc", mode: 'copy', pattern: "*.zip"
    publishDir "${params.outdir}/trim_fastq", mode: 'copy', pattern: "*trim_crop_R1.fastq.gz"
    publishDir "${params.outdir}/trim_fastq", mode: 'copy', pattern: "*trim_crop_R2.fastq.gz"
    input:
        set val(sampleId), file(reads) from fq_ch
    output:
        set val("${sampleId}"), file(reads) into fqc_ch
        set val("${sampleId}"), file('*crop_R1.fastq.gz'), file('*crop_R2.fastq.gz') into trim_ch
        file('*') into trim_fastqc_report
        val 'ok_multiqc' into trimming_ok
    script:
    	if (params.with_trimgalore && params.trimgalore_adapters)
	"""
	trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
                --cores ${task.cpus} --adapter "file:${params.trimgalore_adapters}" --gzip --paired --basename ${sampleId} ${reads}
        mv ${sampleId}_val_1.fq.gz ${sampleId}_trim_R1.fastq.gz
        mv ${sampleId}_val_2.fq.gz ${sampleId}_trim_R2.fastq.gz
        zcat ${sampleId}_trim_R1.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${sampleId}_trim_crop_R1.fastq.gz
        zcat ${sampleId}_trim_R2.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${sampleId}_trim_crop_R2.fastq.gz
        fastqc -t ${task.cpus} ${sampleId}_trim_crop_R1.fastq.gz ${sampleId}_trim_crop_R2.fastq.gz
        """
        else if (params.with_trimgalore && !params.trimgalore_adapters)
        """
        trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
                --cores ${task.cpus} --gzip --paired --basename ${sampleId} ${reads}
	mv ${sampleId}_val_1.fq.gz ${sampleId}_trim_R1.fastq.gz
	mv ${sampleId}_val_2.fq.gz ${sampleId}_trim_R2.fastq.gz
        zcat ${sampleId}_trim_R1.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${sampleId}_trim_crop_R1.fastq.gz
        zcat ${sampleId}_trim_R2.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${sampleId}_trim_crop_R2.fastq.gz
        fastqc -t ${task.cpus} ${sampleId}_trim_crop_R1.fastq.gz ${sampleId}_trim_crop_R2.fastq.gz
	"""
	else
	"""
    	trimmomatic PE -threads ${task.cpus} ${reads} \
                ${sampleId}_trim_R1.fastq.gz R1_unpaired.fastq.gz \
                ${sampleId}_trim_R2.fastq.gz R2_unpaired.fastq.gz \
                ILLUMINACLIP:${params.trimmomatic_adapters}:${params.trim_illuminaclip} SLIDINGWINDOW:${params.trim_slidingwin} \
                MINLEN:${params.trim_minlen} >& ${sampleId}_trim_${outNameStem}_trimmomatic_report.txt 2>&1
        zcat ${sampleId}_trim_R1.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${sampleId}_trim_crop_R1.fastq.gz
        zcat ${sampleId}_trim_R2.fastq.gz | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${sampleId}_trim_crop_R2.fastq.gz
        fastqc -t ${task.cpus} ${sampleId}_trim_crop_R1.fastq.gz ${sampleId}_trim_crop_R2.fastq.gz
        """
}

// MAP FASTQC CHANNEL
fqc_ch
    .map { it -> it[0,1].flatten() }
    .set { fqc_tuple }

// QUALITY CONTROL : RUN FASTQC AND FASTQSCREEN
process fastqc {
    tag "${sampleId}"
    label 'process_low'
    publishDir "${params.outdir}/raw_fastqc", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/raw_fastqc", mode: 'copy', pattern: "*.zip"
    publishDir "${params.outdir}/fastqscreen", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/raw_fastqc", mode: 'copy', pattern: "*.txt"
    input:
        tuple val(sampleId), file(read1), file(read2) from fqc_tuple 
        file(ok) from fqscreen_conf_ok
    output:
	file('*') into fastc_report
        val 'ok' into fastqc_ok
    script:
    ext=("${read1.getExtension()}")
    """
    if [[ ${ext} == "gz" ]] ; then ext="fastq.gz" ; fi
    mv ${read1} ${sampleId}_raw_R1.${ext}   
    mv ${read2} ${sampleId}_raw_R2.${ext}
    fastqc -t ${task.cpus} *raw* 
    fastq_screen --threads ${task.cpus} --force --aligner bwa --conf ${params.outdir}/fastqscreen/conf.fqscreen *raw*
    """
}

// MAPPING : USE BWA aAND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA
process bwaAlign {
    tag "${sampleId}"
    label 'process_long'
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    publishDir "${params.outdir}/bam",  mode: 'copy', pattern: "*.sorted.bam*"
    publishDir "${params.outdir}/bam",  mode: 'copy', pattern: "*.flagstat"
    input:
        set val(sampleId), file(fqR1), file(fqR2) from trim_ch
    output:
        set val(sampleId), file('*.sorted.bam') into sortedbam2filterbam, sortedbam2parseITR
        set file('*.sorted.bam'), file('*.sorted.bam.bai') into sortedAndIndexedBam
        file('*.flagstat') into flagstat_report_bwa
        val 'ok' into bwa_ok
    script:
    """
    ${params.custom_bwa} aln -t ${task.cpus} ${params.genome_fasta} ${fqR1} \
            > ${tmpNameStem}.R1.sai
    ${params.custom_bwa_ra} aln -t ${task.cpus} ${params.genome_fasta} ${fqR2} \
            > ${tmpNameStem}.R2.sai
    ${params.custom_bwa} sampe ${params.genome_fasta} ${tmpNameStem}.R1.sai ${tmpNameStem}.R2.sai ${fqR1} \
            ${fqR2} >${tmpNameStem}.unsorted.sam
    picard SamFormatConverter I=${tmpNameStem}.unsorted.sam O=${tmpNameStem}.unsorted.tmpbam \
            VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    picard SortSam I=${tmpNameStem}.unsorted.tmpbam O=${sampleId}.sorted.bam SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${sampleId}.sorted.bam 
    
    ## CHECK IF BAM FILE IS EMPTY AFTER MAPPING
    samtools flagstat ${sampleId}.sorted.bam > ${sampleId}.sorted.bam.flagstat
    if grep -q "0 + 0 mapped" ${sampleId}.sorted.bam.flagstat
    then
        echo "Bam file ${sampleId}.sorted.bam is empty. Check genome parameter."
        exit 1
    fi

    ## REMOVE MULTIMAPPERS IF OPTION --no_multimap IS TRUE
    ## ! It's important that the final bam files are named *.sorted.bam for the next processes.
    if [[ ${params.no_multimap} ]]
    then
        samtools view -h ${sampleId}.sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${sampleId}.final.bam
        samtools index ${sampleId}.final.bam
        rm ${sampleId}.sorted.bam ; mv ${sampleId}.final.bam ${sampleId}.sorted.bam
        rm ${sampleId}.sorted.bam.bai ; mv ${sampleId}.final.bam.bai ${sampleId}.sorted.bam.bai
    fi
    """
}

// BAM FILTERING
process filterBam {
    tag "${sampleId}"
    label 'process_medium'
    publishDir "${params.outdir}/filterbam",  mode: 'copy', pattern: "*.txt"
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
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.txt"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.bed"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.md.ssDNA_type1.bam*"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.md.ssDNA_type2.bam*"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.md.dsDNA.bam*"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.md.dsDNA_strict.bam*"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.md.unclassified.bam*"
    publishDir "${params.outdir}/parse_itr",  mode: 'copy', pattern: "*.flagstat"
    input:
        set val(sampleId), file(bam) from sortedbam2parseITR 
        //file(bam) from bamAligned
        //file(bai) from bamIDX
    output:
        set val(sampleId), file("${bam}"), file('*bed') into ITRBED
        set val(sampleId), file('*ssDNA_type1.bed') into T1BED, T1BEDrep
        set val(sampleId), file('*dsDNA.bed') into DSBED, DSBEDrep
        set val(sampleId), file('*.md.*bam'),file('*.md.*bam.bai') into BAMwithIDXfr, BAMwithIDXss, BAMwithIDXdt mode flatten
        file ('*.flagstat') into flagstat_report_parseITRs
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

    ## CHECK IF TYPE 1 BAM FILE IS EMPTY AFTER PARSING
    samtools flagstat ${bam}.md.ssDNA_type1.bam > ${bam}.md.ssDNA_type1.bam.flagstat
    if grep -q "0 + 0 mapped" ${bam}.md.ssDNA_type1.bam.flagstat
    then
        echo "Type 1 bam file ${bam}.md.ssDNA_type1.bam is empty. Check genome parameter."
        exit 1
    fi
    """
}


// MAKE DEEPTOOLS BIGWIG
process makeDeeptoolsBigWig { 
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/bigwig",  mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/bigwig",  mode: 'copy', pattern: "*.bigwig"
    publishDir "${params.outdir}/bigwig",  mode: 'copy', pattern: "*.tab"
    input:
        set val(sampleId), file(bam), file(bamidx) from BAMwithIDXdt
    output:
        file('*') into bigwig
        val 'ok' into bigwig_ok
    shell:
    """
    bamCoverage --bam ${bam} --normalizeUsing RPKM --binSize ${params.binsize} \
        --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.RPKM.bigwig 
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
    publishDir "${params.outdir}/samstats",  mode: 'copy', pattern: "*.tab"
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
    publishDir "${params.outdir}/bigwig",  mode: 'copy'
    input:
        set val(sampleId), file(bam), file(bamidx) from BAMwithIDXfr
    output:
        file('*') into frbigwig
	val 'ok' into frbigwig_ok
    script:
        """
        perl ${ssDNA_to_bigwigs_FASTER_LOMEM_script} --bam ${bam} --g ${params.genome} --o ${bam.baseName}.out \
            --s 100 --w 1000 --sc ${params.scratch} --gIdx ${params.fai} -v
        """
}    

// SSDS REPORT
process makeSSreport {
    tag "${sampleId}"
    label 'process_basic'
    publishDir "${params.outdir}/multiqc",  mode: 'copy'
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
        publishDir "${params.outdir}/multiqc",  mode: 'copy'
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
    // CREATE CHANNEL LINKING IP TYPE 1 SSDNA BED WITH CONTROL DSDNA BED
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
        .into { T1BED_shuffle_ch ; T1BED_replicate_ch }
        

    //BED SHUFFLING
    process shufBEDs_ct {
        tag "${id_ip}"
        label 'process_basic'
        publishDir "${params.outdir}/bed_shuffle",  mode: 'copy'
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
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks.be*"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks.xls"
        input:
            tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed) from SQ30BED_ch
            val(shuffle_percent) from satCurvePCs
        output:
            path("*peaks_sc.bed") into allbed
            path("*peaks.xls") optional true into peaks_xls
            path("*peaks.bed") optional true into peaks_bed
            path("*peaks.bedgraph") optional true into peaks_bg
            stdout into gsize
            val 'ok' into callPeaks_ok
        script:
        """
        ## SELECT N LINES FROM IP BED FILE ACCORDING TO satCurve parameter 
        nT=`cat ${ip_bed} |wc -l`
        nPC=`perl -e 'print int('\$nT'*${shuffle_percent})'`
        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## GET GENOME SIZE - BLACKLIST SIZE
        tot_sz=`cut -f3 ${params.fai} |tail -n1`
        bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
        genome_size=`expr \$tot_sz - \$bl_size`
        echo -n \$genome_size
        
        ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter 
        for i in {0..${satCurveReps}}; do
            thisName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
            perl ${pickNlines_script} ${ip_bed} \$nPC >\$nPC.tmp
    
            sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed

	    if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} >> ${params.scratch}/macs2.log
            else
                macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize}  --qvalue ${params.macs_qv} >> ${params.scratch}/macs2.log
            fi

            intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
        
            cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done
  
        ##MERGE BED FILES
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
        
        ##POSTPROCESS IF params.satcurve IS FALSE 
        if [ ${shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
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
        publishDir "${params.outdir}/bed_shuffle",  mode: 'copy'
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
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks.be*"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks_sc.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks.xls"
        input:
            tuple val(id_ip), path(ip_bed) from SQ30BED_ch
            val(shuffle_percent) from satCurvePCs
        output:
            path("*peaks_sc.bed") into allbed
            path("*peaks.xls") optional true into peaks_xls
            path("*peaks.bed") optional true into peaks_bed
            path("*peaks.bedgraph") optional true into peaks_bg
            stdout into gsize
            val 'ok' into callPeaks_ok
        script:
        """
        ## SELECT N LINES FROM IP BED FILE ACCORDING TO satCurve parameter
        nT=`cat ${ip_bed} |wc -l`
        nPC=`perl -e 'print int('\$nT'*${shuffle_percent})'`
        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## GET GENOME SIZE - BLACKLIST SIZE
        tot_sz=`cut -f3 ${params.fai} |tail -n1`
        bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
        genome_size=`expr \$tot_sz - \$bl_size`
        echo -n \$genome_size
    
        ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter
        for i in {0..${satCurveReps}}; do
            thisName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
            perl ${pickNlines_script} ${ip_bed} \$nPC >\$nPC.tmp
    
            sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
            if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} >> ${params.scratch}/macs2.log
            else
                macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --qvalue ${params.macs_qv} >> ${params.scratch}/macs2.log
            fi

            intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
        
            cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done
 
        ##MERGE BED FILES 
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
  
        ##POSTPROCESS IF params.satcurve IS FALSE
        if [ ${shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
            cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
            ## Calculate strength
            perl ${norm_script} --bed ${id_ip}.peaks.bed \
                 --in ${ip_bed} --out ${id_ip}.peaks.bedgraph --rc --rev_src ${reverse_script}
        fi
        """
    }
}

if ( params.with_control && params.with_idr && params.nb_replicates == "2" ) {
    println ("coucou")
    T1BED_replicate_ch
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1].split('_')[0..-2].join('_'), it[2], it[3] ] }
        .groupTuple(by: [0])
        .groupTuple(by: [1])
        .map { it ->  [ it[0], it[1][0], it[2], it[3] ] }
        .map { it -> it[0,1,2,3].flatten() }
        .set { T1BED_replicate_ch }


    process createPseudoReplicates_ct {
        tag "tmp"
        label 'process_basic'
        publishDir "${params.outdir}/pseudo_replicates",  mode: 'copy'
        input:
            tuple val(id_ip), val(id_ct), file(ip_rep1), file(ip_rep2), file(ct_r1), file(ct_r2) from T1BED_replicate_ch
        output:
            tuple val(id_ip), val(id_ct), file(ip_rep1), file(ip_rep2), file('*ct_pool.bed') into TRUEREP
            tuple val(id_ip), val(id_ct), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed'), file('*ct_pool.bed') into PSREP1
            tuple val(id_ip), val(id_ct), file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed'), file('*ct_pool.bed') into PSREP2
            tuple val(id_ip), val(id_ct), file('*pool_r1.bed'), file('*pool_r2.bed'), file('*ct_pool.bed') into PLREP
        script:
        """
        nlines_r1=\$((`cat ${ip_rep1} | wc -l`/2)) 
        nlines_r2=\$((`cat ${ip_rep2} | wc -l`/2))
        shuf ${ip_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
        shuf ${ip_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
        mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
        mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
        mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
        mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed

        nlines_pool=\$((`cat ${ip_rep1} ${ip_rep2} | wc -l`/2))
        cat ${ip_rep1} ${ip_rep2} > ${id_ip}_pool.bed
        shuf ${id_ip}_pool.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
        mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
        mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed

        cat ${ct_r1} ${ct_r2} > ${id_ct}_ct_pool.bed
        """
    }
}
else if ( !params.with_control && params.with_idr && params.nb_replicates == "2" ) {
    T1BEDrep
        .map { it -> [ it[0].split('_')[0..-4].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it -> it[0,1].flatten() }
        .set { T1BED_replicate_ch }
        //.println()   

     process createPseudoReplicates {
        tag "${id_ip}"
        label 'process_basic'
        publishDir "${params.outdir}/pseudo_replicates",  mode: 'copy'
        input:
            tuple val(id_ip), file(ip_rep1), file(ip_rep2) from T1BED_replicate_ch
        output:
            tuple val(id_ip), file(ip_rep1), file(ip_rep2) into TRUEREP
            tuple val("${id_ip}_psrep1"), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed') into PSREP1
            tuple val("${id_ip}_psrep2"), file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed') into PSREP2
            tuple val("${id_ip}_plrep"), file('*pool_r1.bed'), file('*pool_r2.bed') into PLREP
            tuple val("${id_ip}_pool"), file('*poolT.bed') into POOL 
        script:
        """
        nlines_r1=\$((`cat ${ip_rep1} | wc -l`/2)) 
        nlines_r2=\$((`cat ${ip_rep2} | wc -l`/2))
        shuf ${ip_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
        shuf ${ip_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
        mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
        mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
        mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
        mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed

        nlines_pool=\$((`cat ${ip_rep1} ${ip_rep2} | wc -l`/2))
        cat ${ip_rep1} ${ip_rep2} > ${id_ip}_poolT.bed
        shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
        mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
        mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
        """
    }

    process callPeaksForIDR {
        tag "${id_ip}"
        label 'process_basic'
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy'
        input:
            tuple val(id_ip), file(ip_rep1), file(ip_rep2) from TRUEREP
            tuple val(psrep1), file(r1_pseudorep_r1), file(r1_pseudorep_r2) from PSREP1
            tuple val(psrep2), file(r2_pseudorep_r1), file(r2_pseudorep_r2) from PSREP2
            tuple val(plrep), file(pool_r1), file(pool_r2) from PLREP
            tuple val(pool), file(poolT) from POOL
            val(genome_size) from gsize
        output:
            tuple val(id_ip), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak') into PEAKTRUEREP
            tuple val(psrep1), file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak') into PEAKPSREP1
            tuple val(psrep2), file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak') into PEAKPSREP2
            tuple val(plrep), file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak') into PEAKPLREP
            tuple val(pool), file('*poolT*.regionPeak') into PEAKPOOL
        script:
        """
        echo ${genome_size}
        for file in \$(ls *.bed) ; 
        do
            echo \$file
            name=`basename -- \$file .bed`
            macs2 callpeak -g ${genome_size} -t \$file \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.macs_qv}
            npeaks=`cat \$name_peaks.narrowPeak | wc -l`
	    if [ \$npeaks -gt 100000 ]
	    then
		sort -k 8nr,8nr \$name_peaks.narrowPeak | head -n 100000  > \$name.regionPeak
	    else
		npeaks=\$(( \$npeaks - 1 ))
		sort -k 8nr,8nr \$name_peaks.narrowPeak | head -n \$npeaks  > \$name.regionPeak
	    fi
			
 
        done
        """
    }

    // MAP PEAKTRUEREP CHANNEL
    //PEAKTRUEREP
    //    .map { it -> it[0,1].flatten() }
    //    .set { PEAKTRUEREPSEP }


    process IDRanalysis {
        tag "${id_ip}"
        label 'process_basic'
        conda 'idr=2.0.4.2'
        publishDir "${params.outdir}/idrresults",  mode: 'copy'
        input:
            tuple val(id_ip), file(ip_rep1), file(ip_rep2) from PEAKTRUEREP
            tuple val(psrep1), file(r1_pseudorep_r1), file(r1_pseudorep_r2) from PEAKPSREP1
            tuple val(psrep2), file(r2_pseudorep_r1), file(r2_pseudorep_r2) from PEAKPSREP2
            tuple val(plrep), file(pool_r1), file(pool_r2) from PEAKPLREP
            tuple val(pool), file(poolT) from PEAKPOOL
        output:
            file('*.log') into IDRLOG
            file('*.txt') into IDRTXT
        script:
        """
            python ${encode_idr_script} --peak-type regionPeak \
                --idr-thresh ${params.idr_threshold} \
                --idr-rank ${params.idr_rank} \
                --blacklist ${params.blacklist} \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix r1_pseudorep \
                ${r1_pseudorep_r1} ${r1_pseudorep_r2} ${ip_rep1}
    
            python ${encode_idr_script} --peak-type regionPeak \
                --idr-thresh ${params.idr_threshold} \
                --idr-rank ${params.idr_rank} \
                --blacklist ${params.blacklist} \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix r2_pseudorep \
                ${r2_pseudorep_r1} ${r2_pseudorep_r2} ${ip_rep2}

            python ${encode_idr_script} --peak-type regionPeak \
                --idr-thresh ${params.idr_threshold} \
                --idr-rank ${params.idr_rank} \
                --blacklist ${params.blacklist} \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix truerep \
                ${ip_rep1} ${ip_rep2} ${poolT}

             python ${encode_idr_script} --peak-type regionPeak \
                --idr-thresh ${params.idr_threshold} \
                --idr-rank ${params.idr_rank} \
                --blacklist ${params.blacklist} \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix poolrep \
                ${pool_r1} ${pool_r2} ${poolT}
        """
    }
}

 /*           
     
        .println()
//PSREP1.println()
//PSREP2.println()
//PLREP.println()
}         

/*
process makeSatCurve {
    tag "${outNameStem}"
    label 'process_basic'
    conda "${baseDir}/environment_callpeaks.yml"
    publishDir "${params.outdir}/saturation_curve",  mode: 'copy'
    input:
        file(saturation_curve_data) from allbed.collect()
    output:
        path("*satCurve.tab", emit: table) into satcurve_table
        path("*.png", emit: png) into curve
        val 'ok' into makeSatCurve_ok
    script:
    """
    #Option 1 : original perl script
    perl ${getPeaksBedFiles_script} -tf satCurve.tab -dir ${params.outdir}/saturation_curve/peaks
    #Option 2 : playing with fire & awk
    #echo "reads pc  hs" > satCurve.tab
    #wc -l ${params.outdir}/saturation_curve/peaks/*N*_peaks_sc.bed \
    #    | grep -v "total" | sort -k1n,1n \
    #    | sed 's/ \\+\\([0-9]\\+\\) .\\+\\.N\\([0-9]\\+\\)_\\([.0-9]\\+\\)pc.\\+/\\1 \\2 \\3/' 
    #    | awk '{printf "%d\\t%d\\t%d\\n", \$2, \$3*100, \$1}' >> satCurve.tab
    Rscript ${runSatCurve_script} ${satCurveHS_script} satCurve.tab ${outNameStem}    
    """
}

// GENERAL MULTIQC
process general_multiqc {
    tag "${outNameStem}"
    label 'process_basic'
    conda 'bioconda::multiqc=1.9'
    publishDir "${params.outdir}/multiqc",  mode: 'copy'
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
    multiqc -c ${params.multiqc_configfile} -n ${outNameStem}.multiQC \
        ${params.outdir}/raw_fastqc ${params.outdir}/trim_fastqc \
        ${params.outdir}/samstats ${params.outdir}/bigwig 
    """
}
*/
// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    scrdir.deleteDir()
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

