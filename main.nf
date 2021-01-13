#!/usr/bin/env nextflow
/*
========================================================================================
                        SSDS Pipeline version 2.0
			Adapted from Kevin Brick (version 1.8_NF)
                        Pauline Auffret, 2020
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
PROCESS 1  : check_design (CHECK INPUT DESIGN FILE)
PROCESS 2  : makeScreenConfigFile (MAKE CONFIGURATION FILE FOR FASTQSCREEN)
PROCESS 3  : trimming (USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES)
PROCESS 4  : fastqc (QUALITY CONTROL ON RAW READS USING FASTQC AND FASTQSCREEN)
PROCESS 5  : bwaAlign (USE BWA AND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA)
PROCESS 6  : filterBam (MARK DUPLICATES, REMOVE SUPPLEMENTARY ALIGNMENTS, SORT AND INDEX)
PROCESS 7  : parseITRs (PARSE BAM FILES TO GET THE 5 DIFFERENT SSDS TYPES)
PROCESS 8  : makeDeeptoolsBigWig (GENERATES BIGWIG FILES)
PROCESS 9  : samStats (GENERATES SAMSTATS REPORTS)
PROCESS 10 : toFRBigWig (GENERATES FWD/REV BIGWIG FILES)
PROCESS 11 : makeSSreport (GET INFO FROM QC SSDS REPORT)
PROCESS 12 : ssds_multiqc (MAKE A MULTIQC REPORT FOR SSDS FILES)
PROCESS 13 : shufBEDs (BED SHUFFLING)
PROCESS 14 : callPeaks (PEAK CALLING WITH MACS2)
PROCESS 15 : createPseudoReplicates (CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
PROCESS 16 : callPeaksForIDR (CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
PROCESS 17 : IDRanalysis (PERFORM IDR ANALYSIS ON 4 PAIRS OF REPLICATES OR PSEUDOREPLICATES
PROCESS 18 : makeSatCurve (CREATE SATURATION CURVE)
PROCESS 19 : general_multiqc (GENERATES GENERAL MULTIQC REPORT)
----------------------------------------------------------------------------------------
*/

// Construct help message (option --help)
def helpMessage() { 
    log.info"""
=============================================================================
  SSDS Pipeline version 2.0
=============================================================================
    Usage:

    nextflow run main.nf -c conf/igh.config --inputcsv tests/fastq/input.csv  --name "runtest" --trim_cropR1 36 --trim_cropR2 40 --with_trimgalore -profile conda -resume


=============================================================================

Input data parameters:
    --inputcsv                  FILE    PATH TO INPUT CSV FILE (template and default : ${baseDir}/tests/fastq/input.csv)       
    --genomebase		DIR	PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
    --genome    		STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
    --genomedir			DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
    --genome_name		STRING	REFERENCE GENOME NAME (e.g ".mm10", required if your reference genome is not present in your config file)
    --genome_fasta  		FILE	PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
    --fai			FILE	PATH TO GENOME FAI INDEX FILE (required if your reference genome is not present in your config file)
    --genome2screen 		STRING	GENOMES TO SCREEN FOR FASTQC SCREENING (default : ['mm10','hg19','dm3','dm6','hg38','sacCer2','sacCer3'], comma separated list of genomes to screen reads for contamination, names must correspond to existing genomes in your config file)

Output and temporary directory parameters:                            
    --name      		STRING    ANALYSIS NAME (default : "SSDS_pipeline")      
    --outdir    		DIR       PATH TO OUTPUT DIRECTORY (default : name.outdir)           
    --scratch   		DIR       PATH TO TEMPORARY DIRECTORY (default : scratch)

Pipeline dependencies:
    --src	        	DIR	PATH TO SOURCE DIRECTORY (default : accessoryFiles/SSDS/scripts ; contains perl scripts)
    --custom_bwa        	EXE	PATH TO CUSTOM BWA EXEC (default : accessoryFiles/SSDS/bwa_0.7.12)
    --custom_bwa_ra		EXE	PATH TO CUSTOM BWA_SRA EXEC (default : accessoryFiles/SSDS/bwa_ra_0.7.12)
    --hotspots	        	DIR	PATH TO HOTSPOTS FILES DIRECTORY (default : accessoryFiles/SSDS/hotspots)
    --blacklist                 FILE    PATH TO BLACKLIST BED FILE FOR PEAK CALLING (default : accessoryFiles/SSDS/blacklist/mm10/blackList.bed)
    --multiqc_configfile        FILE    OPTIONAL : PATH TO MULTIQC CUSTOM CONFIG FILE (default : ${baseDir}/multiqc_config.yaml)
QC parameters
    --with_ssds_multiqc		BOOL	RUN SSDS MULTIQC (need a functional conda environment, see multiqc-dev_conda-env parameter ; default : false)
    --multiqc_dev_conda_env     DIR	PATH TO MULTIQC-DEV CONDA ENVIRONMENT (used when --with_ssds-multiqc is true ; default : multiqc_dev)

Trimming parameters:
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
    --trimmomatic_adapters      FILE    PATH TO ADAPTERS FILE FOR TRIMMOMATIC (default ${baseDir}/TruSeq2-PE.fa, special formatting see http://www.usadellab.org/cms/?page=trimmomatic)
    --trimgalore_adapters       FILE    OPTIONAL : PATH TO ADAPTERS FILE FOR TRIMGALORE (default : none)

Mapping parameters:
    --no_multimap               BOOL    If you want to remove multimapping reads from bam (default : false)
    --bamPGline			STRING	bam header (default '@PG\\tID:ssDNAPipeline1.8_nxf_KBRICK')

Bigwig parameter:
    --binsize                   INT     Deeptools binsize parameter (default : 50)

Peak calling parameters:
    --with_control              BOOL    If you are running the analysis with input control files (default : false)
    --satcurve                  BOOL    If you want to plot saturation curve (default : true)
    --reps                      INT     Number of iterations for saturation curve (default : 3)
    --bed_trimqual              INT     Mapping quality threshold for bed filtering (default : 30)
    --macs_bw                   INT     Macs2 callpeak bandwidth parameter (default : 1000)
    --macs_slocal               INT     Macs2 callpeak slocal parameter (default : 5000)
    --macs_extsize              INT     Macs2 callpeak extsize parameter (default : 800)
    --macs_qv                   FLOAT   Macs2 callpeak q-value parameter (default : 0.1)
    --macs_pv                   FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : -1)
    --sctype                    STRING  Saturation curve type (either 'minimal', 'standard' or 'expanded' ; default : 'standard')

Optional IDR analysis parameters:
    --with_idr                  BOOL    If you want to perform IDR analysis, only possible if nb_replicates=2 (default : false)
    --nb_replicates             INT     Number of replicates per sample (default : 2)
    --idr_threshold             FLOAT   idr p-value threshold (default : 0.1)
    --idr_rank                  INT     p.value or q.value (default : p.value)
    --idr_filtering_pattern     STRING  Regex for filtering bed files (default :"chr[1-9X]+")
    --idr_macs_qv               FLOAT   Macs2 callpeak q-value parameter (default : -1)
    --idr_macs_pv               FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : 0.1)
    --idr_chrsz                 FILE    Chromosome sizes file, default : ${baseDir}/accessoryFiles/SSDS/mm10/mm10.chrom.sizes (downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes 2021-01-11)


QC parameters:
    --with_ssds_multiqc		BOOL	RUN SSDS MULTIQC (need a functional conda environment, see multiqc-dev_conda-env parameter ; default : false)
    --multiqc_dev_conda_env     DIR	PATH TO MULTIQC-DEV CONDA ENVIRONMENT (used when --with_ssds-multiqc is true ; default : multiqc_dev)
    --multiqc_configfile        FILE    OPTIONAL : PATH TO MULTIQC CUSTOM CONFIG FILE (default : ${baseDir}/multiqc_config.yaml)

Nextflow Tower parameter:
    -with-tower                 BOOL    Enable job monitoring with Nextflow tower (https://tower.nf/)
    --tower_token               STRING  Nextflow tower key token (see https://tower.nf/ to create your account)
=================================================================================================
""".stripIndent()
}

//***************************************************************************//
//           PRELIMINARY SECTION : GLOBAL VARIABLES ANS SETTINGS             //
//***************************************************************************//
//Show pipeline version
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
// Create scratch directory
scrdir = file("${params.scratch}")
result = scrdir.mkdirs()
println result ? "OK" : "Cannot create directory: $scrdir"

// Custom name variables
def outNameStem = "${params.name}.SSDS.${params.genome}"
def tmpNameStem = "${params.name}.tmpFile.${params.genome}"

// External scripts used in the pipeline
def ITR_id_v2c_NextFlow2_script = "${params.src}/ITR_id_v2c_NextFlow2.pl" //Author Kevin Brick
def ssDNA_to_bigwigs_FASTER_LOMEM_script = "${params.src}/ssDNA_to_bigwigs_FASTER_LOMEM.pl" //Author Kevin Brick
def makeSSMultiQCReport_nextFlow_script = "${params.src}/makeSSMultiQCReport_nextFlow.pl" //Author Kevin Brick
def check_design_script = "${params.src}/check_design.py" // Adapted from nf-core chipseq pipeline version 1.2.1
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


//***************************************************************************//
//                                                                           //
//                          BEGINNING PIPELINE                               //
//                                                                           //
// **************************************************************************//

//***************************************************************************//
//                     SECTION 1 : INPUT SETTINGS                            //
//***************************************************************************//

// PROCESS 1 : check_design (CHECK INPUT DESIGN FILE)
// What it does : checks the conformity and integrity of the input csv file
// Input : 1 csv file with 6 columns [group,replicate,fastq_1,fastq_2,antibody,control]
// Output : 2 csv files : 1 for mapping fastq files to their sample ID ; and 1 for mapping chIP files to corresponding control files if needed.
// External tool : python script ${check_design_script} adapted from nf-core chipseq pipeline.
process check_design {
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
    
//CREATE INPUT CHANNEL WITH SAMPLE ID AND CORRESPONDING FASTQ FILES R1 AND R2
//The resulting channel is composed of 3 elements : [sampleID,file_R1.fq(.gz),file_R2.fq(.gz)]
ch_design_reads_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
    .set { fq_ch }
    //.println()

//CREATE INPUT CHANNEL MAPPING chIP SAMPLE ID AND control SAMPLE ID
//The resulting channel is composed of 5 elements : [sampleID,controlID,antibody,replicate(1/0),multiple(1/0)]
ch_design_controls_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(),row.multipleGroups.toBoolean() ] }
    .set { ch_design_controls_csv }
    //.println()

// PROCESS 2 : makeScreenConfigFile (MAKE CONFIGURATION FILE FOR FASTQSCREEN)
// What it does : generates a configuration file for fastqscreen which will contain the list of genomes to be screened during general QC.
// The list of genomes is defined in the parameter --genome2screen
// Output : text file named conf.fqscreen
process makeScreenConfigFile {
    tag "${outNameStem}" 
    label 'process_basic'
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    output:
        file "checkfile.ok" into fqscreen_conf_ok
    script:
        //Create config file and write header
        def glist=params.genome2screen
        File conf  = new File("${params.outdir}/conf.fqscreen")
        conf.write "This is a config file for FastQ Screen\n\n"
        conf << "THREADS ${task.cpus}\n\n"
        //Output list of genomes to screen in config file
        for (item in glist) {
	    if (params.genomes.containsKey(item)) {
            	fasta=params.genomes[ item ].genome_fasta
            	name=params.genomes[ item ].genome_name
            	conf << "DATABASE  ${name}    ${fasta}\n"
	    }
        }
        //Test if file has been created, if not, exit.
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

//***************************************************************************//
//                           SECTION 2 : TRIMMING                            //
//***************************************************************************//

// PROCESS 3 : trimming (USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES)
// What it does : runs trimmomatic (default) or Trim Galore (if option --with_trimagalore is set) for adapter & quality trimming.
// Several parameters can be set for the trimming, see help section.
// For trimmomatic an adapter file need to be set, and for trim galore the choice is yours.
// Finally hard trimming is done using fastx-trimmer on trimmed fastq files and QC is done with fastqc.
// Input : channel with couple of raw fastq files
// Output : channel with couple of trimmed fastq files and all QC reports from fastqc
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

// MAP FASTQC CHANNEL (Need to flat the raw files R1.fq(.gz) and R2.fq(.gz) for process 4)
// The resulting channel is composed of 3 elements : sampleID, fastq1, fastq2
fqc_ch
    .map { it -> it[0,1].flatten() }
    .set { fqc_tuple }

// PROCESS 4 : fastqc (QUALITY CONTROL ON RAW READS USING FASTQC AND FASTQSCREEN)
// What it does : runs fastqc and fastqscreen on raw reads
// Input : raw reads and config file from fastqscreen created in process 2
// Output : QC reports
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
    // Get fastq files extension (for properly rename files)
    ext=("${read1.getExtension()}")
    """
    # Rename raw files so that they contain the sampleID in the name (will be useful for the QC names)
    if [[ ${ext} == "gz" ]] ; then ext="fastq.gz" ; fi
    mv ${read1} ${sampleId}_raw_R1.${ext}   
    mv ${read2} ${sampleId}_raw_R2.${ext}
    # Run fastqc and fastqscreen
    fastqc -t ${task.cpus} *raw* 
    fastq_screen --threads ${task.cpus} --force --aligner bwa --conf ${params.outdir}/fastqscreen/conf.fqscreen *raw*
    """
}

//***************************************************************************//
//                      SECTION 3 : MAPPING AND PARSING                      //
//***************************************************************************//

// PROCESS 5 : bwaAlign (USE BWA AND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA)
// What it does : aligns trimmed ssds reads to the reference genome
// Input : trimmed reads
// Output : sorted and indexed bam file
// External tool : custom bwa (bwa-ra (bwa rigth align) from original pipeline by K. Brick (2012)
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
    
    ## CHECK IF BAM FILE IS EMPTY AFTER MAPPING (if so, exit)
    samtools flagstat ${sampleId}.sorted.bam > ${sampleId}.sorted.bam.flagstat
    if grep -q "0 + 0 mapped" ${sampleId}.sorted.bam.flagstat
    then
        echo "Bam file ${sampleId}.sorted.bam is empty. Check genome parameter."
        exit 1
    fi

    ## REMOVE MULTIMAPPERS IF OPTION --no_multimap IS TRUE
    ## ! It's important that the final bam files are named *.sorted.bam for the next processes.
    if [[ ${params.no_multimap} == "true" ]]
    then
        samtools view -h ${sampleId}.sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${sampleId}.final.bam
        samtools index ${sampleId}.final.bam
        rm ${sampleId}.sorted.bam ; mv ${sampleId}.final.bam ${sampleId}.sorted.bam
        rm ${sampleId}.sorted.bam.bai ; mv ${sampleId}.final.bam.bai ${sampleId}.sorted.bam.bai
    fi
    """
}

// PROCESS 6 : filterBam (MARK DUPLICATES, REMOVE SUPPLEMENTARY ALIGNMENTS, SORT AND INDEX)
// NOT SURE WHY THIS PROCESS IS USED (from original pipeline) #todo
// What it does : filters bam files (remove unmapped and secondary and mark duplicates)
// Input : sorted bam
// Output : filtered sorted and indexed bam; and unmapped & supplementary bam
process filterBam {
    tag "${sampleId}"
    label 'process_medium'
    publishDir "${params.outdir}/filterbam",  mode: 'copy', pattern: "*.txt"
    input:
        set val(sampleId), file(bam) from sortedbam2filterbam
    output:
        set val(sampleId), file('*.unparsed.bam') into bamAligned
        set val(sampleId), file('*.unparsed.bed') into INPUTBED
        set val(sampleId), file('*.unparsed.bam.bai') into bamIDX 
        set val(sampleId), file('*.unparsed.suppAlignments.bam') into bamAlignedSupp
        set val(sampleId), file('*.unparsed.suppAlignments.bam.bai') into bamIDXSupp
        val 'ok' into filterbam_ok
    script:
    """
    # Remove unmapped and supplementary, then mark duplicates and index
    samtools view -F 2048 -hb ${bam} > ${bam.baseName}.ok.bam
    picard MarkDuplicatesWithMateCigar I=${bam.baseName}.ok.bam O=${bam.baseName}.unparsed.bam \
        PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.MDmetrics.txt MINIMUM_DISTANCE=400 \
        CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${bam.baseName}.unparsed.bam

    # Convert bam to bed
    bedtools bamtobed -i ${bam.baseName}.unparsed.bam > ${bam.baseName}.unparsed.bed

    # Get the unmapped and the supplementary, then mark duplicates and index
    samtools view -f 2048 -hb ${bam} > ${bam.baseName}.supp.bam
    picard MarkDuplicatesWithMateCigar I=${bam.baseName}.supp.bam O=${bam.baseName}.unparsed.suppAlignments.bam \
        PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.suppAlignments.MDmetrics.txt \
        MINIMUM_DISTANCE=400 CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT >& ${params.scratch}/picard.out 2>&1
    samtools index ${bam.baseName}.unparsed.suppAlignments.bam
    """
}

// PROCESS 7 : parseITRs (PARSE BAM FILES TO GET THE DIFFERENT SSDS TYPES)
// THIS PROCESS COULD USE SOME CLEANING AND OPTIMIZATION #todo
// What it does : parse the bam file into 5 types (type1 ssds, type2 ssds, ds, ds_strict, unclassified)
// Input : sorted bam file
// Output : bam and bed files of the 5 types
// External tool : perl scripts from K. Brick (original pipeline, 2012) 
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

//***************************************************************************//
//                           SECTION 4 : BIGWIG                              //
//***************************************************************************//
// PROCESS 8 : makeDeeptoolsBigWig (GENERATES BIGWIG FILES)
// What it does : uses deeptools to generates bigwig files for each of the 5 types of bam
// Input : channel BAMwithIDXdt containing indexed bam files of the 5 types of bam
// Output  bigwig files 
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

// PROCESS 9 : samStats (GENERATES SAMSTATS REPORTS)
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

// PROCESS 10 : toFRBigWig (GENERATES FWD/REV BIGWIG FILES)
// External tool : Perl script from K. Brick (original pipeline, 2012) 
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

//***************************************************************************//
//                           SECTION 5 : SSDS QC                             //
//***************************************************************************//

// PROCESS 11 : makeSSreport (GET INFO FROM QC SSDS REPORT)
// External tool : Perl script from K. Brick (original pipeline, 2012)
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
    // PROCESS 12 : ssds_multiqc (MAKE A MULTIQC REPORT FOR SSDS FILES)
    // External tool : custom multiqc python library and conda environment
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



//***************************************************************************//
//                          SECTION 6 : PEAK CALLING                         //
//***************************************************************************//


// Set saturation curve thresholds callPeaks and makeSatCurve processes
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

// CASE 1 : IF INPUT CONTROL ARE PROVIDED
if (params.with_control) {
    // CREATE CHANNEL LINKING IP TYPE 1 SSDNA BED WITH CONTROL DSDNA BED
    // ( The control bed is the input control filtered bam file )

    T1BED
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { T1BED }

    INPUTBED
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { INPUTBED }

    // The resulting channel is composed of 4 elements : sampleID, controlID, chIP type1 bed file, input filtered bam file
    ch_design_controls_csv
        .combine(T1BED)
        .combine(INPUTBED)
        .filter { it[0] == it[5] && it[1] == it[7] }
        .map { it -> it[0,1,6,8].flatten() }
        .into { T1BED_shuffle_ch ; T1BED_replicate_ch }


    // PROCESS 13 : shufBEDs (BED SHUFFLING)
    // What it does : quality trims and shuffles type1 bed files from ITR parsing 
    // Input : type1 bed files from ITR parsing
    // Ouptut : shuffled and filtered type1 bed files
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
            // Define output files names
            def ip_Q30_bed      = ip_bed.name.replaceFirst(".bed",".IP.q30.bed")
            def ip_Q30_shuf_bed = ip_bed.name.replaceFirst(".bed",".IP.sq30.bed")
            def ct_Q30_bed        = ct_bed.name.replaceFirst(".bed",".CT.q30.bed")
            def ct_Q30_shuf_bed   = ct_bed.name.replaceFirst(".bed",".CT.sq30.bed")
            """
            # Quality trimming
            perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= ${params.bed_trimqual} && \$Q[1] >= ${params.bed_trimqual})' ${ip_bed} >${ip_Q30_bed}
            # perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' ${ct_bed} >${ct_Q30_bed}
            perl -lane '@F = split(/\\t/,\$_); print join("\\t",@F) if (\$F[4] >= ${params.bed_trimqual} )' ${ct_bed} >${ct_Q30_bed}
            # Bed shuffling
            shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
            shuf ${ct_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ct_Q30_shuf_bed}
            """
    }

    // CREATE CHANNEL TO COMBINE SQ30BED AND SHUFFLE PERCENT CHANNEL
    SQ30BED_ch
        .combine(satCurvePCs)
        .set { SQ30BED_satcurve_ch }

    // PROCESS 14 : callPeaks (PEAK CALLING WITH MACS2)
    // What it does : call peaks in type1 bed files using macs2. if parameter satCurveReps is > 1,
    // the process will repeat the peak calling satCurveReps-1 times then the resulting bed files will be merged.
    // if satCurve parameter is true, the process will progressively call peaks in downsampled bed files (according to satCurvePCs parameter)
    // to evaluate the saturation of the samples.
    // Input : Shuffled type1 bed files ; Shuffled controlled bed files  and satCurvePCs percent channel
    // Output : Bed/xls files for peaks called in downsampled bed file; and bed/xls files for peaks called in the whole input bed file.
    // Also outputs bedgraph file for the whole input bed file, and genome size parameters (wil be used in post processes) 
    // External tool : perl script from original pipeline by K. Brick (2012)
    process callPeaks_ct {
        tag "${id_ip}"
        label 'process_basic'
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', pattern: "*peaks*.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks*.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*peaks*.xls"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*.macs2.log"
        publishDir "${params.outdir}/FINALPEAKS",            mode: 'copy', pattern: "*finalPeaks_noIDR.bed"
        publishDir "${params.outdir}/FINALPEAKS",            mode: 'copy', pattern: "*finalPeaks_noIDR.xls"
        input:
            tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed), val(shuffle_percent) from SQ30BED_satcurve_ch
        output:
            path("*peaks*.bed") into allbed
            path("*peaks*.xls") optional true into peaks_xls
            path("*peaks.bedgraph") optional true into peaks_bg
            path("*.macs2.log") into macs2log
            stdout into gsize
            val 'ok' into callPeaks_ok
        script:
        """
        ## SELECT N LINES FROM IP BED FILE ACCORDING TO satCurve parameter 
        # nT is the total number of lines in IP bed file
        nT=`cat ${ip_bed} |wc -l`
        # nPC is the percentage of nT corresponding to the current shuffle_percent parameter
        nPC=`perl -e 'print int('\$nT'*${shuffle_percent})'`
        # Select randomly nPC lines in IP bed file and print them into tmp file then sort tmp file
        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## GET GENOME SIZE - BLACKLIST SIZE
        if [ ${params.blacklist} == "None" ]
        then
            genome_size=`cut -f3 ${params.fai} |tail -n1`
        else
            tot_sz=`cut -f3 ${params.fai} |tail -n1`
            bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
            genome_size=`expr \$tot_sz - \$bl_size`
        fi

        # Export the genome_size value for post processes
        echo -n \$genome_size

        ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter 
        for i in {0..${satCurveReps}}; do
            # Create a unique name for output file
            thisName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
            # Call peaks with macs2
	    if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} > IP-${id_ip}_CT-${ct_bed}.macs2.log 2>&1
            else
                macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize}  --qvalue ${params.macs_qv} > IP-${id_ip}_CT-${ct_bed}.macs2.log 2>&1
            fi
    
            # Filter out peaks from blacklist bed file
            if [ ${params.blacklist} != "None" ]
            then
                intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
                # Filter out mitochondrial peaks
                cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            else
                cut -f1-3 \$thisName'_peaks.narrowPeak' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            fi
 
            # Rename excel file
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done

        ##MERGE BED FILES FROM THE LOOP
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
        
        ##POSTPROCESS IF shuffle_percent is 100% ie if all file has been treated 
        if [ ${shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
            cp ${id_ip}.peaks.bed ${id_ip}.finalPeaks_noIDR.bed
            cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
            cp ${id_ip}.peaks.xls ${id_ip}.finalPeaks_noIDR.xls
            ## Calculate strength
            perl ${norm_script} --bed ${id_ip}.peaks.bed \
                 --in ${ip_bed} --out ${id_ip}.peaks.bedgraph --rc --rev_src ${reverse_script}
        fi
        """
    }


// CASE 2 : IF NO INPUT CONTROL IS PROVIDED
} else {

   T1BED
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { T1BED }


    // PROCESS 13 : shufBEDs (BED SHUFFLING)
    // What it does : quality trims and shuffles type1 bed files from ITR parsing 
    // Input : type1 bed files from ITR parsing
    // Ouptut : shuffled and filtered type1 bed files
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
            // Define output files names
            def ip_Q30_bed      = ip_bed.name.replaceFirst(".bed",".IP.q30.bed")
            def ip_Q30_shuf_bed = ip_bed.name.replaceFirst(".bed",".IP.sq30.bed")
            """
            # Quality trimming
            perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) if (\$Q[0] >= ${params.bed_trimqual} && \$Q[1] >= ${params.bed_trimqual})' ${ip_bed} >${ip_Q30_bed}
            # Bed shuffling
            shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
            """
    }

    // CREATE CHANNEL TO COMBINE SQ30BED AND SHUFFLE PERCENT CHANNEL
    SQ30BED_ch
        .combine(satCurvePCs)
        .set { SQ30BED_satcurve_ch }

    // PROCESS 14 : callPeaks (PEAK CALLING WITH MACS2)
    // What it does : call peaks in type1 bed files using macs2. if parameter satCurveReps is > 1,
    // the process will repeat the peak calling satCurveReps-1 times then the resulting bed files will be merged.
    // if satCurve parameter is true, 
    // the process will progressively call peaks in downsampled bed files (according to satCurvePCs parameter)
    // to evaluate the saturation of the samples.
    // Input : Shuffled type1 bed files and satCurvePCs percent channel
    // Output : Bed/xls files for peaks called in downsampled bed file; and bed/xls files for peaks called in the whole input bed file. 
    // Also outputs bedgraph file for the whole input bed file, and genome size parameters (wil be used in post processes)
    // External tool : perl script from original pipeline by K. Brick (2012)
    process callPeaks {
        tag "${id_ip}"
        label 'process_basic'
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/saturation_curve/peaks", mode: 'copy', pattern: "*peaks*.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*.peaks*.bed"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*.peaks*.xls"
        publishDir "${params.outdir}/peaks",                  mode: 'copy', pattern: "*.macs2.log"
        publishDir "${params.outdir}/FINALPEAKS",             mode: 'copy', pattern: "*finalPeaks_noIDR.bed"
        publishDir "${params.outdir}/FINALPEAKS",             mode: 'copy', pattern: "*finalPeaks_noIDR.xls"
        input:
            //tuple val(id_ip), path(ip_bed) from SQ30BED_ch
            //val(shuffle_percent) from satCurvePCs
            tuple val(id_ip), path(ip_bed), val(shuffle_percent) from SQ30BED_satcurve_ch
        output:
            path("*peaks*.bed") into allbed
            path("*peaks*.xls") optional true into peaks_xls
            path("*peaks*.bedgraph") optional true into peaks_bg
            path("*.macs2.log") into macs2log
            stdout into gsize
            val 'ok' into callPeaks_ok
        script:
        """
        ## SELECT N LINES FROM IP BED FILE ACCORDING TO satCurve parameter
        # nT is the total number of lines in IP bed file
        nT=`cat ${ip_bed} |wc -l`
        # nPC is the percentage of nT corresponding to the current shuffle_percent parameter
        nPC=`perl -e 'print int('\$nT'*${shuffle_percent})'`
        # Select randomly nPC lines in IP bed file and print them into tmp file then sort tmp file
        perl ${pickNlines_script} ${ip_bed} \$nPC > \$nPC.tmp
        sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq > \$nPC.IP.bed
    
        ## CALCULATE GENOME SIZE - BLACKLIST SIZE
        if [ ${params.blacklist} == "None" ]
        then
            genome_size=`cut -f3 ${params.fai} |tail -n1`
        else
            tot_sz=`cut -f3 ${params.fai} |tail -n1`
            bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${params.blacklist} |tail -n1`
            genome_size=`expr \$tot_sz - \$bl_size`
        fi

        # Export the genome_size value for post processes
        echo -n \$genome_size
    
        ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter
        for i in {0..${satCurveReps}}; do
            # Create a unique name for output file
            thisName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
            # Call peaks with macs2
            if [ ${params.macs_pv} != -1 ]; then
	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} > IP-${id_ip}.macs2.log 2>&1
            else
                macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$thisName --nomodel --extsize ${params.macs_extsize} --qvalue ${params.macs_qv} > IP-${id_ip}.macs2.log 2>&1
            fi

            # Filter out peaks from blacklist bed file  
            if [ ${params.blacklist} != "None" ]
            then
                intersectBed -a \$thisName'_peaks.narrowPeak' -b ${params.blacklist} -v >\$thisName'.peaks_sc.noBL'
                # Filter out mitochondrial peaks
                cut -f1-3 \$thisName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            else
                cut -f1-3 \$thisName'_peaks.narrowPeak' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName'_peaks_sc.bed'
            fi

            # Rename excel file
            mv \$thisName'_peaks.xls' \$thisName'_peaks_sc.xls'
        done
 
        ##MERGE BED FILES FROM THE LOOP
        sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
  
        ##POSTPROCESS IF shuffle_percent is 100% ie if all file has been treated
        if [ ${shuffle_percent} == 1.00 ]; then
            mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
            cp ${id_ip}.peaks.bed ${id_ip}.finalPeaks_noIDR.bed
            
            cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
            cp ${id_ip}.peaks.xls ${id_ip}.finalPeaks_noIDR.xls
            ## Calculate strength
            perl ${norm_script} --bed ${id_ip}.peaks.bed \
                 --in ${ip_bed} --out ${id_ip}.peaks.bedgraph --rc --rev_src ${reverse_script}
        fi
        """
    }
}


//***************************************************************************//
//                     SECTION 7 : OPTIONAL IDR ANALYSIS                     //
//***************************************************************************//
// THIS ONLY WORKS WITH 2 REPLICATES IN THIS VERSION. #todo



// CASE 1 : IF INPUT CONTROL ARE PROVIDED
if ( params.with_control && params.with_idr && params.nb_replicates == "2" ) {
 
    // This process aimed to rename control input files with a random id tag,
    // because if the same input is used for 2 replicates and they have the same name
    // it will raise an exception in the following processes   
    process renameInputBed {
        tag "${id_ip}"
        label 'process_basic'
        input:
            tuple val(id_ip), val(id_ct), file(ip_bed), file(ct_bed) from T1BED_replicate_ch
        output:
            tuple val(id_ip), val(id_ct), file(ip_bed), file('*renamed*') into T1BED_replicate_ch_renamed
        script:
        """
        random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
        cat ${ct_bed} > ${ct_bed}_\${random_id}_renamed.bed 
        """
        }

    // CREATE CHANNEL TO MAP SAMPLE ID & CONTROL ID WITH ASSOCIATED REPLICATES (BED FILES
    // The resulting channel is composed of 6 elements : sampleID, controlID, T1BED_R1, T1BED_R2, DSBED_R1, DSBED_R2
    T1BED_replicate_ch_renamed
        .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1].split('_')[0..-2].join('_'), it[2], it[3] ] }
        .groupTuple(by: [0])
        .groupTuple(by: [1])
        .map { it ->  [ it[0], it[1][0], it[2], it[3] ] }
        .map { it -> it[0,1,2,3].flatten() }
        .set { T1BED_replicate_ch_renamed }
        //.println()
    

    // PROCESS 15 : createPseudoReplicates (CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
    // What it does : creates 2 pseudo replicates per true replicates, then pool the true replicates, and creates 2 pseudo replicates from this pool.
    // Input : bed files from type1 aligned SSDNA, chip and control
    // Output : The 2 true replicates ; the pool of the 2 true replicates ; the 4 pseudo replicates from true replicates; 
    // the 2 pseudo replicates from the pool of true replicates, and 1 control file (merge of control files if they are different)
    // So 10 files in total.
    process createPseudoReplicates_ct {
        tag "${id_ip}"
        label 'process_medium'
        publishDir "${params.outdir}/pseudo_replicates",  mode: 'copy', pattern: '*.bed'
        input:
            tuple val(id_ip), val(id_ct), file(ip_rep1), file(ip_rep2), file(ct_r1), file(ct_r2) from T1BED_replicate_ch_renamed
        output:
            //tuple val(id_ip), file(ip_rep1), file(ip_rep2) into TRUEREP
            //tuple val("${id_ip}_psrep1"), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed') into PSREP1
            //tuple val("${id_ip}_psrep2"), file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed') into PSREP2
            //tuple val("${id_ip}_plrep"), file('*pool_r1.bed'), file('*pool_r2.bed') into PLREP
            //tuple val("${id_ip}_pool"), file('*poolT.bed') into POOL
            //tuple val("${id_ct}_ct"), file('*ct_pool.bed') into CTPOOL
            tuple val(id_ip), file(ip_rep1), file(ip_rep2), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed'), \
                file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed'), file('*pool_r1.bed'), file('*pool_r2.bed'), \
                file('*poolT.bed'), file('*ct_pool.bed') into ALLREP
            val 'ok' into createPseudoReplicates_ok
        script:
        """
        # Shuffle then split the 2 original replicates in 2 pseudo replicates
        nlines_r1=\$((`cat ${ip_rep1} | wc -l`/2)) 
        nlines_r2=\$((`cat ${ip_rep2} | wc -l`/2))
        shuf ${ip_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
        shuf ${ip_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
        mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
        mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
        mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
        mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed

        # Pool the 2 original replicates then shuffle then split in 2 pseudo replicates
        nlines_pool=\$((`cat ${ip_rep1} ${ip_rep2} | wc -l`/2))
        cat ${ip_rep1} ${ip_rep2} > ${id_ip}_poolT.bed
        shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
        mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
        mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed

        # Test if input files are the same ; if not, merge them to build a new control file for IDR
        if cmp -s ${ct_r1} ${ct_r2} 
        then 
            cat ${ct_r1} > ${id_ct}_ct_pool.bed
        else
            # not sure if the merging process is good #todo
            cat ${ct_r1} ${ct_r2} | sort -n | unique > ${id_ct}_ct_pool.bed
        fi
        """
    }

    
    // PROCESS 16 : callPeaksForIDR (CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
    // What it does : uses macs2 callPeak to perform peak-calling on all replicates (2), pseudo replicates (4), pool (1), pool pseudo replicates (2)
    // Input : all 10 bed files from true replicates and pseudo replicates creation; and genome size from process 14 
    // Output : 9 regionPeak files from peak calling
    process callPeaksForIDR_ct {
        tag "${id_ip}"
        label 'process_medium'
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*narrowPeak*'
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*regionPeak*'
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*.macs2.log'
        input:
            //tuple val(id_ip), file(ip_rep1), file(ip_rep2) from TRUEREP
            //tuple val(psrep1), file(r1_pseudorep_r1), file(r1_pseudorep_r2) from PSREP1
            //tuple val(psrep2), file(r2_pseudorep_r1), file(r2_pseudorep_r2) from PSREP2
            //tuple val(plrep), file(pool_r1), file(pool_r2) from PLREP
            //tuple val(pool), file(poolT) from POOL
            //tuple val(ct), file(ctpool) from CTPOOL
            tuple val(id_ip), file(ip_rep1), file(ip_rep2), file(r1_pseudorep_r1), file(r1_pseudorep_r2),\
                file(r2_pseudorep_r1), file(r2_pseudorep_r2), file(pool_r1), file(pool_r2), \
                file(poolT), file(ctpool) from ALLREP
            val(genome_size) from gsize
        output:
            file('*narrowPeak*') into ALLPEAKS
            //tuple val(id_ip), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak') into PEAKTRUEREP
            //tuple val(psrep1), file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak') into PEAKPSREP1
            //tuple val(psrep2), file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak') into PEAKPSREP2
            //tuple val(plrep), file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak') into PEAKPLREP
            //tuple val(pool), file('*poolT*.regionPeak') into PEAKPOOL
            tuple val(id_ip), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak'),\
                file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak'),\
                file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak'), \
                file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak'), file('*poolT*.regionPeak') into ALLPEAKSREP
            file('*.macs2.log') into MACS2IDRLOG
            val 'ok' into callPeaksForIDR_ok
        script:
        """
        # Get control file basename
        ctname=`basename -- ${ctpool} .bed`
        # Runs macs2 callpeak on all input bed files
        for file in \$(ls *.bed) ; 
        do
            # Get basename of curent bed file
            name=`basename -- \$file .bed`
            # Check if the current file is not the control input file to not call peaks on control file
            if [ \$ctname != \$name ]
            then
                # Call peaks with macs2, using pvalue or qvalue according to paramters
                echo ${genome_size}
                random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
                mkdir ${params.scratch}/\${random_id}
                if [[ ${params.idr_macs_pv} == -1 && ${params.idr_macs_qv} == -1 ]];
                then
                    macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                    --name \$name --nomodel --extsize ${params.macs_extsize} > \$name.macs2.log 2>&1
                elif [ ${params.idr_macs_pv} != -1 ];
                then
                    macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                    --name \$name --nomodel --extsize ${params.macs_extsize} --pvalue ${params.idr_macs_pv} > \$name.macs2.log 2>&1
                elif [ ${params.idr_macs_qv} != -1 ];
                then
                    macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                    --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.idr_macs_qv} > \$name.macs2.log 2>&1
                fi

                # Cut resulting bed files to 10000 lines max (keep that ? #todo)
                npeaks=`cat \${name}_peaks.narrowPeak | wc -l`
                # Check if the number of peaks is not too low
                if [ \$npeaks -lt 20 ]
                then
                    echo "Not enough peaks in \${name}_peaks.narrowPeak file, to perform IDR analysis, please check your bed files and/or idr parameters."
                    exit 1 
	        elif [ \$npeaks -gt 100000 ]
	        then
		    sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n 100000  > \${name}.regionPeak
	        else
		    npeaks=\$(( \$npeaks - 1 ))
		    sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n \$npeaks  > \${name}.regionPeak
	        fi
            fi
        done
        
        """
    }
}

// CASE 2 : IF NO INPUT CONTROL IS PROVIDED
else if ( !params.with_control && params.with_idr && params.nb_replicates == "2" ) {

    // CREATE CHANNEL TO MAP SAMPLE ID WITH ASSOCIATED REPLICATES (BED FILES)
    T1BEDrep
        .map { it -> [ it[0].split('_')[0..-4].join('_'), it[1] ] }
        .groupTuple(by: [0])
        .map { it -> it[0,1].flatten() }
        .set { T1BED_replicate_ch }
        //.println()   


    // PROCESS 15 : createPseudoReplicates (CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
    // What it does : creates 2 pseudo replicates per true replicates, then pool the true replicates, and creates 2 pseudo replicates from this pool.
    // Input : bed files from type1 aligned SSDNA, chip and control
    // Output : The 2 true replicates ; the pool of the 2 true replicates ; the 4 pseudo replicates from true replicates; 
    // the 2 pseudo replicates from the pool of true replicates.
    // So 9 files in total.
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
            val 'ok' into createPseudoReplicates_ok 
        script:
        """
        # Shuffle then split the 2 original replicates in 2 pseudo replicates
        nlines_r1=\$((`cat ${ip_rep1} | wc -l`/2)) 
        nlines_r2=\$((`cat ${ip_rep2} | wc -l`/2))
        shuf ${ip_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
        shuf ${ip_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
        mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
        mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
        mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
        mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed

        # Pool the 2 original replicates then shuffle then split in 2 pseudo replicates
        nlines_pool=\$((`cat ${ip_rep1} ${ip_rep2} | wc -l`/2))
        cat ${ip_rep1} ${ip_rep2} > ${id_ip}_poolT.bed
        shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
        mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
        mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
        """
    }

    // PROCESS 16 : callPeaksForIDR (CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
    // What it does : uses macs2 callPeak to perform peak-calling on all replicates (2), pseudo replicates (4), pool (1), pool pseudo replicates (2)
    // Input : all 9 bed files from true replicates and pseudo replicates creation; and genome size from process 14 
    // Output : 9 regionPeak files from peak calling
    process callPeaksForIDR {
        tag "${id_ip}"
        label 'process_medium'
        conda "${baseDir}/environment_callpeaks.yml"
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*narrowPeak*'
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*regionPeak*'
        publishDir "${params.outdir}/idrpeaks",  mode: 'copy', pattern: '*.macs2.log'
        input:
            tuple val(id_ip), file(ip_rep1), file(ip_rep2) from TRUEREP
            tuple val(psrep1), file(r1_pseudorep_r1), file(r1_pseudorep_r2) from PSREP1
            tuple val(psrep2), file(r2_pseudorep_r1), file(r2_pseudorep_r2) from PSREP2
            tuple val(plrep), file(pool_r1), file(pool_r2) from PLREP
            tuple val(pool), file(poolT) from POOL
            val(genome_size) from gsize
        output:
            file('*narrowPeak*') into ALLPEAKS
            tuple val(id_ip), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak') into PEAKTRUEREP
            tuple val(psrep1), file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak') into PEAKPSREP1
            tuple val(psrep2), file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak') into PEAKPSREP2
            tuple val(plrep), file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak') into PEAKPLREP
            tuple val(pool), file('*poolT*.regionPeak') into PEAKPOOL
            file('*.macs2.log') into MACS2IDRLOG
            val 'ok' into callPeaksForIDR_ok
        script:
        """
        # Runs macs2 callpeak on all input bed files
        for file in \$(ls *.bed) ; 
        do
            # Get basename to name the outputs files according to input file names
            name=`basename -- \$file .bed`
            random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
            mkdir ${params.scratch}/\${random_id}
            # Call peaks with macs2, using pvalue or qvalue according to paramters
            if [[ ${params.idr_macs_pv} == -1 && ${params.idr_macs_qv} == -1 ]];
            then
                macs2 callpeak -g ${genome_size} -t \$file \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                --name \$name --nomodel --extsize ${params.macs_extsize} \$name.macs2.log 2>&1
            elif [ ${params.idr_macs_pv} != -1 ];
            then
                macs2 callpeak -g ${genome_size} -t \$file \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                --name \$name --nomodel --extsize ${params.macs_extsize} --pvalue ${params.idr_macs_pv} \$name.macs2.log 2>&1
            elif [ ${params.idr_macs_qv} != -1 ];
            then
                macs2 callpeak -g ${genome_size} -t \$file \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir ${params.scratch}/\${random_id} \
                    --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.idr_macs_qv} \$name.macs2.log 2>&1
            fi

            # Cut resulting bed files to 10000 lines max (keep that ? #todo)
            npeaks=`cat \${name}_peaks.narrowPeak | wc -l`
            # Check if the number of peaks is not too low
            if [ \$npeaks -lt 20 ]
            then
                echo "Not enough peaks in \${name}_peaks.narrowPeak file, to perform IDR analysis, please check your bed files and/or idr parameters."
                exit 1 
	    elif [ \$npeaks -gt 100000 ]
	    then
		sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n 100000  > \${name}.regionPeak
	    else
		npeaks=\$(( \$npeaks - 1 ))
		sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n \$npeaks  > \${name}.regionPeak
	    fi
        done
        """
    }
}

// PROCESS 17 : IDRanalysis (PERFORM IDR ANALYSIS ON 4 PAIRS OF REPLICATES OR PSEUDOREPLICATES
// What it does : performs IDR (Irreproducible Discovery Rate) statistical analysis on 4 pairs of replicates :
// IDR1 and IDR2 : pseudoreplicates from true replicates ; IDR3 : true replicates ; IDR4 : pseudoreplicates from pooled true replicates ;
// Input : 9 regionpeak files from idr peak calling  
// Ouptut : log and results files from IDR analysis
// External tool : IDR python script from encode chip-seq pipeline
process IDRanalysis {
    tag "${id_ip}"
    label 'process_medium'
    conda 'idr=2.0.4.2 bioconda::ucsc-bedclip=377 bioconda::bedtools=2.29.2 bioconda::ucsc-bedtobigbed=377' 
    publishDir "${params.outdir}/idrresults",  mode: 'copy'
    input:
        tuple val(id_ip), file(ip_rep1), file(ip_rep2), file(r1_pseudorep_r1), file(r1_pseudorep_r2), \
            file(r2_pseudorep_r1), file(r2_pseudorep_r2), file(pool_r1), file(pool_r2), file(poolT) from ALLPEAKSREP
        //tuple val(id_ip), file(ip_rep1), file(ip_rep2) from PEAKTRUEREP
        //tuple val(psrep1), file(r1_pseudorep_r1), file(r1_pseudorep_r2) from PEAKPSREP1
        //tuple val(psrep2), file(r2_pseudorep_r1), file(r2_pseudorep_r2) from PEAKPSREP2
        //tuple val(plrep), file(pool_r1), file(pool_r2) from PEAKPLREP
        //tuple val(pool), file(poolT) from PEAKPOOL
    output:
        //val(id_ip) into IDIP
        file('*unthresholded-peaks*') into UNTHPEAKS
        file('*.log') into IDRLOG
        file('*.png') into IDRPNG
        file('*.bfilt.gz') into BFILT
        tuple val(id_ip), file('*r1_pseudorep*.*Peak.gz'), file('*r2_pseudorep*.*Peak.gz'), \
            file('*truerep*.*Peak.gz'), file('*poolrep*.*Peak.gz') into IDRPEAKS
        //file('*r1_pseudorep*.*Peak.gz') into R1IDRPEAKS
        //file('*r2_pseudorep*.*Peak.gz') into R2IDRPEAKS
        //file('*truerep*.*Peak.gz') into TRUEREPIDRPEAKS
        //file('*poolrep*.*Peak.gz') into POOLREPIDRPEAKS
    when:
        params.with_idr && params.nb_replicates
    script:
    """
        if [ ${params.blacklist} != "None" ]
        then
            blacklist_params="--blacklist ${params.blacklist}"
        else
            blacklist_params=""
        fi

        #IDR1
        python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
            --idr-thresh ${params.idr_threshold} --idr-rank ${params.idr_rank} \
            --chrsz ${params.idr_chrsz} \$blacklist_params \
            --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
            --prefix r1_pseudorep \
            ${r1_pseudorep_r1} ${r1_pseudorep_r2} ${ip_rep1}

        #IDR2
        python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
            --idr-thresh ${params.idr_threshold} --idr-rank ${params.idr_rank} \
            --chrsz ${params.idr_chrsz} \$blacklist_params \
            --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
            --prefix r2_pseudorep \
            ${r2_pseudorep_r1} ${r2_pseudorep_r2} ${ip_rep2}

        #IDR3
        python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
            --idr-thresh ${params.idr_threshold} --idr-rank ${params.idr_rank} \
            --chrsz ${params.idr_chrsz} \$blacklist_params \
            --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
            --prefix truerep \
            ${ip_rep1} ${ip_rep2} ${poolT}

        #IDR4
        python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
            --idr-thresh ${params.idr_threshold} --idr-rank ${params.idr_rank} \
            --chrsz ${params.idr_chrsz} \$blacklist_params \
            --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
            --prefix poolrep \
            ${pool_r1} ${pool_r2} ${poolT}
        
        #rename blackist filtered file for better handling in the output
        if [ ${params.blacklist} != "None" ]
        then
            mv r1_pseudorep.idr${params.idr_threshold}.bfilt.${params.idr_peaktype}.gz r1_pseudorep.bfilt.gz
            mv r2_pseudorep.idr${params.idr_threshold}.bfilt.${params.idr_peaktype}.gz r2_pseudorep.bfilt.gz
            mv truerep.idr${params.idr_threshold}.bfilt.${params.idr_peaktype}.gz truerep.bfilt.gz
            mv poolrep.idr${params.idr_threshold}.bfilt.${params.idr_peaktype}.gz poolrep.bfilt.gz
        else
            echo "no blacklist provided" > log.bfilt
            gzip log.bfilt
        fi

    """
}

process IDRpostprocess {
    tag "${id_ip}"
    label 'process_basic'
        publishDir "${params.outdir}/idrQC",  mode: 'copy', pattern:'*.log'
        publishDir "${params.outdir}/idrQC",  mode: 'copy', pattern:'*.finalPeaks*'
        publishDir "${params.outdir}/FINALPEAKS",  mode: 'copy', pattern:'*.finalPeaks*'
    input:
        tuple val(id_ip), file(r1_pseudorep_peaks), file(r2_pseudorep_peaks), file(truerep_peaks), file(poolrep_peaks) from IDRPEAKS
        //val(id_ip) from IDIP
        //file(r1_pseudorep_peaks) from R1IDRPEAKS
        //file(r2_pseudorep_peaks) from R2IDRPEAKS
        //file(truerep_peaks) from TRUEREPIDRPEAKS
        //file(poolrep_peaks) from POOLREPIDRPEAKS
    output:
        file('*_IDR_QC.log') into IDRQC 
        file('*.finalPeaks*') into FINALPEAKS
    script:
    """
        n1=`zcat ${r1_pseudorep_peaks} | wc -l`
        n2=`zcat ${r2_pseudorep_peaks} | wc -l`
        np=`zcat ${poolrep_peaks} | wc -l`
        nt=`zcat ${truerep_peaks} | wc -l`

        #Compute Self-consistency Ratio = max(N1,N2) / min(N1,N2)
        if [ \$n1 -gt \$n2 ]
        then
            scr=`python -c "print(\$n1/\$n2)"`
        else
            scr=`python -c "print(\$n2/\$n1)"`
        fi
        echo "Self-consistency Ratio = \$scr" >> ${id_ip}_IDR_QC.log

        #Compute Rescue Ratio = max(Np,Nt) / min(Np,Nt)
        if [ \$np -gt \$nt ]
        then
            rr=`python -c "print(\$np/\$nt)"`
            zcat ${poolrep_peaks} > ${id_ip}_${poolrep_peaks}.finalPeaks_IDR.bed
            #filter for chromY, chromM and blacklist ? #todo
        else
            rr=`python -c "print(\$nt/\$np)"`
            zcat ${truerep_peaks} > ${id_ip}_${truerep_peaks}.finalPeaks_IDR.bed
            #filter for chromY, chromM and blacklist ? #todo
        fi
        echo "Rescue Ratio = \$rr" >> ${id_ip}_IDR_QC.log

        #Evaluate reproducibility flag
        if python -c "exit(0 if \$scr < 2 and \$rr < 2 else 1)"; then
            flag='PASS'
        elif python -c "exit(0 if \$scr < 2 or \$rr < 2 else 1)"; then
            flag='BORDERLINE'
        else
            flag='FAIL'
        fi
        echo "Reproducibility flag : \$flag" >> ${id_ip}_IDR_QC.log
        
        
    """
}
        
        


//***************************************************************************//
//                          SECTION 8 : SAT CURVE                            //
//***************************************************************************//


// PROCESS 18 : makeSatCurve (CREATE SATURATION CURVE)
// What it does : 
// Input :
// Ouptut : 
// External tool :
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
    when:
        params.satcurve
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


//***************************************************************************//
//                          SECTION 9 : GENERAL QC                           //
//***************************************************************************//

/*
if ( !params.with_idr ) {
    process createMissingCheckPointsIDR {
        tag "${outNameStem}"
        label 'process_basic'
        output:
            val 'ok' into createPseudoReplicates_ok
            val 'ok' into callPeaksForIDR_ok
            //val 'ok' into IDRanalysis_ok
        script:
        """
        echo "No IDR analysis is run."
        """
    }
}
if ( !params.satcurve ) {
    process createMissingCheckPointsSatCurve {
        tag "${outNameStem}"
        label 'process_basic'
        output:
            val 'ok' into makeSatCurve_ok
        script:
        """
        echo "No saturation curve is run."
        """
    }
} 

// PROCESS 19 : general_multiqc (GENERATES GENERAL MULTIQC REPORT)
// What it does :
// Input :
// Output : 
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
        //val('shufbed_ok') from shufbed_ok.collect()
        //val('callPeaks_ok') from callPeaks_ok.collect()
        //val('createPseudoReplicates_ok') from createPseudoReplicates_ok.collect()
        //val('callPeaksForIDR_ok') from callPeaksForIDR_ok.collect()
       //val('IDRanalysis_ok') from IDRanalysis_ok.collect()
        //val('makeSatCurve_ok') from makeSatCurve_ok.collect()
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


//***************************************************************************//
//                                                                           //
//                          END OF PIPELINE !!                               //
//                                                                           //
// **************************************************************************//

// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    scrdir.deleteDir()
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

