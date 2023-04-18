#!/usr/bin/env nextflow
/*
========================================================================================
                        hotSSDS pipeline version 2.0
			Adapted from Kevin Brick (version 1.8_NF)
                        Pauline Auffret, 2020-2021
========================================================================================
 hotSSDS pipeline
 #### Homepage / Documentation
 Adapted from version 1.8_NF (Kevin Brick)
 Original pipelines :
 https://github.com/kevbrick/SSDSnextflowPipeline
 https://github.com/kevbrick/callSSDSpeaks
 Version 2.0 (Pauline Auffret) :
 https://github.com/jajclement/ssdsnextflowpipeline
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Single-Stranded-DNA-Sequencing (SSDS) Pipeline : Align, Parse and call Peaks in ssDNA
Pipeline overview:

*******************************************
*****   SECTION 1 : INPUT SETTINGS      ***
*******************************************
PROCESS 1  : check_design (CHECK INPUT DESIGN FILE)
PROCESS 2  : makeScreenConfigFile (MAKE CONFIGURATION FILE FOR FASTQSCREEN)
*******************************************
***   SECTION 2 : TRIMMING              ***
*******************************************
PROCESS 3  : crop (HARD TRIMMING AND QUALITY CONTROL ON RAW READS USING FASTQC AND FASTQSCREEN)
PROCESS 4  : trimming (USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM AND REMOVE ADAPTERS FROM RAW SEQUENCES)
*******************************************
***   SECTION 3 : MAPPING AND PARSING   ***
*******************************************
PROCESS 5  : bwaAlign (USE BWA AND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA)
PROCESS 6  : filterBam (MARK DUPLICATES, REMOVE SUPPLEMENTARY ALIGNMENTS, SORT AND INDEX)
PROCESS 7  : parseITRs (PARSE BAM FILES TO GET THE DIFFERENT SSDS TYPES)
PROCESS 8  : makeBigwig (GENERATE NORMALIZED BIGWIG FILES FOR T1 AND T1+T2 BED FILES)
*******************************************
***   SECTION 4 : BIGWIG                ***
*******************************************
PROCESS 9  : makeBigwigReplicates (optional, GENERATE NORMALIZED BIGWIG FILES FOR MERGED REPLICATES T1 AND T1+T2 BED FILES)
PROCESS 10 : makeDeeptoolsBigWig (optional, GENERATES BIGWIG FILES, COVERAGE AND CUMULATIVE COVERAGE PLOTS)
PROCESS 11 : toFRBigWig (optional, GENERATES FWD/REV BIGWIG FILES)
*******************************************
***   SECTION 5 : PEAK CALLING          ***
*******************************************
PROCESS 12 : shufBEDs (BED SHUFFLING)
PROCESS 13 : callPeaks (PEAK CALLING WITH MACS2)
*******************************************
***   SECTION 6 : SSDS QC               ***
*******************************************
PROCESS 14 : samStats (GENERATES SAMSTATS REPORTS)
PROCESS 15 : makeSSreport (COMPUTE STATS FROM SSDS PARSING)
PROCESS 16 : makeFingerPrint (MAKE DEEPTOOLS FINGERPRINT PLOTS)
PROCESS 17 : computeFRIP (Optional, COMPUTE FRIP SCORE FOR PARSED BAM FILE FOR NEW GENOMES)
PROCESS 18 : ssds_multiqc (MAKE MULTIQC REPORT FOR SSDS FILES)
*******************************************
***   SECTION 7 : OPTIONAL IDR ANALYSIS ***
*******************************************
PROCESS 19 : createPseudoReplicates (optional, CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
PROCESS 20 : callPeaksForIDR (optional, CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
PROCESS 21 : IDRanalysis (optional, PERFORM IDR ANALYSIS ON 4 PAIRS OF REPLICATES OR PSEUDOREPLICATES)
PROCESS 22 : IDRpostprocess (optional, IDR PEAKS POST PROCESSING)
*******************************************
***   SECTION 8 : PEAK POST PROCESSING  ***
*******************************************
PROCESS 23 : normalizePeaks (CENTER AND NORMALIZE PEAKS)
PROCESS 24 : mergeFinalPeaks (MERGE PEAKS FROM REPLICATES)
PROCESS 25 : makeSatCurve (optional, CREATE SATURATION CURVE)
*******************************************
***   SECTION 9 : GENERAL QC            ***
*******************************************
PROCESS 26 : general_multiqc (GENERATES GENERAL MULTIQC REPORT)

----------------------------------------------------------------------------------------
*/

// Construct help message (option --help)
def helpMessage() { 
    log.info"""
=============================================================================
  hotSSDS pipeline version 2.0 : Align, parse and call hotspots from SSDNA
=============================================================================
    Usage:

    nextflow run main.nf -c conf/cluster.config --params_file conf/mm10.json --inputcsv tests/fastq/input.csv  --name "runtest" --trim_cropR1 36 --trim_cropR2 40 --with_trimgalore -profile singularity -resume

    Runs with Nextflow DSL1 v20.10.0
=============================================================================
Input data parameters:
    --inputcsv                  FILE    PATH TO INPUT CSV FILE (template and default : ${baseDir}/tests/fastq/input.csv)
    -params_file		FILE	PATH TO PARAMETERS JSON FILE (template and default : ${baseDir}/conf/mm10.json)
    --genomebase                DIR     PATH TO REFERENCE GENOMES 
    --genome                    STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
    --genomedir                 DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
    --genome_name               STRING  REFERENCE GENOME NAME (e.g ".mm10", required if your reference genome is not present in your config file)
    --genome_fasta              FILE    PATH TO GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
    --fai                       FILE    PATH TO GENOME FAI INDEX FILE (required if your reference genome is not present in your config file)
    --genome2screen             STRING  GENOMES TO SCREEN FOR FASTQC SCREENING (default : ['mm10','hg19','dm3','dm6','hg38','sacCer2','sacCer3'], comma separated list of genomes to screen reads for contamination, names must correspond to existing genomes in your config file)
    --chrsize                   FILE    Chromosome sizes file, default : ${baseDir}/data/mm10/mm10.chrom.sizes (downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes 2021-01-11)
    --hotspots                  DIR     PATH TO HOTSPOTS FILES DIRECTORY (set to "None" if none provided ; default :  ${baseDir}/data/hotspots/mm10/hotspots)
    --blacklist                 FILE    PATH TO BLACKLIST BED FILE FOR PEAK CALLING AND IDR (set to "None" if none provided ; default : ${baseDir}/data/blacklist/mm10/blackList.bed)

Output and temporary directory parameters:
    --name                      STRING  ANALYSIS NAME (default : "hotSSDS_pipeline")
    --outdir                    DIR     PATH TO OUTPUT DIRECTORY (default : ${baseDir}/{params.name}.outdir/02_results")
    --publishdir_mode           STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)

Pipeline dependencies:
    --src                       DIR     PATH TO SOURCE DIRECTORY (default : ${baseDir}/bin ; contains perl scripts)
    --custom_bwa                EXE     PATH TO CUSTOM BWA EXEC (default : ${baseDir}/bin/bwa_0.7.12)
    --custom_bwa_ra             EXE     PATH TO CUSTOM BWA_SRA EXEC (default : ${baseDir}/bin/bwa_ra_0.7.12)

Trimming parameters:
    --with_trimgalore           BOOL    Use trim-galore instead of Trimmomatic for quality trimming process (default : false)
    --trimgalore_adapters       FILE    trim-galore : PATH TO ADAPTERS FILE (default : none)
    --trimg_quality             INT     trim-galore : minimum quality (default 10)
    --trimg_stringency          INT     trim-galore : trimming stringency (default 6)
    --trim_minlen               INT     trimmomatic : minimum length of reads after trimming (default 25)
    --trim_cropR1               INT     fastx : Cut the R1 read to that specified length (default 50)
    --trim_cropR2               INT     fastx : Cut the R2 read to that specified length (default 50)
    --trim_slidingwin           STRING  trimmomatic : perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default "4:15")
    --trim_illumina_clip        STRING  trimmomatic : Cut adapter and other illumina-specific sequences from the read (default "2:20:10")
    --trimmomatic_adapters      FILE    PATH TO ADAPTERS FILE FOR TRIMMOMATIC (default ${baseDir}/data/TruSeq2-PE.fa, special formatting see http://www.usadellab.org/cms/?page=trimmomatic)

Mapping parameters:
    --with_multimap             BOOL    Keep multimapping reads from bam (default : false)
    --bamPGline                 STRING  bam header (default '@PG\tID:ssDNAPipeline2.0_PAUFFRET')
    --filtering_flag            INT     SAM flag for filtering bam files (default : 2052 ; see https://broadinstitute.github.io/picard/explain-flags.html)
    --picard_min_distance       INT     Picard parameter for marking duplicates (--MINIMUM_DISTANCE) :  width of the window to search for duplicates of a given alignment, default : -1 (twice the first read's read length)
    --picard_optdup_distance    INT     Picard parameter for marking duplicates (--OPTICAL_DUPLICATE_PIXEL_DISTANCE) : The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform (HiSeq2500). For the patterned flowcell models (Novaseq 6000), 2500 is more appropriate, default : 100
    --get_supp                  BOOL    Publish bam files for supplementary aligments, default : false

Bigwig parameter:
    --bigwig_profile            STRING  Bigwig profile using  bedtools (normalization by total library size) : "T1" will produce bigwig for T1 bed files only, one per replicates ; "T12" will also produce bigwig for merged T1+T2, one per replicates ; "T1rep" will also produce T1 bigwig for merged replicates ; "T12rep" will also produce T1+T2 bigwig for merged replicates (default : "T1")
    --kbrick_bigwig             BOOL    Compute bigwig files ; FR bigwig files and coverage plots as in original pipeline by Kevin Brick using deeptools with FPKM normalization (default : false)
    --binsize                   INT     Deeptools binsize parameter (used only if kbrick_bigwig is TRUE ; default : 50)

Peak calling parameters:
    --with_control              BOOL    Use input control files for peak calling analysis (default : false)
    --satcurve                  BOOL    Plot saturation curve (default : false)
    --sctype                    STRING  Saturation curve type (either 'minimal', 'standard' or 'expanded' ; default : 'standard')
    --reps                      INT     Number of iterations for saturation curve (default : 3)
    --bed_trimqual              INT     Mapping quality threshold for bed filtering (default : 30)
    --macs_bw                   INT     Macs2 callpeak bandwidth parameter (default : 1000)
    --macs_slocal               INT     Macs2 callpeak slocal parameter (default : 5000)
    --macs_extsize              INT     Macs2 callpeak extsize parameter (default : 800)
    --macs_qv                   FLOAT   Macs2 callpeak q-value parameter (default : 0.1)
    --macs_pv                   FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : -1)
    --no_chrY                   BOOL    Filter out chromosomeY peaks from final peak bed files (default : true)

Optional IDR analysis parameters (ENCODE procedure, see https://github.com/ENCODE-DCC/chip-seq-pipeline2) :
    --with_idr                  BOOL    Perform IDR analysis, only possible if nb_replicates=2 (default : false)
    --nb_replicates             INT     Number of replicates per sample (default : 2)
    --idr_peaktype              STRING  The peak file format for IDR (narrowPeak, regionPeak or broadPeak, default : "regionPeak")
    --idr_setup			STRING	Threshold profile for idr. This will define the thresholds for true replicates, pool replicates, self replicates r1 and self replicates r2. Profile "auto" is based on ENCODE guidelines and profile "custom" allows to set custom thresholds (see parameters --idr_threshold_r1 --idr_threshold_r2 --idr_threshold_truerep and --idr_threshold_poolrep ; default : auto)
    --idr_threshold_r1		FLOAT	idr threshold for self replicates r1 (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_r2		FLOAT   idr threshold for self replicates r2 (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_truerep	FLOAT   idr threshold for true replicates (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_poolrep	FLOAT   idr threshold for pooled replicates (used if --idr_setup is "custom" only ; default : 0.01)
    --idr_rank                  INT     p.value or q.value (default : p.value)
    --idr_filtering_pattern     STRING  Regex for filtering bed files (default :"chr[1-9X]+" for mouse ; set ".*" to keep everything)
    --idr_macs_qv               FLOAT   Macs2 callpeak q-value parameter (default : -1)
    --idr_macs_pv               FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : 0.1)

QC parameters:
    --with_ssds_multiqc         BOOL    RUN SSDS MULTIQC (default : true)
    --multiqc_configfile        FILE    OPTIONAL : PATH TO MULTIQC CUSTOM CONFIG FILE (default : ${baseDir}/conf/multiqc_config.yaml)

Nextflow Tower parameter:
    -with-tower                 BOOL    Enable job monitoring with Nextflow tower (https://tower.nf/)

Singularity images parameters:
    --get_sif			BOOL	[REQUIRE INTERNET ACCESS] Check and download singularity images if necessary (if true, pipeline will stops after download. Once downloading has been done, relaunch pipeline with false ; default: false)
    --url_sif                   URL     URL TO PUBLIC SINGULARITY IMAGES REPOSITORY (default : https://zenodo.org/record/7783473/files)

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

// Custom name variables
def global_name = "${params.name}.$workflow.runName"
def tmp_name = "${params.name}.tmpFile"
def macs2_params = "pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}"
def idr_params = "${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}"
def control_status = params.with_control ? "with-input" : "without-input"

// Create empty files for optional results
File emptyfile1  = new File("${params.outdir}/02_results/qc/empty1.txt")
File emptyfile2  = new File("${params.outdir}/02_results/qc/empty2.txt")
File emptyfile3  = new File("${params.outdir}/02_results/qc/empty3.txt")
File emptyfile4  = new File("${params.outdir}/02_results/qc/empty4.txt")
File emptyfile5  = new File("${params.outdir}/02_results/qc/empty5.txt")
File emptyfile6  = new File("${params.outdir}/02_results/qc/empty6.txt")

// Creta logo channel for analysis multiqc report
logo_ch = file(params.logo, checkIfExists: true)

// External scripts used in the pipeline
ITR_id_v2c_NextFlow2_ch = Channel.fromPath("${params.src}/ITR_id_v2c_NextFlow2.pl", checkIfExists: true) //Author Kevin Brick
ssDNA_to_bigwigs_FASTER_LOMEM_ch = Channel.fromPath("${params.src}/ssDNA_to_bigwigs_FASTER_LOMEM.pl", checkIfExists: true) //Author Kevin Brick
makeSSMultiQCReport_nextFlow_ch = Channel.fromPath("${params.src}/makeSSMultiQCReport_nextFlow.pl", checkIfExists: true) //Author Kevin Brick
check_design_ch = Channel.fromPath("${params.src}/check_design.py", checkIfExists: true) // Adapted from nf-core chipseq pipeline version 1.2.1
pickNlines_ch = Channel.fromPath("${params.src}/pickNlines.pl", checkIfExists: true) //Author Kevin Brick
satCurveHS_ch = Channel.fromPath("${params.src}/satCurveHS.R", checkIfExists: true) //Author Kevin Brick
norm_ch = Channel.fromPath("${params.src}/normalizeStrengthByAdjacentRegions.pl", checkIfExists: true) //Author Kevin Brick
reverse_ch = Channel.fromPath("${params.src}/reverseStrandsForOriCalling.pl", checkIfExists: true) //Author Kevin Brick
getPeaksBedFiles_ch = Channel.fromPath("${params.src}/getPeaksBedFiles.pl", checkIfExists: true) //Author Kevin Brick, script adapted by Pauline Auffret
runSatCurve_ch = Channel.fromPath("${params.src}/runSatCurve.R", checkIfExists: true) //Author Pauline Auffret
encode_idr_ch= Channel.fromPath("${params.src}/encode-dcc_chip-seq-pipeline2_src/encode_task_idr.py", checkIfExists: true) //Author Jin Lee from https://github.com/ENCODE-DCC/chip-seq-pipeline2
get_frip_ch = Channel.fromPath("${params.src}/compute_frip.py", checkIfExists: true) //Author Pauline Auffret
plot_ssds_stat_ch = Channel.fromPath("${params.src}/plotSSDSqc.R", checkIfExists: true) //Author Pauline Auffret
runPlotSSDSqc_ch = Channel.fromPath("${params.src}/runPlotSSDSqc.R", checkIfExists: true) //Author Pauline Auffret

// If execution profile is set to test, create inputcsv for test
if (workflow.profile ==~ /.*test.*/) {
    println(workflow.profile)
    File input = new File("${baseDir}/tests/fastq/input.csv")
    input.write "group,replicate,fastq_1,fastq_2,antibody,control\n"
    input << "TEST_IP,1,${baseDir}/tests/fastq/ssdsLong.100k.R1.fastq,${baseDir}/tests/fastq/ssdsLong.100k.R2.fastq,antiDMC1,"
}

// Check if input csv file exists
if (params.inputcsv) { println("Checking input sample file...") ; input_ch = file(params.inputcsv, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

// Check if other input files/directories exist
if (params.chrsize) { 
    println("Checking chrsize input file...") ; 
    Channel
        .fromPath(params.chrsize, checkIfExists: true) 
        .into { chrsize_ch ; chrsize_ch2 }
} else { 
    chrsize_ch = Channel.fromPath("${workDir}") 
    chrsize_ch2 = Channel.fromPath("${workDir}")
} ; println("Ok")

if (params.hotspots && params.hotspots != "None") { 
    println("Checking hotspots directory...") ; 
    Channel
        .fromPath(params.hotspots, checkIfExists: true)
        .into { hotspots_ch ; hotspots_ch2 } 
} else { 
    hotspots_ch = Channel.fromPath("${workDir}") 
    hotspots_ch2 = Channel.fromPath("${workDir}")
} ; println("Ok")

if (params.blacklist && params.blacklist != "None") { 
    println("Checking blacklist file...") ; 
    Channel
        .fromPath(params.blacklist, checkIfExists: true) 
        .into { blacklist_ch ; blacklist_ch2 } ; 
} else { 
    blacklist_ch = Channel.fromPath("${workDir}") 
    blacklist_ch2 = Channel.fromPath("${workDir}")
} ; println("Ok") 

if (params.multiqc_configfile) { 
    println("Checking multiqc_configfile...") ; 
    multiqc_configfile_ch = Channel.fromPath(params.multiqc_configfile, checkIfExists: true) 
} ; else { multiqc_configfile_ch = Channel.fromPath("${workDir}") } ; println("Ok")

if (params.trimgalore_adapters) { 
    println("Checking trimgalore_adapters file...") ; 
    trimgalore_adapters_ch = Channel.fromPath("${params.trimgalore_adapters}", checkIfExists: true)
} ; else { trimgalore_adapters_ch = Channel.fromPath("${workDir}") } ; println("Ok")

if (params.trimmomatic_adapters) { 
    println("Checking trimmomatic_adapters file...") ;
    trimmomatic_adapters_ch = Channel.fromPath("${params.trimmomatic_adapters}", checkIfExists: true)
} ; else { trimmomatic_adapters_ch = Channel.fromPath("${workDir}") } ; println("Ok")

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genome.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Create channel containing Singularity images names
sif_ch = Channel.from(["bam-box-1.0","bam-box_1.0.sif"],["bigwig-box-1.0","bigwig-box-1.0.sif"],["frip-box-1.0","frip-box_1.0.sif"],["idr-box-1.0","idr-box_1.0.sif"],["multiqc-box-1.0","multiqc-box_1.0.sif"],["peak-calling-box-1.0","peak-calling-box_1.0.sif"],["plot-box-1.0","plot-box_1.0.sif"],["python-3.8","python-3.8.sif"],["ssds-qc-box-1.0","ssds-qc-box_1.0.sif"],["trimming-box-1.0","trimming-box_1.0.sif"])

// Define genome variables
params.genome_fasta = params.genome ? params.genomes[ params.genome ].genome_fasta ?: false : false
params.genomedir = params.genome ? params.genomes[ params.genome ].genomedir ?: false : false
params.genome_name = params.genome ? params.genomes[ params.genome ].genome_name ?: false : false
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false

// Check if the bwa index is present in the genome folder
println("Checking bwa index files...") ; fastaindex=file("${params.genome_fasta}" + ".amb", checkIfExists: true) ; println("Bwa index .amb file OK")
println("Checking bwa index files...") ; fastaindex=file("${params.genome_fasta}" + ".ann", checkIfExists: true) ; println("Bwa index .ann file OK")
println("Checking bwa index files...") ; fastaindex=file("${params.genome_fasta}" + ".bwt", checkIfExists: true) ; println("Bwa index .bwt file OK")
println("Checking bwa index files...") ; fastaindex=file("${params.genome_fasta}" + ".pac", checkIfExists: true) ; println("Bwa index .pac file OK")
println("Checking bwa index files...") ; fastaindex=file("${params.genome_fasta}" + ".sa", checkIfExists: true) ; println("Bwa index .sa file OK >>>> All BWAindex OK!!!")

// Define channels containing genome and indexes
bwa_ch = Channel.fromPath(params.custom_bwa, checkIfExists: true)
bwa_ra_ch = Channel.fromPath(params.custom_bwa_ra, checkIfExists: true)

Channel
    .fromPath(params.genome_fasta, checkIfExists: true)
    .into { fasta_ch ; fasta_ch2 ; fasta_ch3 }
Channel
    .fromPath("${params.genomedir}/*.amb", checkIfExists: true)
    .set {index_amb_ch }
Channel
    .fromPath("${params.genomedir}/*.ann", checkIfExists: true)
    .set { index_ann_ch }
Channel
    .fromPath("${params.genomedir}/*.bwt", checkIfExists: true)
    .set { index_bwt_ch }
Channel
    .fromPath("${params.genomedir}/*.pac", checkIfExists: true)
    .set { index_pac_ch }
Channel
    .fromPath("${params.genomedir}/*.sa", checkIfExists: true)
    .set { index_sa_ch }
Channel
    .fromPath(params.fai, checkIfExists: true)
    .into { index_fai_ch ; index_fai_ch2 ; index_fai_ch3 ; index_fai_ch4 }

// Check input parameters conformity
if(params.bigwig_profile != "T1" && params.bigwig_profile != "T12" && params.bigwig_profile != "T1rep" && params.bigwig_profile != "T12rep") {
    println("Error : --bigwig_profile parameter must be either T1 ; T12 ; T1rep or T12rep.")
    exit 0
}
if(params.idr_setup != "auto" && params.idr_setup != "custom") {
    println("Error : --idr_setup parameter must be either auto or custom.")
    exit 0
}
if(params.satcurve) {
    if(params.sctype != "minimal" && params.sctype != "standard" && params.sctype != "expanded") {
        println("Error : --sctype parameter must be either minimal ; standard or expanded.")
        exit 0
    }
}
if(params.publishdir_mode!="copy" && params.publishdir_mode!="symlink" && params.publishdir_mode!="rellink" && params.publishdir_mode!="link" && params.publishdir_mode!="copyNoFollow" && params.publishdir_mode!="move") {
    println("Error : --publishdir_mode must be symlink, rellink, link, copy, copyNoFollow,move, see https://www.nextflow.io/docs/latest/process.html")
    exit 0
}
if(params.satcurve) {
    if(params.reps < 3) {
        println("Error : when --satcurve is true, --reps must be greater than or equal to 3")
        exit 0
    }
}
if(params.with_idr) {
    if(params.nb_replicates != "2") {
        println("Error : when --with_idr is true, --nb_replicates must be equal to 2")
        exit 0
    }
}
 
// Pipeline parameters information
def paramsSection() {
log.info """
==========================================================================
   hotSSDS pipeline version 2.0 : Align, parse and call hotspots from SSDNA  
==========================================================================
** Main parameters ** 
Run name                       : ${params.name}
Input file                     : ${params.inputcsv}
Reference genome               : ${params.genome}
Blacklist file                 : ${params.blacklist}
Hotspots file                  : ${params.hotspots}
Output directory               : ${params.outdir}
Use input control files        : ${params.with_control}
Keep multimappers              : ${params.with_multimap}
Number of replicates           : ${params.nb_replicates}
Bigwig profile                 : ${params.bigwig_profile}
Generate deeptools bigwig      : ${params.kbrick_bigwig}
Exclude chromosome Y peaks     : ${params.no_chrY}
Run IDR analysis               : ${params.with_idr}
Plot saturation curve          : ${params.satcurve}
Use trimgalore                 : ${params.with_trimgalore}
R1 hard trimming               : ${params.trim_cropR1}
R2 hard trimming               : ${params.trim_cropR2}
Use SSDS custom multiQC        : ${params.with_ssds_multiqc}

** Trimming parameters **
Trimmomatic adapter file       : ${params.trimmomatic_adapters}
Trimmomatic quality threshold  : ${params.trimg_quality}
Trimmomatic stringency         : ${params.trimg_stringency}
Trimmomatic sliding window     : ${params.trim_slidingwin}
Trimmomatic illumina clip      : ${params.trim_illuminaclip}
Trimming min length            : ${params.trim_minlen}
Fastq-screen genomes to screen : ${params.genome2screen}
MultiQC configuration file     : ${params.multiqc_configfile}
Trimgalore adapter file        : ${params.trimgalore_adapters}

** Mapping and filtering parameters **
Samtools filtering flag        : ${params.filtering_flag}
Quality threshold              : ${params.bed_trimqual}
Picard min dist for markdup    : ${params.picard_min_distance}
Picard optical dup pixel dist  : ${params.picard_optdup_distance}
Publish supplementary algnmts  : ${params.get_supp}

** Bigwig parameter **
Deeptools bigwig bin size      : ${params.binsize}

** Peak calling parameters **
Macs2 bandwidth parameter      : ${params.macs_bw}
Macs2 slocal parameter         : ${params.macs_slocal}
Macs2 extsize parameter        : ${params.macs_extsize}
Macs2 q-value                  : ${params.macs_qv}
Macs2 p-value                  : ${params.macs_pv}

** Saturation curve parameters **
Saturation curve profile       : ${params.sctype}
Number of iterations           : ${params.reps}

** IDR (Irreproducible Discovery Rate) analysis parameters **
Peak type                      : ${params.idr_peaktype}
Pseudo-rep 1 IDR threshold     : ${params.idr_threshold_r1}
Pseudo-rep 2 IDR threshold     : ${params.idr_threshold_r2}
True rep IDR threshold         : ${params.idr_threshold_truerep}
Pool rep IDR threshold         : ${params.idr_threshold_poolrep}
IDR rank value                 : ${params.idr_rank}
IDR macs2 q-value              : ${params.idr_macs_qv}
IDR macs2 p-value              : ${params.idr_macs_pv}
IDR filtering pattern          : ${params.idr_filtering_pattern}
 
** Pipeline dependencies **
Source directory               : ${params.src}
Chromosome sizes file          : ${params.chrsize}
BWA bin                        : ${params.custom_bwa}
BWA-ra bin                     : ${params.custom_bwa_ra}
BAM custom header              : ${params.bamPGline}
Publish directory mode         : ${params.publishdir_mode}
Download Singularity images ?  : ${params.get_sif}
Singularity images URL	       : ${params.sif_url}         
    
""".stripIndent()
}
// Print parameters in log report
paramsSection()

//***************************************************************************//
//                                                                           //
//                          BEGINNING PIPELINE                               //
//                                                                           //
// **************************************************************************//

// Before running pipeline, check and download singularity images when necessary if option get_sif is true
if (params.get_sif) {
    process get_sif {
        tag "${name}"
        label 'internet_access'
        publishDir "${params.outdir}/01_logs/get_sif", mode: params.publishdir_mode, pattern: '*.log'
        input:
            tuple val(dir), val(name) from sif_ch
        output:
            path '*.log'
            val 'ok' into sif_ok
        script:
        """
        echo "Checking ${baseDir}/containers/${dir}..." > get_sif_${name}.log
        mkdir -p ${baseDir}/containers/${dir}
        if [ -f ${baseDir}/containers/${dir}/${name} ];
        then
            echo "${name} already in ${baseDir}/containers/${dir}." >> get_sif_${name}.log
        else
            echo "Downloading Singularity ${name} images..." >> get_sif_${name}.log
            wget '${params.sif_url}/${name}' -O ${baseDir}/containers/${dir}/${name} -o ./wget_${name}.log || rm -f ${baseDir}/containers/${dir}/${name}
        fi
        if grep -q "ERROR" wget_bam-box_1.0.sif.log || grep -q "failed" wget_bam-box_1.0.sif.log ;
        then
            echo "${name} file download failed. Check sif_url (${params.sif_url}) parameter."
            exit 1
        fi
        """
    }
}

else {

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
        label 'process_basic'
        label "python"
        publishDir "${params.outdir}/02_results/qc/design/pipeline_info",	mode: params.publishdir_mode, pattern: '*.csv'
        publishDir "${params.outdir}/01_logs/check_design",			        mode: params.publishdir_mode, pattern: '*.log'
        input:
            path(design) from input_ch
            path(check_design_script) from check_design_ch 
        output:
            path 'design_reads.csv' into ch_design_reads_csv
            path 'design_controls.csv' into ch_raw_design_controls_csv
            path '*.log'
        script:
        """
        # Run python script to check the design file
        python ${check_design_script} ${design} design_reads.csv design_controls.csv >& ${design}_check_design.log
        # Check if --with_control parameter value is consistent with csv input file
        # Indeed, if no control files provided, nb_line_ctrl_file will only contain one line (the header)
        nb_line_ctrl_file=`cat design_controls.csv | wc -l`
        if [[ \$nb_line_ctrl_file == "1" && "${params.with_control}" == "true" ]]
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
    
    //CREATE INPUT CHANNEL MAPPING ChIP SAMPLE ID AND control SAMPLE ID
    //The resulting channel is composed of 5 elements : [sampleID,controlID,antibody,replicate(1/0),multiple(1/0)]
    ch_raw_design_controls_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(),row.multipleGroups.toBoolean() ] }
        .into { ch_design_controls_csv ; ch_design_controls_csv2 }
        //.println()
    
    // PROCESS 2 : makeScreenConfigFile (MAKE CONFIGURATION FILE FOR FASTQSCREEN)
    // What it does : generates a configuration file for fastqscreen which will contain the list of genomes to be screened during general QC.
    // The list of genomes is defined in the parameter --genome2screen
    // Output : text file named conf.fqscreen
    process makeScreenConfigFile {
        tag "${global_name}" 
        label 'process_basic'
        output:
            path "conf.fqscreen" into conf_fqscreen_ch
            path "checkfile.ok" into fqscreen_conf_ok
        script:
            //Create config file and write header
            def glist=params.genome2screen
            File conf  = new File("${params.outdir}/conf.fqscreen")
            conf.write "This is a config file for FastQ Screen\n\n"
            conf << "THREADS ${task.cpus}\n\n"
            def present = 'false'
            //Output list of genomes to screen in config file
            for (item in glist) {
    	        if (params.genomes.containsKey(item)) {
                	fasta=params.genomes[ item ].genome_fasta
                	name=params.genomes[ item ].genome_name
                	conf << "DATABASE  ${name}    ${fasta}\n"
                }
                if (item == params.genome_name) {
                    present = 'true' 
    	        }
            }
            //Add params.genome to genome to screen if not present in default list
            if (present == 'false') {
                conf << "DATABASE  ${params.genome_name}    ${params.genome_fasta}\n"
            }
            //Test if file has been correctly created, if not, exit pipeline.
    	    """
    	    if [ -f "${params.outdir}/conf.fqscreen" ]; then
        	echo "${params.outdir}/conf.fqscreen exists." > checkfile.ok
                mkdir -p "${params.outdir}/02_results/qc/fastqscreen"
                cp ${params.outdir}/conf.fqscreen .
                mv "${params.outdir}/conf.fqscreen" "${params.outdir}/02_results/qc/fastqscreen/"
    	    else
                echo "The configuration file for fastqscreen could not be generated. Please check genome2screen parameter."
    	    exit 1
    	    fi
    	    """
    }
    
    //***************************************************************************//
    //                           SECTION 2 : TRIMMING                            //
    //***************************************************************************//
    // MAP FASTQ CHANNEL (Need to flat the raw files R1.fq(.gz) and R2.fq(.gz) for process 4)
    // The resulting channel is composed of 3 elements : sampleID, fastq1, fastq2
    fq_ch
        .map { it -> it[0,1].flatten() }
        .set { fq_tuple }
    
    // PROCESS 3 : crop (HARD TRIMMING AND QUALITY CONTROL ON RAW READS USING FASTQC AND FASTQSCREEN)
    // What it does : uses fastx_trimmer to hard trim raw reads then runs fastqc and fastqscreen on cropped reads
    // Input : raw reads and config file from fastqscreen created in process 2
    // Output : QC reports
    process crop {
        tag "${sampleId}"
        label 'process_low'
        label 'trimming'
        publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.html"
        publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.zip"
        publishDir "${params.outdir}/02_results/qc/fastqscreen", mode: params.publishdir_mode, pattern: "*.png"
        publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.txt"
        publishDir "${params.outdir}/01_logs/trimming",	         mode: params.publishdir_mode, pattern: "*.log"
        input:
            tuple val(sampleId), path(read1), path(read2) from fq_tuple
            path(conf_fqscreen) from conf_fqscreen_ch 
            path(ok) from fqscreen_conf_ok
        output:
            tuple val("${sampleId}"), path("*_crop_R1.fastq.gz"), path("*_crop_R2.fastq.gz") into fqcrop_tuple, broute
    	    path('*') into raw_fastqc_ch
            val 'ok' into fastqc_ok
            path('*.log')
        script:
        // Get fastq files extension (for properly rename files)
        ext=("${read1.getExtension()}")
        """
        # Rename raw files so that they contain the sampleID in the name (will be useful for the QC names)
        if [[ ${ext} == "gz" ]]
        then 
            ext="fastq.gz"
            mv ${read1} ${sampleId}_raw_R1.${ext}
            mv ${read2} ${sampleId}_raw_R2.${ext}
        else
            ext="fastq.gz"
            gzip -c ${read1} > ${sampleId}_raw_R1.${ext}
            gzip -c ${read2} > ${sampleId}_raw_R2.${ext}
        fi
        
        # Crop raw reads
        zcat ${sampleId}_raw_R1.${ext} | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${sampleId}_crop_R1.fastq.gz
        zcat ${sampleId}_raw_R2.${ext} | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${sampleId}_crop_R2.fastq.gz
    
        # Run fastqc and fastqscreen
        fastqc -t ${task.cpus} ${sampleId}_crop_R1.fastq.gz ${sampleId}_crop_R2.fastq.gz >& ${sampleId}_raw-fastqc.log 2>&1
     
        fastq_screen --threads ${task.cpus} --force --aligner bwa --conf ${conf_fqscreen} ${sampleId}_crop_R1.fastq.gz ${sampleId}_crop_R2.fastq.gz >& ${sampleId}_fastqscreen.log 2>&1
        """
    }
    
    // CREATE CHANNEL FOR TRIMMING PROCESS
    // The resulting channel is composed of 5 elements : sample_id, raw_fastq_R1, raw_fastq_R2, trim-galore-adapters-file, trimmomatic-adapters-file
    // Note : if no adapter file provided, the elements trim-galore-adapters-file and trimmomatic-adapters-file contains path to current directory
    // (This is to avoid cardinality problem with input channels, and avoid empty channels) 
    trimgalore_adapters_ch
        .combine(trimmomatic_adapters_ch)
        .set { adapters_ch }
    
    fqcrop_tuple
        .combine(adapters_ch)
        .set { fqcrop_tuple_adapters }
    
    // PROCESS 4 : trimming (USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM AND REMOVE ADAPTERS FROM RAW SEQUENCES)
    // What it does : runs trimmomatic (default) or Trim Galore (if option --with_trimagalore is set) for adapter & quality trimming.
    // Several parameters can be set for the trimming, see help section.
    // For trimmomatic an adapter file need to be set, and for trim galore the choice is yours.
    // Finally QC is done on trimmed reads with fastqc.
    // Input : channel fqcrop_tuple with couple of cropped fastq files
    // Output : channel trim_ch with couple of trimmed fastq files and all QC reports from fastqc
    process trimming {
        tag "${sampleId}" 
        label 'process_low'
        label 'trimming'
        publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*_report.txt"
        publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*trim*.html"
        publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*trim*.zip"
        publishDir "${params.outdir}/02_results/trimming/trim_fastq",  mode: params.publishdir_mode, pattern: "*crop_trim_R1.fastq.gz"
        publishDir "${params.outdir}/02_results/trimming/trim_fastq",  mode: params.publishdir_mode, pattern: "*crop_trim_R2.fastq.gz"
        publishDir "${params.outdir}/01_logs/trimming",		           mode: params.publishdir_mode, pattern: "*.log"
        input:
            tuple val(sampleId), path(cropread1), path(cropread2), path(trimgalore_adapters), path(trimmomatic_adapters) from fqcrop_tuple_adapters
        output:
            tuple val("${sampleId}"), path('*crop_trim_R1.fastq.gz'), path('*crop_trim_R2.fastq.gz') into trim_ch, toot
            path('*') into trim_fastqc_ch 
            val 'ok' into trimming_ok
            path('*.log')
        script:
        if (params.with_trimgalore && params.trimgalore_adapters)
    	"""
    	trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
                    --cores ${task.cpus} --adapter "file:${trimgalore_adapters}" --gzip --paired --basename ${sampleId} ${cropread1} ${cropread2} >& ${sampleId}_trimgalore.log 2>&1
        mv ${sampleId}_val_1.fq.gz ${sampleId}_crop_trim_R1.fastq.gz
        mv ${sampleId}_val_2.fq.gz ${sampleId}_crop_trim_R2.fastq.gz
            
        fastqc -t ${task.cpus} ${sampleId}_crop_trim_R1.fastq.gz ${sampleId}_crop_trim_R2.fastq.gz >& ${sampleId}_trim-fastqc.log 2>&1
        """
        else if (params.with_trimgalore && !params.trimgalore_adapters)
        """
        trim_galore --quality ${params.trimg_quality} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
                --cores ${task.cpus} --gzip --paired --basename ${sampleId} ${cropread1} ${cropread2} >& ${sampleId}_trimgalore.log 2>&1
        mv ${sampleId}_val_1.fq.gz ${sampleId}_crop_trim_R1.fastq.gz
        mv ${sampleId}_val_2.fq.gz ${sampleId}_crop_trim_R2.fastq.gz
    
        fastqc -t ${task.cpus} ${sampleId}_crop_trim_R1.fastq.gz ${sampleId}_crop_trim_R2.fastq.gz >& ${sampleId}_trim-fastqc.log 2>&1
    	"""
    	else
    	"""		
        trimmomatic PE -threads ${task.cpus} ${cropread1} ${cropread2} \
                ${sampleId}_crop_trim_R1.fastq.gz R1_unpaired.fastq.gz \
                ${sampleId}_crop_trim_R2.fastq.gz R2_unpaired.fastq.gz \
                ILLUMINACLIP:${trimmomatic_adapters}:${params.trim_illuminaclip} SLIDINGWINDOW:${params.trim_slidingwin} \
                MINLEN:${params.trim_minlen} >& ${sampleId}_crop_trim_${global_name}_trimmomatic_report.txt 2>&1
            
        fastqc -t ${task.cpus} ${sampleId}_crop_trim_R1.fastq.gz ${sampleId}_crop_trim_R2.fastq.gz >& ${sampleId}_trim-fastqc.log 2>&1
        """
    }
    
    
    //***************************************************************************//
    //                      SECTION 3 : MAPPING AND PARSING                      //
    //***************************************************************************//
    // CREATE CHANNEL FOR MAPPING PROCESS
    // The resulting channel contains 13 elements :
    // sample_id, trimmed_faq_R1, trimmed_fq_R2, bwa, bwa_ra, genome_fasta, 7 bwa index files
    
    bwa_ch
      .combine(bwa_ra_ch)
      .set { bwa_all_ch }
    
    fasta_ch
      .combine(index_amb_ch)
      .combine(index_ann_ch)
      .combine(index_bwt_ch)
      .combine(index_pac_ch)
      .combine(index_sa_ch)
      .combine(index_fai_ch)
      .set { all_bwa_index_ch }
    
    trim_ch
      .combine(bwa_all_ch)
      .combine(all_bwa_index_ch)
      .set { bwa_trim_ch }
    
    // PROCESS 5 : bwaAlign (USE BWA AND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA)
    // What it does : aligns trimmed ssds reads to the reference genome
    // Input : channel trim_ch with trimmed reads
    // Output : sorted and indexed bam file
    // External tool : custom bwa (bwa-ra (bwa rigth align) from original pipeline by K. Brick (2012))
    process bwaAlign {
        tag "${sampleId}"
        label 'bam'
        publishDir "${params.outdir}/02_results/bwa/bam",      mode: params.publishdir_mode, pattern: "*.sorted.bam*"
        publishDir "${params.outdir}/02_results/qc/flagstat",  mode: params.publishdir_mode, pattern: "*.flagstat"
        publishDir "${params.outdir}/01_logs/bwa",  	       mode: params.publishdir_mode, pattern: "*.log"
        publishDir "${params.outdir}/01_logs/bwa",  	       mode: params.publishdir_mode, pattern: "*.out"
        input:
            tuple val(sampleId), path(fqR1), path(fqR2), path(bwa), path(bwa_ra), path(fasta), path(amb), path(ann), path(bwt), path(pac), path(sa), path(fai) from bwa_trim_ch 
        output:
            tuple val(sampleId), path('*.sorted.bam') into SORTEDBAM
            path('*.flagstat') into bwa_flagstat_ch
            path('*.log')
    	    path('*.out')
            val 'ok' into bwa_ok
        script:
        """
        # Align R1 reads (fully complementary to the genome)  with bwa aln
        ./${bwa} aln -t ${task.cpus} ${fasta} ${fqR1} \
                > ${tmp_name}.R1.sai 2>${tmp_name}.R1.sai.log
        # Align R2 reads (potentially contain fill-in ITR part at the end of the 5')
        # with bwa-ra aln (custom version of bwa that search for the longest mappable suffix in the query)
        ./${bwa_ra} aln -t ${task.cpus} ${fasta} ${fqR2} \
                > ${tmp_name}.R2.sai 2>${tmp_name}.R2.sai.log
        # Generate alignments in the SAM format given R1 and R2 alignments
        ./${bwa} sampe ${fasta} ${tmp_name}.R1.sai ${tmp_name}.R2.sai ${fqR1} \
                ${fqR2} >${tmp_name}.unsorted.sam 2> ${tmp_name}.unsorted.sam.log
        # Convert SAM to BAM file
        picard SamFormatConverter I=${tmp_name}.unsorted.sam O=${tmp_name}.unsorted.tmpbam \
                VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${tmp_name}.unsorted.picardSFC.out 2>&1
        # Sort and index sam file
        picard SortSam I=${tmp_name}.unsorted.tmpbam O=${sampleId}.sorted.bam SO=coordinate \
                VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${sampleId}.sorted.picardSS.out 2>&1
        samtools index ${sampleId}.sorted.bam 
        
        ## CHECK IF BAM FILE IS EMPTY AFTER MAPPING (if so, exit pipeline)
        samtools flagstat ${sampleId}.sorted.bam > ${sampleId}.sorted.bam.flagstat
        if grep -q "^0 + 0 mapped" ${sampleId}.sorted.bam.flagstat
        then
            echo "Bam file ${sampleId}.sorted.bam is empty. Check genome parameter."
            exit 1
        fi
        
        ## MARK DUPLICATES TO HAVE STATISTICS ABOUT DUPLICATION
        picard -Xms8g -Xmx8g MarkDuplicatesWithMateCigar I=${sampleId}.sorted.bam O=${sampleId}.sorted.md.bam \
            PG=Picard2.9.2_MarkDuplicates M=${sampleId}.sorted.md.txt MINIMUM_DISTANCE=${params.picard_min_distance} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.picard_optdup_distance} REMOVE_DUPLICATES=false \
            CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${sampleId}.sorted.md.out 2>&1
        samtools flagstat ${sampleId}.sorted.md.bam > ${sampleId}.sorted.md.bam.flagstat
     
        ## REMOVE MULTIMAPPERS IF OPTION --with_multimap IS FALSE
        ## ! It's important that the final bam files are named *.sorted.bam for the next processes.
        if [[ ${params.with_multimap} == "false" ]]
        then
            # from https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
            # samtools view -h ${sampleId}.sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${sampleId}.final.bam
            samtools view -h ${sampleId}.sorted.bam | grep -v "XT:A:R" | samtools view -b | samtools sort -o ${sampleId}.final.bam
            samtools index ${sampleId}.final.bam
            rm ${sampleId}.sorted.bam ; mv ${sampleId}.final.bam ${sampleId}.sorted.bam
            rm ${sampleId}.sorted.bam.bai ; mv ${sampleId}.final.bam.bai ${sampleId}.sorted.bam.bai
        fi
        """
    }
    
    // PROCESS 6 : filterBam (MARK DUPLICATES, REMOVE SUPPLEMENTARY ALIGNMENTS, SORT AND INDEX)
    // What it does : filters bam files (remove unmapped and secondary and mark duplicates)
    // Input : sorted bam
    // Output : filtered sorted and indexed bam; and unmapped & supplementary bam
    process filterBam {
        tag "${sampleId}"
        label 'process_medium'
        label 'bam'
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/bed",  mode: params.publishdir_mode, pattern: "*.bed*"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/bam",  mode: params.publishdir_mode, pattern: "*.bam*"
        publishDir "${params.outdir}/01_logs/filterbam/flag_${params.filtering_flag}/log",         mode: params.publishdir_mode, pattern: "*.out"
        input:
            tuple val(sampleId), file(bam) from SORTEDBAM
        output:
            tuple val(sampleId), file('*.unparsed.bam') into FILTEREDBAM
            path('*.unparsed.bam') into UNPARSEDBAM
            path('*.unparsed.bam.bai') into UNPARSEDBAI
            tuple val(sampleId), file('*.unparsed.bed') into INPUTBED
            set val(sampleId), file('*.unparsed.bam'), file('*.unparsed.bam.bai') 
            set val(sampleId), file('*.unparsed.suppAlignments.bam') into bamAlignedSupp optional true
            set val(sampleId), file('*.unparsed.suppAlignments.bam.bai') into bamIDXSupp optional true
            val 'ok' into filterbam_ok
            val 'ok' into filterbam_tofingerprint
    	    file('*.txt')
    	    file('*.out')
        script:
        """
        # Remove unmapped and supplementary, then mark duplicates and index
        samtools view -F ${params.filtering_flag} -f 2 -q ${params.bed_trimqual} -hb ${bam} > ${bam.baseName}.ok.bam
        picard -Xms8g -Xmx8g MarkDuplicatesWithMateCigar I=${bam.baseName}.ok.bam O=${bam.baseName}.unparsed.bam \
            PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.MDmetrics.txt MINIMUM_DISTANCE=${params.picard_min_distance} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.picard_optdup_distance} REMOVE_DUPLICATES=true \
            CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam.baseName}.unparsed.picardMD.out 2>&1
        samtools index ${bam.baseName}.unparsed.bam
    
        # Convert bam to bed (usefull for the control input files which will not be parsed into 5 bed files through parseITRs process)
        bedtools bamtobed -i ${bam.baseName}.unparsed.bam > ${bam.baseName}.unparsed.bed
    
        if ${params.get_supp} :
        then
            # Get the supplementary, then mark duplicates and index
            samtools view -f 2048 -hb ${bam} > ${bam.baseName}.supp.bam
            picard -Xms8g -Xmx8g MarkDuplicatesWithMateCigar I=${bam.baseName}.supp.bam O=${bam.baseName}.unparsed.suppAlignments.bam \
                PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.suppAlignments.MDmetrics.txt \
                MINIMUM_DISTANCE=${params.picard_min_distance} CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.picard_optdup_distance} REMOVE_DUPLICATES=true \
                VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam.baseName}.supp.picardMD.out 2>&1
            samtools index ${bam.baseName}.unparsed.suppAlignments.bam
        fi
        """
    }
    
    // IF RUNNING WITH INPUT CONTROL FILES, CREATE CHANNEL WITH ONLY BAM FILES FOR IP SAMPLES AND NOT CONTROL SAMPLES
    // SO THAT CONTROL SAMPLES ARE NOT PROCESSED THROUGH PARSEITRS PROCESS
    if (params.with_control) {
        ch_design_controls_csv2
            .combine(FILTEREDBAM)
            .map { it -> [ it[0], it[1], it[5], it[5].split('_')[0..-2].join('_'), it[6] ] } 
            .filter { it[0] == it[3] }
            .map { it -> it[2,4] }
            .set { FILTEREDBAM  }
    }
    
    // CREATE CHANNEL FOR parseITRs process
    fasta_ch2
      .combine(ITR_id_v2c_NextFlow2_ch)
      .set { fasta_and_ITR_script }
    
    FILTEREDBAM
      .combine(fasta_and_ITR_script)
      .set { FILTEREDBAM_fasta_script }
    
    
    
    // PROCESS 7 : parseITRs (PARSE BAM FILES TO GET THE DIFFERENT SSDS TYPES)
    // THIS PROCESS COULD USE SOME CLEANING AND OPTIMIZATION #todo
    // What it does : parse the bam file into 5 types (type1 ssds, type2 ssds, ds, ds_strict, unclassified)
    // Input : filtered bam file
    // Output : bam and bed files of the 5 types
    // External tool : perl scripts from K. Brick (original pipeline, 2012) 
    process parseITRs {
        tag "${sampleId}"
        label 'bam'
        label 'process_high'
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/type1/bed",       mode: params.publishdir_mode, pattern: "*.ssDNA_type1.bed"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/type1/bam",       mode: params.publishdir_mode, pattern: "*f.ssDNA_type1.bam*"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/type2/bed",       mode: params.publishdir_mode, pattern: "*.ssDNA_type2.bed"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/type2/bam",       mode: params.publishdir_mode, pattern: "*f..ssDNA_type1.bam*"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/dsDNA/bed",       mode: params.publishdir_mode, pattern: "*.dsDNA*.bed"    
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/dsDNA/bed",       mode: params.publishdir_mode, pattern: "*f.dsDNA*.bam*"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/unclassified/bed",mode: params.publishdir_mode, pattern: "*.unclassified.bed"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/unclassified/bam",mode: params.publishdir_mode, pattern: "*f.unclassified.bam*"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/flagstat",        mode: params.publishdir_mode, pattern: "*.flagstat"
        publishDir "${params.outdir}/02_results/qc/flagstat",							 mode: params.publishdir_mode, pattern: "*.flagstat"
        publishDir "${params.outdir}/01_logs/parseITRs/flag_${params.filtering_flag}",         mode: params.publishdir_mode, pattern: "*.out"
        publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/norm_factors",    mode: params.publishdir_mode, pattern: "*norm_factors.txt"
        input:
            tuple val(sampleId), path(bam), path(fasta), path(ITR_id_v2c_NextFlow2_script) from FILTEREDBAM_fasta_script 
        output:
            tuple val(sampleId), file("${bam}"), file('*.ssDNA_type1.bed'), file('*.ssDNA_type2.bed'), \
                file('*.dsDNA.bed'), file('*.dsDNA_strict.bed'), file('*.unclassified.bed')  into ITRBED
            tuple val(sampleId), file("${bam}"), file("${bam}.bai"), file('*.f.ssDNA_type1.bam'), file('*.f.ssDNA_type1.bam.bai'), file('*.f.ssDNA_type2.bam'), \
                file('*.f.ssDNA_type2.bam.bai'), file('*.f.dsDNA.bam'), file('*.f.dsDNA.bam.bai'), file('*.f.dsDNA_strict.bam'), \
                file('*.f.dsDNA_strict.bam.bai'), file('*.f.unclassified.bam'), file('*.f.unclassified.bam.bai')  into ITRBAM
            tuple val(sampleId), file('*ssDNA_type1.bed') into T1BED, T1BEDrep
            tuple val(sampleId), file('*.f.*.bam'),file('*.f.*.bam.bai') into BAMwithIDXfr, BAMwithIDXss, BAMwithIDXdt mode flatten
            tuple val(sampleId), file('*.f.ssDNA_type1.bam'), file('*.f.ssDNA_type1.bam.bai') into T1BAMwithIDXfr, T1BAMwithIDXdt mode flatten
            tuple val(sampleId), file('*norm_factors.txt'), file('*.ssDNA_type1.bed'), \
                    file('*.ssDNA_type12.bed') into BEDtoBW, BEDtoBWrep
            file('*.flagstat') into itr_flastat_ch
    	file('*.out')
            val 'ok' into parseitr_ok
        script:
        """
        # Parse filtered bam file into 5 types
        perl ${ITR_id_v2c_NextFlow2_script} ${bam} ${fasta} >& ${sampleId}_parseITR.log 2>&1
        # Sort resulting bed files (the bed files will be used for peak calling)
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
    
        # Get library size
        libsizeT1=`wc -l  *.ssDNA_type1.bed | grep -v "_random"| awk '{print \$1}'`
        libsizeT2=`wc -l  *.ssDNA_type2.bed | grep -v "_random"| awk '{print \$1}'`
        libsizeDS=`wc -l  *.dsDNA.bed | grep -v "_random"| awk '{print \$1}'`
        libsizeDSstrict=`wc -l  *.dsDNA_strict.bed | grep -v "_random"| awk '{print \$1}'`
        libsizeUN=`wc -l  *.unclassified.bed | grep -v "_random"| awk '{print \$1}'`
        # Keep les unclassified ? #todo
        libsizeTotal=`expr \$libsizeT1 + \$libsizeT2 + \$libsizeDS + \$libsizeUN`
        echo \$libsizeTotal > ${sampleId}_norm_factors.txt
    
        # Merge T1 and T2 bed files
        cat ${bam}.ssDNA_type1.bed ${bam}.ssDNA_type2.bed | sort -k1,1 -k2,2n > ${bam}.ssDNA_type12.bed
     
        # Convert sam to bam files (the bam files with be used for bigwig/bedgrah generation if params.kbrick_bigwig is true)
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
        
        # Sort bam files
        picard SortSam I=${bam}.ssDNA_type1.US.bam  O=${bam}.f.ssDNA_type1.bam \
            SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam}.ssDNA_type1.picardSS.out 2>&1
        picard SortSam I=${bam}.ssDNA_type2.US.bam  O=${bam}.f.ssDNA_type2.bam \
            SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam}.ssDNA_type2.picardSS.out 2>&1
        picard SortSam I=${bam}.dsDNA.US.bam O=${bam}.f.dsDNA.bam \
            SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam}.dsDNA.picardSS.out 2>&1
        picard SortSam I=${bam}.dsDNA_strict.US.bam O=${bam}.f.dsDNA_strict.bam \
            SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam}.dsDNA_strict.picardSS.out 2>&1
        picard SortSam I=${bam}.unclassified.US.bam O=${bam}.f.unclassified.bam \
            SO=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam}.unclassified.picardSS.out 2>&1
        
        # Index bam files
        samtools index ${bam}
        samtools index ${bam}.f.ssDNA_type1.bam
        samtools index ${bam}.f.ssDNA_type2.bam
        samtools index ${bam}.f.dsDNA.bam
        samtools index ${bam}.f.dsDNA_strict.bam
        samtools index ${bam}.f.unclassified.bam
        
        # Remove temporary bam and sam files
        rm -rf *.US.bam
        rm -rf *.sam
    
        ## CHECK IF TYPE 1 BAM FILE IS EMPTY AFTER PARSING
        samtools flagstat ${bam}.f.ssDNA_type1.bam > ${bam}.f.ssDNA_type1.bam.flagstat
        if grep -q "^0 + 0 mapped" ${bam}.f.ssDNA_type1.bam.flagstat
        then
            echo "Type 1 bam file ${bam}.f.ssDNA_type1.bam is empty. Check genome parameter."
            exit 1
        fi
        """
    }
    
    //***************************************************************************//
    //                           SECTION 4 : BIGWIG                              //
    //***************************************************************************//
    // CREATE CHANNEL FOR makeBigwig process
    BEDtoBW
      .combine(chrsize_ch)
      .set { BEDtoBW_chrsize_ch }
    
    // PROCESS 8 : makeBigwig (GENERATE NORMALIZED BIGWIG FILES FOR T1 AND T1+T2 BED FILES)
    // What it does : First compute the normalization factors using total library sizes,
    // then use bedtools to generates bedgraph then bigwig files for T1 bed file, and
    // for merged T1+T2 bed file if parameter bigwig_profile is set to T12 or T12rep
    // Input : channel BEDtoBW containing total libray size, T1 and T1+T2 bed files 
    // from parseITRs process
    // Output : Bigwig files for T1 (and T1+T2 optional)
    process makeBigwig {
        tag "${sampleId}"
        label 'process_basic'
        label 'bigwig'
        publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",      	mode: params.publishdir_mode, pattern: "*.bw"
        publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",      	mode: params.publishdir_mode, pattern: "*.bedgraph"
        publishDir "${params.outdir}/01_logs/bigwig/${params.bigwig_profile}",		mode: params.publishdir_mode, pattern: "*.log"
        input:
            tuple val(sampleId), path(libsizeTotal), path(ssDNA_type1_bed), path(ssDNA_type12_bed), path(chrsize) from BEDtoBW_chrsize_ch
        output:
            file('*.bw')
            file('*.bedgraph')
            file('*.log')
            val('ok') into makeBigwig_ok
        script:
        """
        ## Compute normalization factor (i.e, for T1 : normfactor=librarySizeT1/librarySizeTotal*1000000)
        # Get total library sizes
        libsizeTot=`cat ${libsizeTotal}`
    	
        # Get T1 library sizes
        libsizeT1=`wc -l ${ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
    
        # Get normalization factors
        normfactT1=`python -c "print(round((1000000/(\$libsizeT1)),5))"`
        normfactTot=`python -c "print(round((1000000/(\$libsizeTot)),5))"`
    
        # Compute bedgraph and bigwig with normfactT1
        genomeCoverageBed -i ${ssDNA_type1_bed} -g ${chrsize} -scale \$normfactT1 -bga > ${sampleId}_ssDNA_type1_RPM-T1.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_ssDNA_type1_RPM-T1.bedgraph > ${sampleId}_ssDNA_type1_RPM-T1_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_ssDNA_type1_RPM-T1_sorted.bedgraph ${chrsize} \
            ${sampleId}_ssDNA_type1_RPM-T1.bw >& ${sampleId}_ssDNA_type1.bedGraphToBigWig.normfactT1.log 2>&1
        
        # Compute bedgraph and bigwig with normfactTot
        genomeCoverageBed -i ${ssDNA_type1_bed} -g ${chrsize} -scale \$normfactTot -bga > ${sampleId}_ssDNA_type1_RPM-Tot.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_ssDNA_type1_RPM-Tot.bedgraph > ${sampleId}_ssDNA_type1_RPM-Tot_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_ssDNA_type1_RPM-Tot_sorted.bedgraph ${chrsize} \
            ${sampleId}_ssDNA_type1_RPM-Tot.bw >& ${sampleId}_ssDNA_type1.bedGraphToBigWig.normfactTot.log 2>&1
     
        if [[ ${params.bigwig_profile} == "T12" || ${params.bigwig_profile} == "T12rep" ]]
        then  
            # Get normalization factors
            libsizeT12=`wc -l ${ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
            normfactT12=`python -c "print(round((1000000/(\$libsizeT12)),5))"`
    
            # Compute bedgraph and bigwig with normfactT12
            genomeCoverageBed -i ${ssDNA_type12_bed} -g ${chrsize} -scale \$normfactT12 -bga > ${sampleId}_ssDNA_type12_RPM-T12.bedgraph
            sort -k1,1 -k2,2n ${sampleId}_ssDNA_type12_RPM-T12.bedgraph > ${sampleId}_ssDNA_type12_RPM-T12_sorted.bedgraph
            bedGraphToBigWig ${sampleId}_ssDNA_type12_RPM-T12_sorted.bedgraph ${chrsize} \
                ${sampleId}_ssDNA_type12_RPM-T12.bw >& ${sampleId}_ssDNA_type12.bedGraphToBigWig.normfactT12.log 2>&1
     
            # Compute bedgraph and bigwig with normfactTot
            genomeCoverageBed -i ${ssDNA_type12_bed} -g ${chrsize} -scale \$normfactTot -bga > ${sampleId}_ssDNA_type12_RPM-Tot.bedgraph
            sort -k1,1 -k2,2n ${sampleId}_ssDNA_type12_RPM-Tot.bedgraph > ${sampleId}_ssDNA_type12_RPM-Tot_sorted.bedgraph
            bedGraphToBigWig ${sampleId}_ssDNA_type12_RPM-Tot_sorted.bedgraph ${chrsize} \
                ${sampleId}_ssDNA_type12_RPM-Tot.bw >& ${sampleId}_ssDNA_type12.bedGraphToBigWig.normfactTot.log 2>&1
        fi
        """
    }
     
    //IF RUNNING WITH REPLICATES AND BIGWIG PROFILES T1rep or T12rep
    if ( params.bigwig_profile == "T1rep" || params.bigwig_profile == "T12rep") {
    
            if ( params.nb_replicates == "2" ) {
    
            // Create channel to group samples by replicates
            // The resulting channel is composed of 7 elements : sampleID, norm factor file for replicate 1 (R1), 
            // norm factor file for replicate 2 (R2), T1 bed file for R1, T1 bed for file for R2,
            // T1+T2 bed file for R1, T1+T2 bed file for R2
            BEDtoBWrep
                .map { it -> [ it[0].split('_')[0..-3].join('_'), it[1], it[2], it[3] ] }
                .groupTuple(by: [0])
                .map { it -> it[0,1,2,3].flatten() }
                .set { BEDtoBWrep }
                //.println()
    
            // PROCESS 9 : makeBigwigReplicates (GENERATE NORMALIZED BIGWIG FILES FOR MERGED REPLICATES T1 AND T1+T2 BED FILES)
            // What it does : Merge T1 bed files between replicates, compute the normalization factors using total library sizes,
            // then use bedtools to generates bedgraph and bigwig files for T1 bed file, and merge T1+T2 files between replicates
            // and generate bigwig file if parameter bigwig_profile is set to T12rep
            // Input : channel BEDtoBW containing total libray size, T1 and T1+T2 bed files for the 2 replicates
            // Output : Bigwig files for merged T1 (and merged T1+T2 optional)
            process makeBigwigReplicates {
                tag "${sampleId}"
                label 'process_basic'
                label 'bigwig'
                publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",      	mode: params.publishdir_mode, pattern: "*.bw"
                publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",      	mode: params.publishdir_mode, pattern: "*.bedgraph"
                publishDir "${params.outdir}/01_logs/bigwig/${params.bigwig_profile}", 	mode: params.publishdir_mode, pattern: "*.log"
                input:
                    tuple val(sampleId), file(R1_libsizeTotal), file(R2_libsizeTotal), file(R1_ssDNA_type1_bed), file(R2_ssDNA_type1_bed), \
                        file(R1_ssDNA_type12_bed), file(R2_ssDNA_type12_bed) from BEDtoBWrep
                output:
                    file('*.bw')
    		    file('*.bedgraph')
                    file('*.log')
                    val('ok') into makeBigwigReplicates_ok
                script:
                """
                ## Compute normalization factor (i.e, for T1 : normfactor=1000000/librarySizeTotal)
                # Get total library sizes
                R1libsizeTot=`cat ${R1_libsizeTotal}`
                R2libsizeTot=`cat ${R2_libsizeTotal}`
                
                # Get T1 library sizes
                R1libsizeT1=`wc -l ${R1_ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
                R2libsizeT1=`wc -l ${R2_ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
                
                # Compute normalization factor
                R1R2normfactT1=`python -c "print(round(1000000/(\$R1libsizeT1+\$R2libsizeT1),5))"`
    	        R1R2normfactTot=`python -c "print(round(1000000/(\$R1libsizeTot+\$R2libsizeTot),5))"`
                
                ## Compute bedgraphs and bigwig for T1 and T12
                # Merge T1 bed files from the 2 replicates
                cat ${R1_ssDNA_type1_bed} ${R2_ssDNA_type1_bed} | grep -v "_random"| sort -k1,1 -k2,2n > ${sampleId}_R1R2_ssDNA_type1.bed
    
                # Compute bedgraph then bigwig for merged T1 bed files
                # normalized with R1R2normfactT1
    	        genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type1.bed -g ${params.chrsize} -scale \$R1R2normfactT1 \
                    -bga > ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bedgraph
                sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bedgraph > ${sampleId}_R1R2_ssDNA_type1_RPM-T1_sorted.bedgraph
                bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type1_RPM-T1_sorted.bedgraph \
                    ${params.chrsize} ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bw >& ${sampleId}_R1R2_ssDNA_type1.bedGraphToBigWig.normfactT1.log 2>&1
                
    	        # normalized with R1R2normfactoTot
        	genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type1.bed -g ${params.chrsize} -scale \$R1R2normfactTot \
                    -bga > ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bedgraph
                sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bedgraph > ${sampleId}_R1R2_ssDNA_type1_RPM-Tot_sorted.bedgraph
                bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type1_RPM-Tot_sorted.bedgraph \
                    ${params.chrsize} ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bw >& ${sampleId}_R1R2_ssDNA_type1.bedGraphToBigWig.normfactTot.log 2>&1
                
                # If bigwig_profile is "T12rep", compute the bigwigs for the merged T1 and T2 bed files
                if [[ ${params.bigwig_profile} == "T12rep" ]]
                then
                    # Get T1+T2 library sizes
                    R1libsizeT12=`wc -l ${R1_ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
                    R2libsizeT12=`wc -l ${R2_ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
                    
                    # Compute normalization factor
                    R1R2normfactT12=`python -c "print(round(1000000/(\$R1libsizeT12+\$R2libsizeT12),5))"`
                    
                    # Merge T1+T2 bed files from the 2 replicates
                    cat ${R1_ssDNA_type12_bed} ${R2_ssDNA_type12_bed} | grep -v "_random" | sort -k1,1 -k2,2n > ${sampleId}_R1R2_ssDNA_type12.bed
                
                    # Compute bedgraph then bigwig for merged T1+T2 bed files
                    # normalized with R1R2normfactT12
    		    genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type12.bed -g ${params.chrsize} -scale \$R1R2normfactT12 \
                        -bga > ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bedgraph
                    sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bedgraph > ${sampleId}_R1R2_ssDNA_type12_RPM-T12_sorted.bedgraph
                    bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type12_RPM-T12_sorted.bedgraph ${params.chrsize} \
                        ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bw >& ${sampleId}_R1R2_ssDNA_type12.bedGraphToBigWig.normfactT12.log 2>&1
                
    	            # normalized with R1R2normfactTot
                    genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type12.bed -g ${params.chrsize} -scale \$R1R2normfactTot \
                        -bga > ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bedgraph
                    sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bedgraph > ${sampleId}_R1R2_ssDNA_type12_RPM-Tot_sorted.bedgraph
                    bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type12_RPM-Tot_sorted.bedgraph ${params.chrsize} \
                        ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bw >& ${sampleId}_R1R2_ssDNA_type12.bedGraphToBigWig.normfactTot.log 2>&1
                
    	        fi
                """
            }
        }
        else {
            println("Selected bigwig_profile is ${params.bigwig_profile} but nb_replicates is ${params.nb_replicates}, check the parameters.")
            exit 0
        }
    }
    else {
        makeBigwigReplicates_ok = Channel.value( 'ok' )
    }        
                
    if ( params.kbrick_bigwig ) {
          
        // PROCESS 10 : makeDeeptoolsBigWig (GENERATES BIGWIG FILES, COVERAGE AND CUMULATIVE COVERAGE PLOTS)
        // What it does : uses deeptools to generate bigwig files for each of the 5 types of bam, then
        // plot the coverage for each bam files and plot the cumulative coverage (fingerprint)
        // Input : channel BAMwithIDXdt containing indexed bam files of the 5 types of bam
        // Output : bigwig files 
        process makeDeeptoolsBigWig { 
            tag "${sampleId}"
            label 'process_basic'
            label 'bigwig'
            publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}/plot",  	mode: params.publishdir_mode, pattern: "*.png"
            publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}",       	mode: params.publishdir_mode, pattern: "*.bigwig"
            publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}/tab",   	mode: params.publishdir_mode, pattern: "*.tab"
            publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}",mode: params.publishdir_mode, pattern: "*log*"
        input:
            set val(sampleId), file(bam), file(bamidx) from BAMwithIDXdt
        output:
            file('*')
            val 'ok' into kb_bigwig_ok
            path('*.log')
        shell:
        """
        # Use deeptools bamCoverage to generate bigwig file with RPKM normalization
        bamCoverage --bam ${bam} --normalizeUsing RPKM --binSize ${params.binsize} \
            --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.RPKM.bigwig >& ${sampleId}_makeDeeptoolsBigWig_bamCoverage.log 2>&1
        # Use deeptools plotCoverage to plot the coverage plot for each bam file (to assess the sequencing depth) 
        plotCoverage --bamfiles ${bam} --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.coveragePlot.png >& ${sampleId}_makeDeeptoolsBigWig_plotCoverage.log 2>&1
        # Use deeptools plotFingerprint to plot a profile of cumulative read coverages (quality control to assess ChIP signal from background)
        plotFingerprint --bamfiles ${bam} --labels ${bam} --numberOfProcessors ${task.cpus} \
            --minMappingQuality 30 --skipZeros --plotFile ${bam.baseName}.deeptools.fingerprints.png \
            --outRawCounts ${bam.baseName}.deeptools.fingerprints.tab >& ${sampleId}_makeDeeptoolsBigWig_plotFingerprint.log 2>&1
        """ 
        }
    
        BAMwithIDXfr
            .combine(ssDNA_to_bigwigs_FASTER_LOMEM_ch)
            .set { BAMwithIDXfr_ssDNA_to_bigwigs_ch }
    
        fasta_ch3
            .combine(index_fai_ch2)
            .set { fasta_fai_ch }
    
        BAMwithIDXfr_ssDNA_to_bigwigs_ch
            .combine(fasta_fai_ch)
            .set {  BAMwithIDXfr_ssDNA_to_bigwigs_fasta_ch }
    
        // PROCESS 11 : toFRBigWig (GENERATES FWD/REV BIGWIG FILES)
        // What it does : Generates forward/reverse bigwig files for each of the 5 types of bam
        // Input : channel BAMwithIDXfr containing indexed bam files of the 5 types of bam
        // Output : FR bigwig files
        // External tool : Perl script from K. Brick (original pipeline, 2012) 
        process toFRBigWig {
            tag "${sampleId}"
            label 'process_basic'
            label 'bigwig'
            publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/FRbigwig", 	mode: params.publishdir_mode, pattern: '*.out'
            publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/FRbigwig",  	mode: params.publishdir_mode, pattern: '*.log'
            publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/FRbigwig",      	mode: params.publishdir_mode, pattern: '*.bigwig'
            input:
                tuple val(sampleId), path(bam), path(bamidx), path(ssDNA_to_bigwigs_FASTER_LOMEM_script), path(fasta), path(fai) from BAMwithIDXfr_ssDNA_to_bigwigs_fasta_ch
            output:
                file('*')
                file('*.log')
    	    val 'ok' into fr_bigwig_ok
            script:
            """
            # Use K. Brick perl script to generate forward/reverse bigwig files
            # (Normalization is done relatively to input bam file size and not total library size)
            perl ${ssDNA_to_bigwigs_FASTER_LOMEM_script} --bam ${bam} --g ${fasta} --o ${bam.baseName}.out \
                --s 100 --w 1000 --sc "." --gIdx ${fai} -v >& ${sampleId}_toFRBigwig.log 2>&1
            """
        }    
    }
    else {
        kb_bigwig_ok = Channel.value( 'ok' )
        fr_bigwig_ok = Channel.value( 'ok' )
    }
    
    
    //***************************************************************************//
    //                          SECTION 5 : PEAK CALLING                         //
    //***************************************************************************//
    
    // Set saturation curve thresholds for callPeaks and makeSatCurve processes
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
        // CREATE CHANNEL LINKING IP TYPE 1 SSDNA BED WITH CONTROL BED
        // ( The control bed is the input control filtered bam file )
        // The resulting channel will be  composed of 4 elements : sampleID, controlID, chIP type1 bed file, input filtered bam file
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
    
        ch_design_controls_csv
            .combine(T1BED)
            .combine(INPUTBED)
            .filter { it[0] == it[5] && it[1] == it[7] }
            .map { it -> it[0,1,6,8].flatten() }
            .into { T1BED_shuffle_ch ; T1BED_replicate_ch }
            //.println()
    
        // PROCESS 12 : shufBEDs (BED SHUFFLING)
        // What it does : quality trims and shuffles type1 bed files from ITR parsing 
        // Input : type1 bed files from ITR parsing
        // Ouptut : shuffled and filtered type1 bed files
        process shufBEDs_ct {
            tag "${id_ip}"
            label 'process_medium'
            label 'peakcalling'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/bed_shuffle/trim_q${params.bed_trimqual}",  mode: params.publishdir_mode
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
                perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) \
                    if (\$Q[0] >= ${params.bed_trimqual} && \$Q[1] >= ${params.bed_trimqual})' ${ip_bed} >${ip_Q30_bed}
                perl -lane '@F = split(/\\t/,\$_); print join("\\t",@F) if (\$F[4] >= ${params.bed_trimqual} )' ${ct_bed} >${ct_Q30_bed}
                # Bed shuffling
                if [ ${params.genome_name} == "mm10" ]
                then
                    shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
                    shuf ${ct_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ct_Q30_shuf_bed}
                else
                    shuf ${ip_Q30_bed} >${ip_Q30_shuf_bed}
                    shuf ${ct_Q30_bed} >${ct_Q30_shuf_bed}
                fi
                """
        }
    
        // CREATE CHANNEL TO COMBINE SQ30BED AND SHUFFLE PERCENT CHANNEL
        SQ30BED_ch
            .combine(satCurvePCs)
            .set { SQ30BED_satcurve_ch }
    
        SQ30BED_satcurve_ch
            .combine(blacklist_ch)
            .set {SQ30BED_satcurve_blk_ch }
    
        index_fai_ch3
            .combine(pickNlines_ch)
            .set { fai_pickNlines }
    
        SQ30BED_satcurve_blk_ch
            .combine(fai_pickNlines)
            .set { SQ30BED_satcurve_blk_fai_pickNlines_ch }
    
        // PROCESS 13 : callPeaks (PEAK CALLING WITH MACS2)
        // What it does : call peaks in type1 bed files using macs2. if parameter satCurveReps is > 1,
        // the process will repeat the peak calling satCurveReps-1 times then the resulting bed files will be merged.
        // if satCurve parameter is true, the process will progressively call peaks in downsampled bed files (according to satCurvePCs parameter)
        // to evaluate the saturation of the samples.
        // Input : Shuffled type1 bed files ; Shuffled controlled bed files  and satCurvePCs percent channel
        // Output : Bed/xls files for peaks called in downsampled bed file; and bed/xls files for peaks called in the whole input bed file.
        // Also outputs genome size parameters (wil be used in post processes) 
        process callPeaks_ct {
            tag "${id_ip}"
            label 'process_basic'
            label 'peakcalling'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/saturation_curve/${params.sctype}", mode: params.publishdir_mode, pattern: "*peaks_sc.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/bed",         mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks_sc.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/bed",         mode: params.publishdir_mode, pattern: "*summits.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/xls",         mode: params.publishdir_mode, pattern: "*peaks*.xls"
            publishDir "${params.outdir}/01_logs/peaks/${control_status}/macs2/${macs2_params}",     	    mode: params.publishdir_mode, pattern: "*.macs2.log"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks",                        mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks",                        mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.xls"
    	publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/narowpeaks",  mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks.narrowPeak"
            input:
                tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed), val(shuffle_percent), path(blacklist), path(fai), path(pickNlines_script) from SQ30BED_satcurve_blk_fai_pickNlines_ch
            output:
                path("*peaks*.bed") into allbed
                path("*peaks_sc.bed") into satcurvebed_ch optional true
                path("*peaks*.xls") optional true
                path("*peaks.bedgraph") optional true
                path("*.macs2.log")
                tuple val(id_ip), file(ip_bed), file("*peaks.bed") optional true into ALLPEAKSTOPP, ALLPEAKSTOPP2
                stdout into gsize
                val 'ok' into callPeaks_ok
    	    path("*")
                tuple val(id_ip), path("*finalPeaks_noIDR.bed") optional true into FINALPEAKSNOIDR
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
                genome_size=`cut -f3 ${fai} |tail -n1`
            else
                nb_line_bl=`cat ${blacklist} | wc -l`
                tot_sz=`cut -f3 ${fai} |tail -n1`
                if [ \$nb_line_bl == "0" ]
                then
                    genome_size=`expr \$tot_sz - 0`
                else
                    bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${blacklist} |tail -n1`
                    genome_size=`expr \$tot_sz - \$bl_size`
                fi
            fi
    
            # Export the genome_size value for post processes
            echo -n \$genome_size
    
            ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter 
            for i in {0..${satCurveReps}}; do
                # Create a unique name for output file
                peakFileName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
                # Call peaks with macs2
    	    if [ ${params.macs_pv} != -1 ]; then
    	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                        --name \$peakFileName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} > IP-${id_ip}_CT-${ct_bed}.macs2.log 2>&1
            else
                macs2 callpeak -g \$genome_size -t \$nPC.IP.bed -c ${ct_bed} \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                        --name \$peakFileName --nomodel --extsize ${params.macs_extsize}  --qvalue ${params.macs_qv} > IP-${id_ip}_CT-${ct_bed}.macs2.log 2>&1
            fi
        
            # Filter out peaks from blacklist bed file
            if [ ${params.blacklist} != "None" ]
            then
                intersectBed -a \$peakFileName'_peaks.narrowPeak' -b ${blacklist} -v >\$peakFileName'.peaks_sc.noBL'
                # Filter out mitochondrial peaks
                cut -f1-3 \$peakFileName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
            else
                cut -f1-3 \$peakFileName'_peaks.narrowPeak' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
            fi
            # Filter out chrY peaks if params.no_chrY is true
            if [ ${params.no_chrY} ]
            then
                grep -v 'chrY' \$peakFileName'_peaks_sc.noM' |sort -k1,1 -k2n,2n > \$peakFileName'_peaks_sc.bed'
            else
                mv \$peakFileName'_peaks_sc.noM' \$peakFileName'_peaks_sc.bed'
            fi
     
            # Rename excel file
            mv \$peakFileName'_peaks.xls' \$peakFileName'_peaks_sc.xls'
            done
    
            ##MERGE BED FILES FROM THE LOOP
            sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
            
            ##POSTPROCESS IF shuffle_percent is 100% ie if all file has been treated 
            if [ ${shuffle_percent} == 1.00 ]; then
                mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
                cp ${id_ip}.peaks.bed ${id_ip}.finalPeaks_noIDR.bed
                cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
                cp ${id_ip}.peaks.xls ${id_ip}.finalPeaks_noIDR.xls
            fi
            rm *.tmp
            """
        }
    }
    
    // CASE 2 : IF NO INPUT CONTROL IS PROVIDED
    else {
    
       T1BED
            .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
            .groupTuple(by: [0])
            .map { it ->  [ it[0], it[1].flatten() ] }
            .set { T1BED }
    
        // PROCESS 12 : shufBEDs (BED SHUFFLING)
        // What it does : quality trims and shuffles type1 bed files from ITR parsing 
        // Input : type1 bed files from ITR parsing
        // Ouptut : shuffled and filtered type1 bed files
        process shufBEDs {
            tag "${id_ip}"
            label 'process_medium'
            label 'peakcalling'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/bed_shuffle/trim_q${params.bed_trimqual}",  mode: params.publishdir_mode
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
                perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) \
                        if (\$Q[0] >= ${params.bed_trimqual} && \$Q[1] >= ${params.bed_trimqual})' ${ip_bed} >${ip_Q30_bed}
                # Bed shuffling
                if [ ${params.genome_name} == "mm10" ]
                then
                    shuf ${ip_Q30_bed} |grep -P '^chr[0123456789IVLXYZW]+\\s' >${ip_Q30_shuf_bed}
                else
                    shuf ${ip_Q30_bed}>${ip_Q30_shuf_bed}
                fi
                """
        }
    
        // CREATE CHANNEL TO COMBINE SQ30BED AND SHUFFLE PERCENT CHANNEL
        SQ30BED_ch
            .combine(satCurvePCs)
            .set { SQ30BED_satcurve_ch }
    
        SQ30BED_satcurve_ch
            .combine(blacklist_ch)
            .set {SQ30BED_satcurve_blk_ch }
    
        index_fai_ch3
            .combine(pickNlines_ch)
            .set { fai_pickNlines }
    
        SQ30BED_satcurve_blk_ch
            .combine(fai_pickNlines)
            .set { SQ30BED_satcurve_blk_fai_pickNlines_ch }
    
    
        // PROCESS 13 : callPeaks (PEAK CALLING WITH MACS2)
        // What it does : call peaks in type1 bed files using macs2. if parameter satCurveReps is > 1,
        // the process will repeat the peak calling satCurveReps-1 times then the resulting bed files will be merged.
        // if satCurve parameter is true, 
        // the process will progressively call peaks in downsampled bed files (according to satCurvePCs parameter)
        // to evaluate the saturation of the samples.
        // Input : Shuffled type1 bed files and satCurvePCs percent channel
        // Output : Bed/xls files for peaks called in downsampled bed file; and bed/xls files for peaks called in the whole input bed file. 
        // Also outputs genome size parameters (wil be used in post processes)
        process callPeaks {
            tag "${id_ip}"
            label 'process_medium'
            label 'peakcalling'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/saturation_curve/${params.sctype}/peaks", mode: params.publishdir_mode, pattern: "*peaks_sc.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/bed",               mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks_sc.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/bed",               mode: params.publishdir_mode, pattern: "*.summits.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/xls",               mode: params.publishdir_mode, pattern: "*.peaks*.xls"
            publishDir "${params.outdir}/01_logs/peaks/${control_status}/macs2/${macs2_params}",	                  mode: params.publishdir_mode, pattern: "*.macs2.log"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks",                              mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.bed"
            publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks",                              mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.xls"
    	publishDir "${params.outdir}/02_results/peaks/${control_status}/macs2/${macs2_params}/narrowPeak",	  mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks.narrowPeak"
            input:
                tuple val(id_ip), path(ip_bed), val(shuffle_percent), path(blacklist), path(fai), path(pickNlines_script) from SQ30BED_satcurve_blk_fai_pickNlines_ch
            output:
                path("*peaks*.bed") into allbed
                path("*peaks_sc.bed") into satcurvebed_ch optional true
                path("*peaks*.xls") optional true
                path("*peaks*.bedgraph") optional true
                path("*.macs2.log")
    	    path("*")
                tuple val(id_ip), file(ip_bed), file("*peaks.bed") optional true into ALLPEAKSTOPP, ALLPEAKSTOPP2
                stdout into gsize
                val 'ok' into callPeaks_ok
                tuple val(id_ip), path("*finalPeaks_noIDR.bed") optional true into FINALPEAKSNOIDR
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
                genome_size=`cut -f3 ${fai} |tail -n1`
            else
    	    nb_line_bl=`cat ${blacklist} | wc -l`
                tot_sz=`cut -f3 ${fai} |tail -n1`
    	    if [ \$nb_line_bl == "0" ]
    	    then
    		genome_size=`expr \$tot_sz - 0`
    	    else
                	bl_size=`perl -lane '\$tot+=(\$F[2]-\$F[1]); print \$tot' ${blacklist} |tail -n1`
                	genome_size=`expr \$tot_sz - \$bl_size`
    	    fi
            fi
    
            # Export the genome_size value for post processes
            echo -n \$genome_size
        
            ## CALL PEAKS WITH MACS2 N TIMES ACCORDING TO params.rep parameter
            for i in {0..${satCurveReps}}; do
                # Create a unique name for output file
                peakFileName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
                # Call peaks with macs2
                if [ ${params.macs_pv} != -1 ]; then
    	        macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                        --name \$peakFileName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} > IP-${id_ip}.macs2.log 2>&1
                else
                    macs2 callpeak -g \$genome_size -t \$nPC.IP.bed \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                        --name \$peakFileName --nomodel --extsize ${params.macs_extsize} --qvalue ${params.macs_qv} > IP-${id_ip}.macs2.log 2>&1
                fi
    
                # Filter out peaks from blacklist bed file  
                if [ ${params.blacklist} != "None" ]
                then
                    intersectBed -a \$peakFileName'_peaks.narrowPeak' -b ${blacklist} -v >\$peakFileName'.peaks_sc.noBL'
                    # Filter out mitochondrial peaks
                    cut -f1-3 \$peakFileName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
                else
                    cut -f1-3 \$peakFileName'_peaks.narrowPeak' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
                fi
                # Filter out chrY peaks if params.no_chrY is true
                if [ ${params.no_chrY} ]
                then
                    grep -v 'chrY' \$peakFileName'_peaks_sc.noM' |sort -k1,1 -k2n,2n > \$peakFileName'_peaks_sc.bed'
                else
                    mv \$peakFileName'_peaks_sc.noM' \$peakFileName'_peaks_sc.bed'
                fi
     
                # Rename excel file
                mv \$peakFileName'_peaks.xls' \$peakFileName'_peaks_sc.xls'
            done
     
            ##MERGE BED FILES FROM THE LOOP
            sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed
      
            ##POSTPROCESS IF shuffle_percent is 100% ie if all file has been treated
            if [ ${shuffle_percent} == 1.00 ]; then
                mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
                cp ${id_ip}.peaks.bed ${id_ip}.finalPeaks_noIDR.bed
                
                cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
                cp ${id_ip}.peaks.xls ${id_ip}.finalPeaks_noIDR.xls
            fi
            rm *.tmp
            """
        }
    }
    
    //***************************************************************************//
    //                           SECTION 6 : SSDS QC                             //
    //***************************************************************************//
    
    // PROCESS 14 : samStats (GENERATES SAMSTATS REPORTS)
    // What it does : Use Samtools to report comprehensive statistics from alignment file
    // Input : channel BAMwithIDXss containing indexed bam files of the 5 types of bam
    // Ouptut : Comprehensive and summary of statistical reports for bam files
    process samStats {
        tag "${sampleId}"
        label 'process_basic'
        label 'ssds'
        publishDir "${params.outdir}/02_results/qc/samstats",  mode: params.publishdir_mode, pattern: "*.tab"
        input:
            set val(sampleId), file(bam), file(bamidx) from BAMwithIDXss
        output:
            file '*stats.tab' into samstats_ch
            val 'ok' into samstats_ok
        script:
            """
            samtools idxstats ${bam} > ${bam.baseName}.idxstats.tab
            samtools stats ${bam} > ${bam.baseName}.samstats.tab
            """
    }
    
    // CREATE CHANNEL FOR makeSSreport PROCESS
    ITRBED
        .combine(hotspots_ch)
        .set { ITRBED_hotspots_ch }
    
    ITRBED_hotspots_ch
        .combine(makeSSMultiQCReport_nextFlow_ch)
        .set { ITRBED_hotspots_makeSSMultiQCReport_ch }
    
    // PROCESS 15 : makeSSreport (COMPUTE STATS FROM SSDS PARSING)
    // What it does : Compile alignement stats from the 5 types of bed for a given sample
    // and compute FRIP scores (the fraction of reads that fall into a peak)
    // Input : Total bam file and parsed bed files (T1, T2, ds, ds_strict and unclassified)
    // Output : text report
    // External tool : Perl script from K. Brick (original pipeline, 2012)
    process makeSSreport {
        tag "${sampleId}"
        label 'process_basic'
        label 'ssds'
        publishDir "${params.outdir}/02_results/qc/multiqc",  mode: params.publishdir_mode, pattern: '*SSDSreport*'
        publishDir "${params.outdir}/01_logs/makeSSreport", mode: params.publishdir_mode, pattern: '*.log'
        input:
            tuple val(sampleId), path(bam), path(T1), path(T2), path(ds), path(dss), path(unclassified), path(hotspots), path(makeSSMultiQCReport_nextFlow_script) from ITRBED_hotspots_makeSSMultiQCReport_ch
        output:
            set val(sampleId), file('*SSDSreport*') into SSDSreport2ssdsmultiqc, SSDSreport2ssdsmultiqc2
            path('*SSDSreport*') into ssdsreport_ch
            val 'ok' into ssreport_ok
            file('*.log')
        script:
        """
        perl ${makeSSMultiQCReport_nextFlow_script} ${bam} ${T1} ${T2} ${ds} ${dss} ${unclassified} \
            --g ${params.genome} --h ${hotspots} >& ${sampleId}_makeSSreport.log 2>&1
        """
    }
    
    
    // PROCESS 16 : makeFingerPrint (MAKE DEEPTOOLS FINGERPRINT PLOTS)
    // What it does : plot a profile of cumulative read coverages (quality control to assess ChIP signal) for all samples
    // Input : filtered & indexed bam files for all samples (control and IP)
    // Output : Plot and tab files for metrics 
    
    allbam_ch = UNPARSEDBAM.collect()
    allbai_ch = UNPARSEDBAI.collect()
    
    process makeFingerPrint {
        tag "${global_name}"
        label 'process_low'
        label 'bigwig'
        publishDir "${params.outdir}/02_results/qc/fingerprint",  mode: params.publishdir_mode, pattern: '*.png'
        publishDir "${params.outdir}/02_results/qc/fingerprint",  mode: params.publishdir_mode, pattern: '*.tab'
        publishDir "${params.outdir}/01_logs/fingerprint",  mode: params.publishdir_mode, pattern: '*.log'
        input:
            val('filterbam_ok') from filterbam_tofingerprint.collect()
            path(all_bam) from allbam_ch
            path(all_bai) from allbai_ch
        output:
            file('*.png') into png_fingerprint_ch
            file('*.tab') into tab_fingerprint_ch
            val('ok') into makefingerprint_ok
            path('*.log')
        script:
        """
        #Use deeptools plotFingerprint to plot a profile of cumulative read coverages (quality control to assess ChIP signal from background)
        plotFingerprint --bamfiles ${all_bam} --smartLabels --numberOfProcessors ${task.cpus} \
            --minMappingQuality 30 --skipZeros --plotFile ${global_name}.deeptools.fingerprints.png \
            --outRawCounts ${global_name}.deeptools.fingerprints.rawcounts.tab --outQualityMetrics ${global_name}.deeptools.fingerprints.qual.tab \
    	>& ${global_name}_makeFingerPrint.log
        """
    }
    
    
    // PROCESS 17 : computeFRIP (COMPUTE FRIP SCORE FOR PARSED BAM FILE FOR NEW GENOMES)
    // What it does : Use deeptools to compute DE NOVO FRIP scores (the fraction of reads from bam that fall into a peak)
    // Input : Parsed bed files (T1, T2, ds, ds_strict and unclassified); final but not merged peak bed file and SSDS QC report from process 13
    // Output : Extended QC report 
    // External tool : Python script (Pauline Auffret, 2021) 
    if (params.hotspots == "None") {
    
        // First, create channel to combine bam files and peak file and SSDS QC report
        ITRBAM
            .combine(FINALPEAKSNOIDR)
            .combine(SSDSreport2ssdsmultiqc2)
            .map { it -> [ it[0].split('_')[0..-2].join('_'),it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[10], it[11], it[12], it[13], it[14], it[15].split('_')[0..-2].join('_'), it[16] ] }
            .filter { it[0] == it[13] && it[0] == it[15] }
            .map { it -> it[0,1,2,3,4,5,6,7,8,9,10,11,12,14,16].flatten() }
            .set { ITRBAMANDPEAKS }
    	//.into { ITRBAMANDPEAKS ; test2 }
    
        ITRBAMANDPEAKS
            .combine(get_frip_ch)
            .set { ITRBAMANDPEAKS_get_frip_ch }
    
        process computeFRIP {
            tag "${sampleId}"
            label 'process_low'
            label 'frip'
            publishDir "${params.outdir}/02_results/qc/multiqc",	mode: params.publishdir_mode, patten: '*SSDSreport*'
    	publishDir "${params.outdir}/02_results/qc/frip",	mode: params.publishdir_mode, patten: '*'
            publishDir "${params.outdir}/01_logs/frip",		mode: params.publishdir_mode, patten: '*.log'
            input:
                tuple val(sampleId), path(bam), path(bai), path(T1), path(T1bai), path(T2), path(T2bai), path(ds), path(dsbai), \
                    path(dss), path(dssbai), path(unc), path(uncbai), path(peaks), path(report), path(get_frip_script) from ITRBAMANDPEAKS_get_frip_ch
            output:
                set val(sampleId), file('*SSDSreport*') into SSDSreport2ssdsmultiqcdenovo
                path('*SSDSreport*') into ssds_report_denovo_ch
    	    path("*")
    	    val('ok') into computefrip_ok
            script:
            """
            # For each bam type, execute python script to compute FRIP score with Deeptools and print results to text file
            python ${get_frip_script} ${peaks} ${bam} ${task.cpus} total ${params.genome_name} ${sampleId}_total.frip #Need to edit custom ssds multiQC library to include this in SSDS multiqc  #todo 
            python ${get_frip_script} ${peaks} ${T1} ${task.cpus} ssType1 ${params.genome_name} ${sampleId}_type1.frip >& ${sampleId}_computeFRIP_type1.log 2>&1
            python ${get_frip_script} ${peaks} ${T2} ${task.cpus} ssType2 ${params.genome_name} ${sampleId}_type2.frip >& ${sampleId}_computeFRIP_type2.log 2>&1
            python ${get_frip_script} ${peaks} ${ds} ${task.cpus} dsDNA ${params.genome_name} ${sampleId}_ds.frip >& ${sampleId}_computeFRIP_ds.log 2>&1
            #python ${get_frip_script} ${peaks} ${dss} ${task.cpus} dsDNA_strict ${params.genome_name} ${sampleId}_dss.frip #Need to edit ssds multiQC library to include this in SSDS multiqc  #todo
    	    python ${get_frip_script} ${peaks} ${dss} ${task.cpus} dsDNA ${params.genome_name} ${sampleId}_dss.frip >& ${sampleId}_computeFRIP_dss.log 2>&1
            python ${get_frip_script} ${peaks} ${unc} ${task.cpus} unclassified ${params.genome_name} ${sampleId}_unc.frip >& ${sampleId}_computeFRIP_unc.log 2>&1
    
            # Concatenate results FRIP text files with SSDS QC report from process 13
            #cat ${report} ${sampleId}_total.frip ${sampleId}_type1.frip ${sampleId}_type2.frip ${sampleId}_ds.frip ${sampleId}_dss.frip ${sampleId}_unc.frip > ${sampleId}_denovo.SSDSreport.tab 
            cat ${report} ${sampleId}_type1.frip ${sampleId}_type2.frip ${sampleId}_ds.frip ${sampleId}_dss.frip ${sampleId}_unc.frip > ${sampleId}_denovo.SSDSreport.tab
    	"""
        }
    
    
        if (params.with_ssds_multiqc) {
    
            runPlotSSDSqc_ch
                .combine(plot_ssds_stat_ch)
                .set { ssds_qc_scripts_ch }
    
            SSDSreport2ssdsmultiqcdenovo
                .combine(ssds_qc_scripts_ch)
                .set { SSDSreport2ssdsmultiqcdenovo_script_ch }
    
            // PROCESS 18 : ssds_multiqc_denovo (MAKE MULTIQC REPORT FOR SSDS FILES)
            // What it does : For each sample, wompute a multiqc quality control report
            // Input : QC report from SSDSreport2ssdsmultiqcdenovo
            // Ouptut : multiQC HTML report
            // External tool : custom multiqc python library and conda environment
            process ssds_multiqc_denovo {
                tag "${sampleId}"
        	    label 'process_basic'
                label 'plot'
                publishDir "${params.outdir}/02_results/qc/multiqc",  mode: params.publishdir_mode
                publishDir "${params.outdir}/01_logs/multiqc",  mode: params.publishdir_mode, pattern: '*.log'
                input:
                    tuple val(sampleId), path(report), path(runPlotSSDSqc_script), path(plot_ssds_stat_script) from SSDSreport2ssdsmultiqcdenovo_script_ch
                output:
                    file('ssds/*.png') into multiqc_denovo_ch
    		path('ssds/*.pdf')
    		val('ok') into multiqc_denovo_ok
                script:
                """
                #multiqc -m ssds -n ${sampleId}.multiQC . 
                mkdir -p ssds
                Rscript ${runPlotSSDSqc_script} ${plot_ssds_stat_script} ${report} ${sampleId} "./ssds" >& ${sampleId}_ssds_multiqc_denovo.log
                cd ssds ; for file in \$(ls *.png); do mv \$file \${file%.*}_mqc.png; done ; cd ..
                """
            }
        }
    }
    
    if (!params.with_ssds_multiqc)  {
        computefrip_ok = Channel.value( 'ok' )
        multiqc_denovo_ok = Channel.value( 'ok' )
        ssds_report_denovo_ch = Channel.fromPath(emptyfile1)
        multiqc_denovo_ch = Channel.fromPath(emptyfile2)
    }
    
    
    if (params.with_ssds_multiqc && params.hotspots != "None") {
    
        runPlotSSDSqc_ch
            .combine(plot_ssds_stat_ch)
            .set { ssds_qc_scripts_ch }
    
        SSDSreport2ssdsmultiqc
            .combine(ssds_qc_scripts_ch)
            .set { SSDSreport2ssdsmultiqc_script_ch }
    
        // PROCESS 18 : ssds_multiqc (MAKE MULTIQC REPORT FOR SSDS FILES)
        // What it does : For each sample, wompute a multiqc quality control report
        // Input : QC report from SSDSreport2ssdsmultiqc
        // Ouptut : multiQC HTML report
        // External tool : custom multiqc python library and conda environment
        process ssds_multiqc {
            tag "${sampleId}"
            label 'process_basic'
            label 'plot'
            publishDir "${params.outdir}/02_results/qc/multiqc",		mode: params.publishdir_mode
            publishDir "${params.outdir}/01_logs/multiqc",			mode: params.publishdir_mode, pattern: '*.log'
            input:
                tuple val(sampleId), path(report), path(runPlotSSDSqc_script), path(plot_ssds_stat_script) from SSDSreport2ssdsmultiqc_script_ch
            output:
                file('ssds/*.png') into multiqc_ch
                path('ssds/*.pdf')
    	    val('ok') into ssds_multiqc_ok
                val('ok') into computefrip_ok
                val('ok') into multiqc_denovo_ok
            script:
            """
            #multiqc -m ssds -n ${sampleId}.multiQC .
            mkdir -p ssds
            Rscript ${runPlotSSDSqc_script} ${plot_ssds_stat_script} ${report} ${sampleId} "./ssds" >& ${sampleId}_ssds_multiqc.log
            cd ssds ; for file in \$(ls *.png); do mv \$file \${file%.*}_mqc.png; done ; cd ..
            """
            }
        multiqc_denovo_ok = Channel.value( 'ok' )
        ssds_report_denovo_ch = Channel.fromPath(emptyfile3)
        multiqc_denovo_ch = Channel.fromPath(emptyfile4)
    }
    else {
        ssds_multiqc_ok = Channel.value( 'ok' )
        multiqc_denovo_ok = Channel.value( 'ok' )
        ssds_report_denovo_ch = Channel.fromPath(emptyfile3)
        multiqc_denovo_ch = Channel.fromPath(emptyfile4)
        multiqc_ch = Channel.fromPath(emptyfile5)
    }
    
    //***************************************************************************//
    //                     SECTION 7 : OPTIONAL IDR ANALYSIS                     //
    //***************************************************************************//
    // THIS ONLY WORKS WITH 2 REPLICATES IN THIS VERSION. #todo
    
    if (params.with_idr && params.nb_replicates == "2" ) {
        // CASE 1 : IF INPUT CONTROL ARE PROVIDED
        if (params.with_control) {
     
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
    
            // CREATE CHANNEL TO GROUP SAMPLE ID & CONTROL ID WITH ASSOCIATED REPLICATES (BED FILES)
            // The resulting channel is composed of 6 elements : sampleID, controlID, T1BED_R1, T1BED_R2, DSBED_R1, DSBED_R2
            T1BED_replicate_ch_renamed
                .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1].split('_')[0..-2].join('_'), it[2], it[3] ] }
                .groupTuple(by: [0])
                .groupTuple(by: [1])
                .map { it ->  [ it[0], it[1][0], it[2], it[3] ] }
                .map { it -> it[0,1,2,3].flatten() }
                .set { T1BED_replicate_ch_renamed }
                //.println()
    
    
            // PROCESS 19 : createPseudoReplicates (CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
            // What it does : creates 2 pseudo replicates per true replicates, then pool the true replicates, and creates 2 pseudo replicates from this pool.
            // Input : bed files from type1 aligned SSDNA, chip and control
            // Output : The 2 true replicates ; the pool of the 2 true replicates ; the 4 pseudo replicates from true replicates; 
            // the 2 pseudo replicates from the pool of true replicates, and 1 control file (merge of control files if they are different)
            // So 10 files in total.
            process createPseudoReplicates_ct {
                tag "${id_ip}"
                label 'process_medium'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/pseudo_replicates",  mode: params.publishdir_mode, pattern: '*.bed'
                input:
                    tuple val(id_ip), val(id_ct), file(t1_rep1), file(t1_rep2), file(ct_r1), file(ct_r2) from T1BED_replicate_ch_renamed
                output:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed'), \
                        file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed'), file('*pool_r1.bed'), file('*pool_r2.bed'), \
                        file('*poolT.bed'), file('*ct_pool.bed') into ALLREP
                    val 'ok' into createPseudoReplicates_ok
                script:
                """
                # Shuffle then split the 2 original replicates in 2 pseudo replicates
                nlines_r1=\$((`cat ${t1_rep1} | wc -l`/2)) 
                nlines_r2=\$((`cat ${t1_rep2} | wc -l`/2))
                shuf ${t1_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
                shuf ${t1_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
                mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
                mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
                mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
                mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed
    
                # Pool the 2 original replicates then shuffle then split in 2 pseudo replicates
                nlines_pool=\$((`cat ${t1_rep1} ${t1_rep2} | wc -l`/2))
                cat ${t1_rep1} ${t1_rep2} > ${id_ip}_poolT.bed
                shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
                mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
                mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
    
                # Test if input files are the same ; if not, merge them to build a new control file for IDR
                if cmp -s ${ct_r1} ${ct_r2} 
                then 
                    cat ${ct_r1} > ${id_ct}_ct_pool.bed
                else
                    # check the merging process is good #todo
                    cat ${ct_r1} ${ct_r2} | sort -n | unique > ${id_ct}_ct_pool.bed
                fi
                """
            }
    
        
            // PROCESS 20 : callPeaksForIDR (CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
            // What it does : uses macs2 callPeak to perform peak-calling on all replicates (2), pseudo replicates (4), pool (1), pool pseudo replicates (2)
            // Input : all 10 bed files from true replicates and pseudo replicates creation; and genome size from process 14 
            // Output : 9 regionPeak files from peak calling
            process callPeaksForIDR_ct {
                tag "${id_ip}"
                label 'process_medium'
                label 'peakcalling'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2",      mode: params.publishdir_mode, pattern: '*narrowPeak*'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2",      mode: params.publishdir_mode, pattern: '*regionPeak*'
                publishDir "${params.outdir}/01_logs/idr/${control_status}/${idr_params}/macs2",  mode: params.publishdir_mode, pattern: '*.macs2.log'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2",      mode: params.publishdir_mode, pattern: '*'
                input:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2), file(r1_pseudorep_r1), file(r1_pseudorep_r2),\
                        file(r2_pseudorep_r1), file(r2_pseudorep_r2), file(pool_r1), file(pool_r2), \
                        file(poolT), file(ctpool) from ALLREP
                    val(genome_size) from gsize
                output:
                    file('*narrowPeak*')
                    tuple val(id_ip),file(t1_rep1), file(t1_rep2), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak'),\
                        file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak'),\
                        file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak'), \
                        file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak'), file('*poolT*.regionPeak') into ALLPEAKSREP
                    file('*.macs2.log')
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
                        mkdir \${random_id}
                        if [[ ${params.idr_macs_pv} == -1 && ${params.idr_macs_qv} == -1 ]];
                        then
                            macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                            --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                            --name \$name --nomodel --extsize ${params.macs_extsize} > \$name.macs2.log 2>&1
                        elif [ ${params.idr_macs_pv} != -1 ];
                        then
                            macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                            --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                            --name \$name --nomodel --extsize ${params.macs_extsize} --pvalue ${params.idr_macs_pv} > \$name.macs2.log 2>&1
                        elif [ ${params.idr_macs_qv} != -1 ];
                        then
                            macs2 callpeak -g ${genome_size} -t \$file -c ${ctpool} \
                            --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                            --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.idr_macs_qv} > \$name.macs2.log 2>&1
                        fi
    
                        # Cut resulting bed files to 10000 lines max for better handling in IDR statistical analysis
                        npeaks=`cat \${name}_peaks.narrowPeak | wc -l`
                        # Check if the number of peaks is not too low
                        if [ \$npeaks -lt 20 ]
                        then
                            echo "Not enough peaks in \${name}_peaks.narrowPeak file to perform IDR analysis, please check your bed files and/or idr parameters."
                            exit 1 
                        elif [ \$npeaks -gt ${params.idr_maxpeaks} ]
                        then
                            sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n ${params.idr_maxpeaks}  > \${name}.regionPeak
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
        else {
    
        // CREATE CHANNEL TO GROUP SAMPLE ID WITH ASSOCIATED REPLICATES (BED FILES)
            T1BEDrep
                .map { it -> [ it[0].split('_')[0..-3].join('_'), it[1] ] }
                .groupTuple(by: [0])
                .map { it -> it[0,1].flatten() }
                .set { T1BED_replicate_ch }
                //.println()   
    
    
            // PROCESS 19 : createPseudoReplicates (CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
            // What it does : creates 2 pseudo replicates per true replicates, then pool the true replicates, and creates 2 pseudo replicates from this pool.
            // Input : bed files from type1 aligned SSDNA, chip and control
            // Output : The 2 true replicates ; the pool of the 2 true replicates ; the 4 pseudo replicates from true replicates; 
            // the 2 pseudo replicates from the pool of true replicates.
            // So 9 files in total.
            process createPseudoReplicates {
                tag "${id_ip}"
                label 'process_basic'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/pseudo_replicates",  mode: params.publishdir_mode
                input:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2) from T1BED_replicate_ch
                output:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2), file('*r1_pseudorep_r1.bed'), file('*r1_pseudorep_r2.bed'), \
                        file('*r2_pseudorep_r1.bed'), file('*r2_pseudorep_r2.bed'), file('*pool_r1.bed'), file('*pool_r2.bed'), \
                        file('*poolT.bed') into ALLREP
                    val 'ok' into createPseudoReplicates_ok 
                script:
                """
                # Shuffle then split the 2 original replicates in 2 pseudo replicates
                nlines_r1=\$((`cat ${t1_rep1} | wc -l`/2)) 
                nlines_r2=\$((`cat ${t1_rep2} | wc -l`/2))
                shuf ${t1_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
                shuf ${t1_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
                mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
                mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
                mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
                mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed
    
                # Pool the 2 original replicates then shuffle then split in 2 pseudo replicates
                nlines_pool=\$((`cat ${t1_rep1} ${t1_rep2} | wc -l`/2))
                cat ${t1_rep1} ${t1_rep2} > ${id_ip}_poolT.bed
                shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
                mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
                mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
                """
            }
    
            // PROCESS 20 : callPeaksForIDR (CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
            // What it does : uses macs2 callPeak to perform peak-calling on all replicates (2), pseudo replicates (4), pool (1), pool pseudo replicates (2)
            // Input : all 9 bed files from true replicates and pseudo replicates creation; and genome size from process 14 
            // Output : 9 regionPeak files from peak calling
            process callPeaksForIDR {
                tag "${id_ip}"
                label 'process_medium'
                label 'peakcalling'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2", mode: params.publishdir_mode, pattern: '*narrowPeak*'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2", mode: params.publishdir_mode, pattern: '*regionPeak*'
                publishDir "${params.outdir}/01_logs/idr/${control_status}/${idr_params}/macs2",	mode: params.publishdir_mode, pattern: '*.macs2.log'            
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/macs2", mode: params.publishdir_mode, pattern: '*'
                input:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2), file(r1_pseudorep_r1), file(r1_pseudorep_r2),\
                        file(r2_pseudorep_r1), file(r2_pseudorep_r2), file(pool_r1), file(pool_r2), \
                        file(poolT) from ALLREP
                    val(genome_size) from gsize
                output:
                    tuple val(id_ip), file(t1_rep1), file(t1_rep2), file('*_R1*type1*.regionPeak'), file('*_R2*type1*.regionPeak'),\
                        file('*r1_pseudorep_r1*.regionPeak'), file('*r1_pseudorep_r2*.regionPeak'),\
                        file('*r2_pseudorep_r1*.regionPeak'), file('*r2_pseudorep_r2*.regionPeak'), \
                        file('*pool_r1*.regionPeak'), file('*pool_r2*.regionPeak'), file('*poolT*.regionPeak') into ALLPEAKSREP
                    file('*narrowPeak*')
                    file('*.macs2.log')
                    val 'ok' into callPeaksForIDR_ok
                script:
                """
                # Runs macs2 callpeak on all input bed files
                for file in \$(ls *.bed) ; 
                do
                    # Get basename to name the outputs files according to input file names
                    name=`basename -- \$file .bed`
                    random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
                    mkdir \${random_id}
                    # Call peaks with macs2, using pvalue or qvalue according to paramters
                    if [[ ${params.idr_macs_pv} == -1 && ${params.idr_macs_qv} == -1 ]];
                    then
                        macs2 callpeak -g ${genome_size} -t \$file \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                        --name \$name --nomodel --extsize ${params.macs_extsize} > \$name.macs2.log 2>&1
                    elif [ ${params.idr_macs_pv} != -1 ];
                    then
                        macs2 callpeak -g ${genome_size} -t \$file \
                        --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                        --name \$name --nomodel --extsize ${params.macs_extsize} --pvalue ${params.idr_macs_pv} > \$name.macs2.log 2>&1
                    elif [ ${params.idr_macs_qv} != -1 ];
                    then
                        macs2 callpeak -g ${genome_size} -t \$file \
                            --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                            --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.idr_macs_qv} > \$name.macs2.log 2>&1
                    fi
    
                    # Cut resulting bed files to 10000 lines max for better handling in IDR statistical analysis
                    npeaks=`cat \${name}_peaks.narrowPeak | wc -l`
                    # Check if the number of peaks is not too low
                    if [ \$npeaks -lt 20 ]
                    then
                        echo "Not enough peaks in \${name}_peaks.narrowPeak file, to perform IDR analysis, please check your bed files and/or idr parameters."
                        exit 1 
                    elif [ \$npeaks -gt ${params.idr_maxpeaks} ]
                    then
                        sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n ${params.idr_maxpeaks}  > \${name}.regionPeak
                    else
                        npeaks=\$(( \$npeaks - 1 ))
                        sort -k 8nr,8nr \${name}_peaks.narrowPeak | head -n \$npeaks  > \${name}.regionPeak
                    fi
                done
                """
            }
        }
    }
    
    if (params.with_idr && params.nb_replicates == "2" ) {
    
        blacklist_ch2
            .combine(chrsize_ch2)
            .set { blacklist_chrsize_ch }
    
        blacklist_chrsize_ch
            .combine(encode_idr_ch)
            .set { blacklist_chrsize_encode_idr_ch }
    
        blacklist_chrsize_encode_idr_ch
            .combine(index_fai_ch4)
            .set { blacklist_chrsize_encode_idr_fai_ch } 
    
        ALLPEAKSREP
            .combine(blacklist_chrsize_encode_idr_fai_ch)
            .set { ALLPEAKSREP_pack }
    
        // PROCESS 21 : IDRanalysis (PERFORM IDR ANALYSIS ON 4 PAIRS OF REPLICATES OR PSEUDOREPLICATES)
        // What it does : performs IDR (Irreproducible Discovery Rate) statistical analysis on 4 pairs of replicates :
        // IDR1 and IDR2 : pseudoreplicates from true replicates ; IDR3 : true replicates ; IDR4 : pseudoreplicates from pooled true replicates ;
        // Input : 9 regionpeak files from idr peak calling  
        // Ouptut : log and results files from IDR analysis
        // External tool : IDR python script from encode chip-seq pipeline
        process IDRanalysis {
            tag "${id_ip}"
            label 'process_medium'
            label 'idr'
    	publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/plot",               mode: params.publishdir_mode, pattern: '*.png'
            publishDir "${params.outdir}/01_logs/idr/${control_status}/${idr_params}",            		 mode: params.publishdir_mode, pattern: '*.log'
            publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/unthresholded-peaks",mode: params.publishdir_mode, pattern: '*unthresholded-peaks*'
            publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/bfilt",              mode: params.publishdir_mode, pattern: '*.bfilt.gz'
            publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/peaks",              mode: params.publishdir_mode, pattern: '*Peak.gz'
            input:
                tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(ip_rep1), path(ip_rep2), path(r1_pseudorep_r1), path(r1_pseudorep_r2), \
                    path(r2_pseudorep_r1), path(r2_pseudorep_r2), path(pool_r1), path(pool_r2), path(poolT), path(blacklist), path(chrsize), \
                    path(encode_idr_script), path(fai) from ALLPEAKSREP_pack
            output:
                file('*unthresholded-peaks*')
                file('*.log')
                file('*.png')
                file('*.bfilt.gz')
                file('*Peak.gz')
                tuple val(id_ip), file(t1_rep1), file(t1_rep2), file('*r1_pseudorep*.*Peak.gz'), file('*r2_pseudorep*.*Peak.gz'), \
                    file('*truerep*.*Peak.gz'), file('*poolrep*.*Peak.gz') into IDRPEAKS
                val 'ok' into IDRanalysis_ok
            script:
            """
            #Set up blacklist parameter if required
            if [ ${params.blacklist} != "None" ]
            then
                blacklist_params="--blacklist ${blacklist}"
            else
                blacklist_params=""
            fi
            
            #Get genome size and number of relaxed peaks to set up idr p-values thresholds 
            #(see https://sites.google.com/site/anshulkundaje/projects/idr/deprecated?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1)
            genome_size=`cut -f3 ${fai} |tail -n1`
    
            npeaksR1=`cat ${ip_rep1} | wc -l`
            npeaksR2=`cat ${ip_rep2} | wc -l`
            echo \$genome_size
            echo \$npeaksR1 
            echo \$npeaksR2
            if [ \$genome_size -lt 500000000 ] && [ "${params.idr_setup}" = "auto" ]
            then
                idr_threshold_r1=0.01
                idr_threshold_r2=0.01
                idr_threshold_truerep=0.01
                idr_threshold_poolrep=0.01
            elif [ \$genome_size -gt 500000000 ] && [ "${params.idr_setup}" = "auto" ] 
            then
                if [[ \$npeaksR1 -gt 150000 && \$npeaksR2 -gt 150000 ]]
                then 
                    idr_threshold_r1=0.01
                    idr_threshold_r2=0.01
                    idr_threshold_truerep=0.01
                    idr_threshold_poolrep=0.0025
                else  
                    idr_threshold_r1=0.05
                    idr_threshold_r2=0.05
                    idr_threshold_truerep=0.05
                    idr_threshold_poolrep=0.01   
                fi
            elif [ "${params.idr_setup}" = "custom" ]
            then
                idr_threshold_r1=${params.idr_threshold_r1}
                idr_threshold_r2=${params.idr_threshold_r2}
                idr_threshold_truerep=${params.idr_threshold_truerep}
                idr_threshold_poolrep=${params.idr_threshold_poolrep}
            fi
            
            #IDR1
            if [ -s ${r1_pseudorep_r1} ] && [ -s ${r1_pseudorep_r2} ] && [ -s ${ip_rep1} ]
            then
            python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
                --idr-thresh \$idr_threshold_r1 --idr-rank ${params.idr_rank} \
                --chrsz ${chrsize} \$blacklist_params \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix ${id_ip}_r1_pseudorep \
                ${r1_pseudorep_r1} ${r1_pseudorep_r2} ${ip_rep1} > ${id_ip}_r1_pseudorep.encodeidr.out.log 2>&1
            fi
            #IDR2
            if [ -s ${r1_pseudorep_r1} ] && [ -s ${r1_pseudorep_r2} ] && [ -s ${ip_rep2} ]
            then
            python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
                --idr-thresh \$idr_threshold_r2 --idr-rank ${params.idr_rank} \
                --chrsz ${chrsize} \$blacklist_params \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix ${id_ip}_r2_pseudorep \
                ${r2_pseudorep_r1} ${r2_pseudorep_r2} ${ip_rep2} > ${id_ip}_r2_pseudorep.encodeidr.out.log 2>&1
            fi
            #IDR3
            if [ -s ${ip_rep1} ] && [ -s ${ip_rep2} ] && [ -s ${poolT} ]
            then
            python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
                --idr-thresh \$idr_threshold_truerep --idr-rank ${params.idr_rank} \
                --chrsz ${chrsize} \$blacklist_params \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix ${id_ip}_truerep \
                ${ip_rep1} ${ip_rep2} ${poolT} > ${id_ip}_truerep.encodeidr.out.log 2>&1
            fi
            #IDR4
            if [ -s ${pool_r1} ] && [ -s ${pool_r2} ] && [ -s ${poolT} ]
            then
            python ${encode_idr_script} --peak-type ${params.idr_peaktype} \
                --idr-thresh \$idr_threshold_poolrep --idr-rank ${params.idr_rank} \
                --chrsz ${chrsize} \$blacklist_params \
                --regex-bfilt-peak-chr-name ${params.idr_filtering_pattern} \
                --prefix ${id_ip}_poolrep \
                ${pool_r1} ${pool_r2} ${poolT} > ${id_ip}_poolrep.encodeidr.out.log 2>&1
            fi
    
            #rename blackist filtered file for better handling in the output
            find . -name "${id_ip}_r1_pseudorep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_r1_pseudorep.bfilt.gz \\;
            find . -name "${id_ip}_r2_pseudorep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_r2_pseudorep.bfilt.gz \\;
            find . -name "${id_ip}_truerep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {}  ${id_ip}_truerep.bfilt.gz \\;
            find . -name "${id_ip}_poolrep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_poolrep.bfilt.gz \\;
    
            """
        }
    
        // PROCESS 22 : IDRpostprocess (IDR PEAKS POST PROCESSING)
        // What it does : Compute rescue ratio and self consistency ratio to evaluate reproducibility
        // and export reproducible set of peaks
        // Input : peaks called in IDR process
        // Output : final set of peaks post IDR and text file with reproducibility ratio and flags
        process IDRpostprocess {
            tag "${id_ip}"
            label 'process_basic'
            label 'idr'
                publishDir "${params.outdir}/01_logs/idr/${control_status}/${idr_params}/qc",   		    mode: params.publishdir_mode, pattern:'*.log'
                publishDir "${params.outdir}/02_results/idr/${control_status}/${idr_params}/peaks",             mode: params.publishdir_mode, pattern:'*.finalPeaks*'
                publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks/idr/${idr_params}",  mode: params.publishdir_mode, pattern:'*.finalPeaks*'
            input:
                tuple val(id_ip), file(t1_rep1), file(t1_rep2), file(r1_pseudorep_peaks), file(r2_pseudorep_peaks), \
                    file(truerep_peaks), file(poolrep_peaks) from IDRPEAKS
            output:
                tuple val(id_ip), file(t1_rep1), file(t1_rep2), file('*.finalPeaks_IDR.bed') into ALLPEAKSTOPPIDR
                file('*_IDR_QC.log') 
                file('*.finalPeaks*')
                val 'ok' into IDRpostprocess_ok
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
                zcat ${poolrep_peaks} > ${poolrep_peaks}.finalPeaks_IDR.bed
                #filter for chromY, chromM and blacklist ? #todo
            else
                rr=`python -c "print(\$nt/\$np)"`
                zcat ${truerep_peaks} > ${truerep_peaks}.finalPeaks_IDR.bed
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
    }
    
    if (!params.with_idr) {
        createPseudoReplicates_ok = Channel.value( 'ok' )
        callPeaksForIDR_ok = Channel.value( 'ok' )
        IDRanalysis_ok = Channel.value( 'ok' )
        IDRpostprocess_ok = Channel.value( 'ok' )
    }
    
    //***************************************************************************//
    //                     SECTION 8 : PEAK POST PROCESSING                      //
    //***************************************************************************//        
    
    // PROCESS 23 : normalizePeaks (CENTER AND NORMALIZE PEAKS)
    // What it does : perform the peak centering using Kevin Brick's method :
    // 1. Recenters the peaks by the median of the F/R dists
    // 2. Calculate the in-peak background in a number of ways.
    // 3. Output recentered peaks with strength. Strength is normalized by the number
    //    of wrong-direction fragments (i.e. REV to left, FWD to right) and expressed
    //    as RPKM, using the total normalized in-hotspot tag count as a denominator.
    // Input : final peak set bed file and mapped T1 bed file
    // Output : bedgraph and tab
    // External tool : Perl script from K. Brick (original pipeline, 2012)
    if (!params.with_idr) {
    
        norm_ch
            .combine(reverse_ch)
            .set { norm_reverse_ch }
    
        ALLPEAKSTOPP
            .combine(norm_reverse_ch)
            .set { ALLPEAKSTOPP_ch }
    
        process normalizePeaks {
            tag "${id_ip}"
            label 'process_basic'
            label 'bigwig'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/normalized/no-idr",        mode: params.publishdir_mode, pattern: '*.bedgraph'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/normalized/no-idr/tab",    mode: params.publishdir_mode, pattern: '*.tab'
            publishDir "${params.outdir}/01_logs/peaks/${control_status}/normalized/no-idr",mode: params.publishdir_mode, pattern: '*.log'
            input:
                tuple val(id_ip), path(ip_bed), path(peaks_bed), path(norm_script), path(reverse_script) from ALLPEAKSTOPP_ch
            output:
               tuple val(id_ip), file("*normpeaks.bedgraph"), file("*normpeaks.tab")
                file('*.normpeaks.log')
                val('ok') into normalizePeaks_ok
            script:
            """
            # Normalize and recenter peaks
            perl ${norm_script} --bed ${peaks_bed} \
                --in ${ip_bed} --out ${id_ip}.no-idr.normpeaks.bedgraph \
                --rc --rev_src ${reverse_script} > ${id_ip}.no-idr.normpeaks.log 2>&1
            """
        }
    
        //IF RUNNING WITH REPLICATES, THEN MERGE FINAL PEAKS FROM  REPLICATES
        if (params.nb_replicates == "2") {
    
        //CREATE CHANNEL TO GROUP SAMPLES BY REPLICATES
        //The final channel is composed of (1+nb_replicates) elements : sampleID, nb_replicates*(peak bed file)
            ALLPEAKSTOPP2
                .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1], it[2] ] }
                .groupTuple(by: [0])
                .map { it -> it[0,2].flatten() }
                .set { ALLPEAKSTOPP2 }
                //.println()
    
            // PROCESS 24 : mergeFinalPeaks (MERGE PEAKS FROM REPLICATES)
            // What it does : Merge peaks bed files from replicates
            // Input :  peak bed files, grouped by replicates, from call_peaks process
            // Output : merged bed file for all replicates
            process mergeFinalPeaks {
                tag "${id_ip}"
                label 'peakcalling'
                label 'process_low'
                publishDir "${params.outdir}/02_results/peaks/${control_status}/finalpeaks/no-idr",      mode: params.publishdir_mode, pattern: '*.bed'
                publishDir "${params.outdir}/01_logs/peaks/${control_status}/finalpeaks/no-idr",  mode: params.publishdir_mode, pattern: '*.log'
                input:
                    tuple val(id_ip), file(R1peaks), file(R2peaks) from ALLPEAKSTOPP2
                output:
                    file('*mergePeaks.bed')
                    file('*.log')
                    val('ok') into mergeFinalPeaks_ok
                script:
                """
                # Concatenate, sort, then merge peak sets from replicates
                cat ${R1peaks} ${R2peaks} | sort -k1,1 -k2,2n > ${id_ip}_R1R2_no-idr.mergePeaks.bed
                mergeBed -i ${id_ip}_R1R2_no-idr.mergePeaks.bed >& ${id_ip}_R1R2_no-idr.mergeBed.log 2>&1
                """
            }
        }
        else {
            mergeFinalPeaks_ok = Channel.value( 'ok' )
        }
    }
    
    else if ( params.with_idr && params.nb_replicates == "2" ) {
    
        norm_ch
            .combine(reverse_ch)
            .set { norm_reverse_ch }
    
        ALLPEAKSTOPPIDR
            .combine(norm_reverse_ch)
            .set { ALLPEAKSTOPPIDR_ch }
    
    // PROCESS 23 : normalizePeaks (CENTER AND NORMALIZE PEAKS)
    // What it does : perform the peak centering using Kevin Brick's method :
    // 1. Recenters the peaks by the median of the F/R dists
    // 2. Calculate the in-peak background in a number of ways.
    // 3. Output recentered peaks with strength. Strength is normalized by the number
    //    of wrong-direction fragments (i.e. REV to left, FWD to right) and expressed
    //    as RPKM, using the total normalized in-hotspot tag count as a denominator.
    // Input : final peak set bed file and mapped T1 bed file
    // Output : bedgraph and tab
    // External tool : Perl script from K. Brick (original pipeline, 2012)
       process normalizePeaks_idr {
            tag "${id_ip}"
            label 'process_basic'
            label 'bigwig'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/normalized/idr/${idr_params}",             mode: params.publishdir_mode, pattern: '*.bedgraph'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/normalized/idr/${idr_params}/tab",         mode: params.publishdir_mode, pattern: '*.tab'
            publishDir "${params.outdir}/01_logs/peaks/${control_status}/normalized/idr/${idr_params}",     mode: params.publishdir_mode, pattern: '*.log'
            input:
                tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(peaks_bed), path(norm_script), path(reverse_script) from ALLPEAKSTOPPIDR_ch
            output:
                tuple val(id_ip), file("*normpeaks*.bedgraph"), file("*normpeaks*.tab")
                file('*.normpeaks*.log')
                val('ok') into normalizePeaks_idr_ok
                val('ok') into mergeFinalPeaks_ok
    	    val('ok') into normalizePeaks_ok
            script:
            """
            # Merge the 2 T1BED replicates and sort
            cat ${t1_rep1} ${t1_rep2} | sort -k1,1 -k2,2n > ${id_ip}_rep1_rep2.bed
            
            # Normalize and recenter peaks
            perl ${norm_script} --bed ${peaks_bed} \
                --in ${id_ip}_rep1_rep2.bed --out ${id_ip}.idr.normpeaks.bedgraph \
                --rc --rev_src ${reverse_script} > ${id_ip}.idr.normpeaks.log 2>&1
    	perl ${norm_script} --bed ${peaks_bed} \
    	    --in ${t1_rep1} --out ${id_ip}.idr.normpeaks.R1.bedgraph \
    	    --rc --rev_src ${reverse_script} > ${id_ip}.idr.normpeaks.R1.log 2>&1
    	perl ${norm_script} --bed ${peaks_bed} \
    	     --in ${t1_rep2} --out ${id_ip}.idr.normpeaks.R2.bedgraph \
    	     --rc --rev_src ${reverse_script} > ${id_ip}.idr.normpeaks.R2.log 2>&1
            """
        }
    }
    
    else {
        normalizePeaks_ok = Channel.value( 'ok' )
        mergeFinalPeaks_ok = Channel.value( 'ok' )
    }
    
    //***************************************************************************//
    //                       SECTION 8 : SATURATION  CURVE                       //
    //***************************************************************************//
    // Create channel containing peaks called in progressively downsampled samples from previous process
    //N_peak_ch = satcurvebed_ch.collect()
    
    getPeaksBedFiles_ch
        .combine(runSatCurve_ch)
        .combine(satCurveHS_ch)
        .set { satcurve_scripts_ch }
    
    // PROCESS 25 : makeSatCurve (CREATE SATURATION CURVE)
    // What it does : Compute saturation curve of the samples 
    // The saturation curve plots the number of peaks called function of the number of reads in the samples
    // Input : All bed files from callPeaks process corresponding to peaks called in progressively downsampled samples
    // Ouptut : png of the saturation curve and data table
    // External tool : Perl script from K. Brick (original pipeline 2012) and R script (adapted from original pipeline 2012)
    if (params.satcurve) {
        process makeSatCurve {
            tag "${global_name}"
            label 'process_basic'
            label 'peakcalling'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/saturation_curve/${params.sctype}",  mode: params.publishdir_mode, pattern: '*.tab'
            publishDir "${params.outdir}/02_results/peaks/${control_status}/saturation_curve/${params.sctype}",  mode: params.publishdir_mode, pattern: '*.png'
            publishDir "${params.outdir}/01_logs/peaks/${control_status}/saturation_curve/${params.sctype}",     mode: params.publishdir_mode, pattern: '*.log'
            input:
                path(saturation_curve_data) from allbed.collect()
                tuple path(getPeaksBedFiles_script), path(runSatCurve_script), path(satCurveHS_script) from satcurve_scripts_ch
            output:
                path("*.tab") into satcurve_tab_ch
                path("*.png") into satcurve_ch
                val 'ok' into makeSatCurve_ok
                path("*.log")
            script:
            """
            # Get number of peaks in samples in the ${params.outdir}/saturation_curve/peaks directory (from callPeaks process)
            perl ${getPeaksBedFiles_script} -tf satCurve.tab -dir "." >& ${global_name}_makeSatCurve.log 2>&1
            # Plot saturation curve
            Rscript ${runSatCurve_script} ${satCurveHS_script} satCurve.tab ${global_name} >& ${global_name}_makeSatCurve_Rscript.log 2>&1
            #mv ${global_name}.saturationCurve.png ${global_name}.saturationCurve_mqc.png
            """
        }
    }
    else {
        makeSatCurve_ok = Channel.value( 'ok' )
        satcurve_ch = Channel.fromPath(emptyfile6)
    }
    
    //***************************************************************************//
    //                          SECTION 9 : GENERAL QC                           //
    //***************************************************************************//
    // PROCESS 26 : general_multiqc (GENERATES GENERAL MULTIQC REPORT)
    // What it does : Compute multiCQ report for the analysis based on output folder content
    // Input : all termination tags of processes so that this process only runs at the end of the pipeline
    // Output : multiQC HTML report
    process general_multiqc {
        tag "${global_name}"
        label 'process_basic'
        label 'multiqc'
        label 'internet_access'
        publishDir "${params.outdir}/02_results/qc/multiqc",  mode: 'copy'
        publishDir "${params.outdir}/01_logs/multiqc",  	  mode: 'copy', pattern: '*.log'
        input:
            path(logo) from logo_ch
            path(multiqc_configfile) from multiqc_configfile_ch
            path(png_fingerprint) from png_fingerprint_ch
            path(tab_fingerprint) from tab_fingerprint_ch
            path(bwa_flagstat) from bwa_flagstat_ch
            path(itr_flastat) from itr_flastat_ch
            path(ssdsreport) from ssdsreport_ch
            path(ssds_report_denovo) from ssds_report_denovo_ch
            path(multiqc_denovo) from multiqc_denovo_ch
            path(multiqc) from multiqc_ch
            path(samstats) from samstats_ch
            path(raw_fastqc) from raw_fastqc_ch
            path(trim_fastqc) from trim_fastqc_ch
            path(satcurve) from satcurve_ch
            val('trimming_ok') from trimming_ok.collect().ifEmpty([])
       	    val('fastqc_ok') from fastqc_ok.collect().ifEmpty([])
    	    val('bwa_ok') from bwa_ok.collect().ifEmpty([])
    	    val('filterbam_ok') from filterbam_ok.collect().ifEmpty([])
    	    val('parseitr_ok') from parseitr_ok.collect().ifEmpty([])
    	    val('makeBigwig_ok') from makeBigwig_ok.collect().ifEmpty([])
            val('makeBigwigReplicates_ok') from makeBigwigReplicates_ok.collect().ifEmpty([])
            val('kb_bigwig_ok') from kb_bigwig_ok.collect().ifEmpty([])
    	    val('fr_bigwig_ok') from fr_bigwig_ok.collect().ifEmpty([])
    	    val('shufbed_ok') from shufbed_ok.collect().ifEmpty([])
    	    val('callPeaks_ok') from callPeaks_ok.collect().ifEmpty([])
    	    val('samstats_ok') from samstats_ok.collect().ifEmpty([])
            val('ssreport_ok') from ssreport_ok.collect().ifEmpty([])
    	    val('makefingerprint_ok') from makefingerprint_ok.collect().ifEmpty([])
    	    val('computefrip_ok') from computefrip_ok.collect().ifEmpty([])
    	    val('multiqc_denovo_ok') from multiqc_denovo_ok.collect().ifEmpty([])
    	    val('ssds_multiqc_ok') from ssds_multiqc_ok.collect().ifEmpty([])
    	    val('createPseudoReplicates_ok') from createPseudoReplicates_ok.collect().ifEmpty([])
    	    val('callPeaksForIDR_ok') from callPeaksForIDR_ok.collect().ifEmpty([])
    	    val('IDRanalysis_ok') from IDRanalysis_ok.collect().ifEmpty([])
    	    val('IDRpostprocess_ok') from IDRpostprocess_ok.collect().ifEmpty([])
    	    val('normalizePeaks_ok') from normalizePeaks_ok.collect().ifEmpty([])
    	    val('mergeFinalPeaks_ok') from mergeFinalPeaks_ok.collect().ifEmpty([])
    	    val('makeSatCurve_ok') from makeSatCurve_ok.collect().ifEmpty([])
        output:
    	    file('*')
        script:
        """
        multiqc --export -c ${multiqc_configfile} \
    		-n ${global_name}.multiQC.quality-control.report "." \
            	>& ${global_name}.multiQC.quality-control.log 2>&1
        
        """
    }
}

//***************************************************************************//
//                                                                           //
//                          END OF PIPELINE !!                               //
//                                                                           //
//***************************************************************************//

// PRINT LOG MESSAGE ON COMPLETION        
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Command line: $workflow.commandLine"
    println "Execution profile: $workflow.profile"
    println "Script ID: $workflow.scriptId"
    println "Run name: $workflow.runName"

}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Error report: ${workflow.errorReport}"
}




