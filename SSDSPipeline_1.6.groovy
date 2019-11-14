#!/usr/bin/env nextflow

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "SSDS PIPELINE (Version 1.6_NF)                                "
  log.info "Author: Kevin Brick                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "

  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/pipelines/pipeIt \\"
  log.info " "
  log.info "/data/RDCO/code/pipelines/pipeIt --h for help \\ "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run $NXF_PIPEDIR/SSDSPipeline_1.6.groovy \\"
  log.info " --bam           <bam file>"
  log.info " --sra           <sra name (SRRXXXXX)>"
  log.info " --obj           <object store regex>"
  log.info " --fq1           <fastq read1 file>"
  log.info " --fq2           <fastq read2 file>"
  log.info " --r1Len         <read1 length>"
  log.info " --r2Len         <read2 length>"
  log.info " --genome        <string> \\"
  log.info " --threads       <number of threads: default = 6> \\"
  log.info " --outName       <string default=bam file stem> \\"
  log.info " --outdir        <string: default = outName> \\"
  log.info " --sample_name   <string: default = outName> \\"
  log.info " --platform      <string: default = ILLUMINA > \\"
  log.info " --platform_unit <string: default=HISEQ2500> \\"
  log.info " --library       <string> \\"
  log.info " --rundate       <YYMMDD: default = 000101> \\"
  log.info " --outdir_tmp    <string>/tmp \\"
  log.info " -with-trace -with-timeline"
  log.info " "
  log.info "HELP: nextflow run $NXF_PIPEDIR/SSDSPipeline_1.6.groovy --help"
  log.info " "
  log.info "=========================================================================="
  log.info "Required Arguments:"
  log.info " "
  log.info "          --bam				STRING     BAM FILE                     OR"
  log.info "          --sra				STRING     SRR SAMPLE NAME (SRR only)   OR"
  log.info "          --obj				STRING     Regex for object store       OR"
  log.info "          --fq1				STRING     READ1 FASTQ FILE"
  log.info "          --fq2				STRING     READ2 FASTQ FILE"
  log.info "          --genome          STRING     FASTA FILE (GENOME REFERENCE)"
  log.info " "
  log.info "Output and Tempory directory Arguments"
  log.info " "
  log.info "          --outName      STRING        Output BAM file name"
  log.info "          --outdir          DIR        Output directory"
  log.info "          --outdir_tmp      DIR        Tempory directory"
  log.info " "
  log.info "SAM Read Group Arguments:"
  log.info " "
  log.info "          --rg_id           STRING     SAM Read Group Tag RG:ID (defaults to outName)"
  log.info "          --sample_name     STRING     SAM Read Group Tag RG:SM"
  log.info "          --platform        STRING     SAM Read Group Tag RG:PL"
  log.info "          --platform_unit   STRING     SAM Read Group Tag RG:PU"
  log.info "          --library         STRING     SAM Read Group Tag RG:LB"
  log.info "          --rundate         STRING     SAM Read Group Tag RG:DT"
  log.info ""
  log.info "====================================================================================="
  log.info "Optional Pipeline Arguments:"
  log.info "====================================================================================="
  log.info " "
  log.info "Optional BWA Alignment Arguments:"
  log.info " "
  log.info "          --bwa_args        STRING     Optional bwa arguments eg: \"-I 250,50\""
  log.info " "
  log.info "====================================================================================="
  exit 1

}

// Define arg defaults

//number of threads
params.threads = 16

//genome version
params.genome_fasta = "$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"

//sample information
params.bam = ""
params.sra = ""
params.obj = ""
params.fq1 = ""
params.fq2 = ""
params.r1Len = "36"
params.r2Len = "40"
params.name = ""
params.outName = "${params.name}.SSDS.${params.genome}"
params.tmpName = "${params.name}.tmpFile.${params.genome}"
params.rg_id = "${params.name}"
params.sample_name=""
params.library=""
params.rundate=""
params.platform="ILLUMINA"
params.platform_unit="HISEQ2500"
params.mem="16G"
params.bwaSplitSz = 20000000

//output and tempory directories
params.outdir = "./${params.outName}"
params.outdir_tmp = "/tmp"

params.bamPGline = '@PG	ID:ssDNAPipeline1.6_KBRICK	VN:1.6nxf'

//input BAM file
if (params.bam){
	params.bai = "${params.bam}.bai"
}

if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq1){inputType = 'fastqSR'}
if (params.fq2){inputType = 'fastqPE'}

//log.info
log.info "===================================================================="
log.info "SSDS (bam) PIPELINE : Align ssDNA"
log.info "===================================================================="
log.info "ref genome         : ${params.genome}"
log.info "genome fasta       : ${params.genome_fasta}"
log.info "bam                : ${params.bam}"
log.info "sra                : ${params.sra}"
log.info "obj                : ${params.obj}"
log.info "fq1                : ${params.fq1}"
log.info "fq2                : ${params.fq2}"
log.info "r1Len              : ${params.r1Len}"
log.info "r2Len              : ${params.r2Len}"
log.info "RG:SM              : ${params.sample_name}"
log.info "RG:PL              : ${params.platform}"
log.info "RG:PU              : ${params.platform_unit}"
log.info "RG:LB              : ${params.library}"
log.info "RG:ID              : ${params.rg_id}"
log.info "RG:DT              : ${params.rundate}"
log.info "outName            : ${params.outName}"
log.info "outdir             : ${params.outdir}"
log.info "temp_dir           : ${params.outdir_tmp}"
log.info "threads            : ${params.threads}"
log.info "mem                : ${params.mem}"
log.info "max reads for bwa  : ${params.bwaSplitSz}"

process getFQs{
	scratch '/lscratch/$SLURM_JOBID'
		clusterOptions ' --gres=lscratch:500 --partition=norm'
		echo true
		cpus 12
		memory '32g'

		time { 24.hour * task.attempt }
    errorStrategy { 'retry' }
    maxRetries 1

		module 'nextflow/0.30.2'

		//publishDir params.outdir, mode: 'copy', overwrite: false
		publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*zip"
		publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*html"
		publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*png"
		publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*txt"

		output:
			//file '*.R1.fastq' optional true into initFQ1
			//file '*.R2.fastq' optional true into initFQ2
			set val(params.outName),file('*.R1.fastq'),file('*.R2.fastq') into FQx1,FQx1b
			file '*.R2.fastq'   optional true into initFQ2
			file '*fastqc*'     optional true into fastQCOut
			file '*screen*'     optional true into fastQScreenOut

		script:
			switch (inputType) {
        		case 'sra':
					"""
					nextflow run -c \$NXF_PIPEDIR/nextflow.local.config \$NXF_PIPEDIR/getFQ.groovy \
					--genome ${params.genome} --sra ${params.sra} --outdir . --withFQC false
					"""
	            	break

				case 'obj':
					"""
					nextflow run -c \$NXF_PIPEDIR/nextflow.local.config \$NXF_PIPEDIR/getFQ.groovy \
					--genome ${params.genome} --obj ${params.obj} --outdir . --withFQC false
					"""
					break

				case 'bam':
					"""
					nextflow run -c \$NXF_PIPEDIR/nextflow.local.config \$NXF_PIPEDIR/getFQ.groovy \
					--genome ${params.genome} --bam ${params.bam} --outdir . --withFQC false
					"""
					break

				case 'fastqSR':
					"""
					nextflow run -c \$NXF_PIPEDIR/nextflow.local.config \$NXF_PIPEDIR/getFQ.groovy \
					--genome ${params.genome} --fq1 ${params.fq1} --outdir . --withFQC false
					"""
					break

				case 'fastqPE':
					"""
					nextflow run -c \$NXF_PIPEDIR/nextflow.local.config \$NXF_PIPEDIR/getFQ.groovy \
					--genome ${params.genome} --fq1 ${params.fq1} --fq2 ${params.fq2} --outdir . --withFQC false
					"""
					break
			}
  }

process runFASTQC {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 16
  memory '16g'
	time '8h'

	module 'fastqtools/0.8'
	module 'fastqc/0.11.8'
	module 'bwa/0.7.17'

	tag { params.outName }

	publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*zip"
	publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*html"
	publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*png"
	publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*txt"

	input:
	  set val(myNM),file(fqr1),file(fqr2) from FQx1b

	output:
		file '*zip'  into fqcZip
		file '*html' into repHTML
		file '*png'  into repPNG
		file '*txt'  into repTXT

	script:
		"""
		head -n 10000000 ${fqr1} >${params.outName}.R1.fastq
		head -n 10000000 ${fqr2} >${params.outName}.R2.fastq

		fastqc -t ${params.threads} ${params.outName}.R1.fastq
		fastqc -t ${params.threads} ${params.outName}.R2.fastq

		perl \$NXF_PIPEDIR/accessoryFiles/SSDS/scripts/generateFastQCScreenConfig.pl hg19,hg38,mm10,rn6,saccer3,phix,illuminaadapters,univec,bsub,ecoli,canfam3,mondom 5 >fastq_screen.conf

		ln -s ${params.outName}.R1.fastq ${params.outName}.fastq
	   	\$NXF_PIPEDIR/accessoryFiles/SSDS/fastq_screen_v0.11.4/fastq_screen --threads ${params.threads} --force \
	    	                                               --aligner bwa ${params.outName}.fastq \
	    	                                               --conf fastq_screen.conf

		"""
  }

process trimFASTQs {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:300'

  cpus 16
  memory '4g'
  time '8h'

	module 'trimgalore/0.4.5'

	tag { params.outName }
	publishDir params.outdir,  mode: 'copy', pattern: "*_report.*"

	input:
	  set val(myNM),file(fqr1),file(fqr2) from FQx1

	output:
		set val(params.outName),file('*R1.fastq'),file('*R2.fastq') into FQx2, FQx3
		file '*trimming_report.txt' into trimGaloreReport

	script:
		"""
		trimFQ1=`echo "$fqr1" |perl -pi -e 's/.fastq/_val_1.fq/'`
		trimFQ2=`echo "$fqr2" |perl -pi -e 's/.fastq/_val_2.fq/'`

    ## KB 05-30-19: change minimum length from 15->25 ... prevent BWA error
		trim_galore -q 10 --paired --dont_gzip --stringency 6 --length 25 $fqr1 $fqr2

		mv \$trimFQ1 ${params.outName}.R1.fastq
		mv \$trimFQ2 ${params.outName}.R2.fastq

		mv ${fqr1}_trimming_report.txt ${params.outName}.R1_trimgalore_trimming_report.txt
		mv ${fqr2}_trimming_report.txt ${params.outName}.R2_trimgalore_trimming_report.txt

		"""
  }

FQx2.splitFastq( by: params.bwaSplitSz, pe:true, file:true).set { splitFQ }

//Start pipeline
//1st process will convert BAM -> FASTQ then align
process bwaAlign {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 16
	memory 32
	time '60h'

	module 'java/1.8.0_92'
	module 'picard/2.9.2'
	module 'fastxtoolkit/0.0.14'
	// KB 2018-06-04: Biowulf Centos7 update (June 04 2018)
	//module 'samtools/1.5'
	module 'samtools/1.8'
	//END 180604

	tag { params.outName }

	input:
	set val(myNM),file(fqR1),file(fqR2) from splitFQ

	output:
	  file 'bar*bam' into multiBAMaln, multiBAM2merge

	script:
	  // generate random filenames
  	  def nm = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  	  def bamOut = 'bar_' + nm + '.bam'

		"""
		fastx_trimmer -f 1 -l ${params.r1Len} -i $fqR1 -o ${params.tmpName}.R1.fastq
		fastx_trimmer -f 1 -l ${params.r2Len} -i $fqR2 -o ${params.tmpName}.R2.fastq

		\$SSPIPELINEPATH/bwa_0.7.12 aln \
		-t ${task.cpus} \
		${params.genome_fasta} \
		${params.tmpName}.R1.fastq >${params.tmpName}.R1.sai

		\$SSPIPELINEPATH/bwa_ra_0.7.12 aln \
		-t ${task.cpus} \
		${params.genome_fasta} \
		${params.tmpName}.R2.fastq >${params.tmpName}.R2.sai

		\$SSPIPELINEPATH/bwa_0.7.12 sampe \
		${params.genome_fasta} \
		${params.tmpName}.R1.sai \
		${params.tmpName}.R2.sai \
		${params.tmpName}.R1.fastq \
        ${params.tmpName}.R2.fastq >${params.tmpName}.unsorted.sam

		java -jar \$PICARDJAR SamFormatConverter \
                   I=${params.tmpName}.unsorted.sam \
                   O=${params.tmpName}.unsorted.tmpbam \
                   VALIDATION_STRINGENCY=LENIENT 2>pica.err

		java -jar \$PICARDJAR SortSam \
                   I=${params.tmpName}.unsorted.tmpbam \
                   O=${bamOut} \
                   SO=coordinate \
                   VALIDATION_STRINGENCY=LENIENT 2>pica.err

		samtools index ${bamOut}
		ls -l >list.tab
		"""
		//samtools view -Shb /dev/stdin >${params.tmpName}.unsorted.tmpbam
		//-R '@RG\\tID:${params.rg_id}\\tSM:${params.sample_name}\\tPU:${params.platform_unit}\\tPL:${params.platform}\\tLB:${params.library}\\tDT:${params.rundate}'
  }

//Start pipeline
//merge BAMs and merk duplicates
process mergeInitBAMs {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 8
  memory '32g'

	time '6h'

	module 'java/1.8.0_92'
	module 'picard/2.9.2'
	module 'samtools/1.8'

	tag { params.outName }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  file bam_files from multiBAM2merge.collect()

	output:
	  file '*SSDS*.unparsed.bam' into bamAligned
	  file '*SSDS*.unparsed.bam.bai' into bamIDX
	  file '*SSDS*.unparsed.suppAlignments.bam' into bamAlignedSupp
	  file '*SSDS*.unparsed.suppAlignments.bam.bai' into bamIDXSupp
	  file '*MDmetrics.txt' into listz

	script:
	  // get INPUT files as string
	  def input_args = bam_files.collect{ "I=$it" }.join(" ")

		"""
		java -jar \$PICARDJAR MergeSamFiles \
                   $input_args \
                   O=allREADS.bam \
                   AS=true \
                   SO=coordinate \
                   VALIDATION_STRINGENCY=LENIENT 2>>pica.err


    samtools view -F 2048 -hb allREADS.bam >allREADS.ok.bam

		java -jar \$PICARDJAR MarkDuplicatesWithMateCigar \
                   I=allREADS.ok.bam \
                   O=${params.outName}.unparsed.bam \
                   PG=Picard2.9.2_MarkDuplicates \
                   M=${params.outName}.MDmetrics.txt \
                   MINIMUM_DISTANCE=400 \
				   CREATE_INDEX=false \
				   ASSUME_SORT_ORDER=coordinate \
                   VALIDATION_STRINGENCY=LENIENT 2>>pica.err

		samtools index ${params.outName}.unparsed.bam

    samtools view -f 2048 -hb allREADS.bam >allREADS.supp.bam

		java -jar \$PICARDJAR MarkDuplicatesWithMateCigar \
                   I=allREADS.supp.bam \
                   O=${params.outName}.unparsed.suppAlignments.bam \
                   PG=Picard2.9.2_MarkDuplicates \
                   M=${params.outName}.suppAlignments.MDmetrics.txt \
                   MINIMUM_DISTANCE=400 \
				   CREATE_INDEX=false \
				   ASSUME_SORT_ORDER=coordinate \
                   VALIDATION_STRINGENCY=LENIENT 2>>pica.err

		samtools index ${params.outName}.unparsed.suppAlignments.bam
		"""
  }

//Parse ITRs and output SSDS BAMs and BEDs
process parseITRs {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 16
	memory '32g'
	module 'picard/2.9.2'
	// KB 2018-06-04: Biowulf Centos7 update (June 04 2018)
	//module 'samtools/1.5'
	module 'samtools/1.8'
	//END 180604

	time '6h'

  time { 8.hour * task.attempt }
  errorStrategy 'retry'
  maxRetries 2

	tag { params.outName }
	//publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  //file bamAligned
	  //each file(chr_bam) from chrBAM
	  each file(split_bam) from multiBAMaln

	output:
	  file '*e1.bam'         into bamT1
	  file '*e2.bam'         into bamT2
	  file '*dsDNA.bam'      into bamD
	  file '*strict.bam'     into bamDS
	  file '*ied.bam'        into bamUN

	  file '*e1.bam.bai'     into baiT1
	  file '*e2.bam.bai'     into baiT2
	  file '*dsDNA.bam.bai'  into baiD
	  file '*strict.bam.bai' into baiDS
	  file '*ied.bam.bai'    into baiUN

	  file '*e1.bed'         into bedT1
	  file '*e2.bed'         into bedT2
	  file '*dsDNA.bed'      into bedD
	  file '*strict.bed'     into bedDS
	  file '*ied.bed'        into bedUN
	  val 'OK'               into ssdsDone

	script:
    //perl \$SSPIPELINEPATH/ITR_id_v2c_NextFlow2.pl ${split_bam} ${params.genome} >PIPE.out 2>PIPE.err
		"""
		perl \$NXF_PIPEDIR/accessoryFiles/SSDS/scripts/ITR_id_v2c_NextFlow2.pl ${split_bam} ${params.genome}

		sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${split_bam}.ssDNA_type1.bed  -o ${split_bam}.ssDNA_type1.bed
		sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${split_bam}.ssDNA_type2.bed  -o ${split_bam}.ssDNA_type2.bed
		sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${split_bam}.dsDNA.bed        -o ${split_bam}.dsDNA.bed
		sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${split_bam}.dsDNA_strict.bed -o ${split_bam}.dsDNA_strict.bed
		sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${split_bam}.unclassified.bed -o ${split_bam}.unclassified.bed

		samtools view -H ${split_bam}  >header.txt
		echo ${params.bamPGline}      >>header.txt

		cat header.txt ${split_bam}.ssDNA_type1.sam  >${split_bam}.ssDNA_type1.RH.sam
		cat header.txt ${split_bam}.ssDNA_type2.sam  >${split_bam}.ssDNA_type2.RH.sam
		cat header.txt ${split_bam}.dsDNA.sam        >${split_bam}.dsDNA.RH.sam
		cat header.txt ${split_bam}.dsDNA_strict.sam >${split_bam}.dsDNA_strict.RH.sam
		cat header.txt ${split_bam}.unclassified.sam >${split_bam}.unclassified.RH.sam

		samtools view -Shb ${split_bam}.ssDNA_type1.RH.sam  >${split_bam}.ssDNA_type1.US.bam
		samtools view -Shb ${split_bam}.ssDNA_type2.RH.sam  >${split_bam}.ssDNA_type2.US.bam
		samtools view -Shb ${split_bam}.dsDNA.RH.sam        >${split_bam}.dsDNA.US.bam
		samtools view -Shb ${split_bam}.dsDNA_strict.RH.sam >${split_bam}.dsDNA_strict.US.bam
		samtools view -Shb ${split_bam}.unclassified.RH.sam >${split_bam}.unclassified.US.bam

		java -jar \$PICARDJAR SortSam I=${split_bam}.ssDNA_type1.US.bam  O=${split_bam}.ssDNA_type1.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
		java -jar \$PICARDJAR SortSam I=${split_bam}.ssDNA_type2.US.bam  O=${split_bam}.ssDNA_type2.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
		java -jar \$PICARDJAR SortSam I=${split_bam}.dsDNA.US.bam        O=${split_bam}.dsDNA.bam        SO=coordinate VALIDATION_STRINGENCY=LENIENT
		java -jar \$PICARDJAR SortSam I=${split_bam}.dsDNA_strict.US.bam O=${split_bam}.dsDNA_strict.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
		java -jar \$PICARDJAR SortSam I=${split_bam}.unclassified.US.bam O=${split_bam}.unclassified.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT

		samtools index ${split_bam}.ssDNA_type1.bam
		samtools index ${split_bam}.ssDNA_type2.bam
		samtools index ${split_bam}.dsDNA.bam
		samtools index ${split_bam}.dsDNA_strict.bam
		samtools index ${split_bam}.unclassified.bam
		"""
  }

Channel.from(["ssDNA_type1","ssDNA_type2","dsDNA","dsDNA_strict","unclassified"])
       .set {parseType}

process gatherOutputs {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'

	cpus 16
	memory '32g'
	module 'picard/2.9.2'
	// KB 2018-06-04: Biowulf Centos7 update (June 04 2018)
	//module 'samtools/1.5'
	module 'samtools/1.8'
	//END 180604

	time '4h'

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
  java -jar \$PICARDJAR MergeSamFiles O=${params.outName}.${pType}.US.BAM `ls *${pType}.bam | sed 's/.*\$/I=& /'`
  java -jar \$PICARDJAR SortSam I=${params.outName}.${pType}.US.BAM   O=${params.outName}.${pType}.S.BAM   SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR MarkDuplicatesWithMateCigar I=${params.outName}.${pType}.S.BAM  O=${params.outName}.${pType}.bam  PG=Picard2.9.2_MarkDuplicates M=${params.outName}.ssDNA_type1.MDmetrics.txt  CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT

  samtools index ${params.outName}.${pType}.bam

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 -m *${pType}.bed  >${params.outName}.${pType}.bed
  rm -f ${params.outName}.${pType}.US.BAM
  rm -f ${params.outName}.${pType}.S.BAM
    """
  }

//4th process will make a deeptools 150bp smoothed RPKM bigwig
process makeDeeptoolsBigWig {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 16
	memory '32g'

	module 'samtools/1.8'
	module 'deeptools/3.0.1'
	module 'ucsc/365'

	time '6h'

	tag { ibam }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  //file(ibam), from mergeBAMwithIDXdt
	  set file(ibam),file(ibamidx) from mergeBAMwithIDXdt

	output:
	  file '*deeptools.*' into dtDeepTools
	  val 'OK' into dtDone
	  val 'OK' into bwDone

	shell:
		bigwig          = ibam.name.replaceAll(/.bam/,'.deeptools.RPKM.bigwig')
		covPlot         = ibam.name.replaceAll(/.bam/,'.deeptools.coveragePlot.png')
		fingerprintPlot = ibam.name.replaceAll(/.bam/,'.deeptools.fingerprints.png')
		fingerprintData = ibam.name.replaceAll(/.bam/,'.deeptools.fingerprints.tab')
		//samtools index $ibam;
		"""
		touch empty.deeptools.flag
		bamCoverage --bam $ibam --normalizeUsing RPKM --numberOfProcessors ${params.threads} -o $bigwig
		plotCoverage --bamfiles $ibam --numberOfProcessors ${params.threads} -o $covPlot
		plotFingerprint --bamfiles $ibam --labels $ibam --numberOfProcessors ${params.threads} --minMappingQuality 30 --skipZeros --plotFile $fingerprintPlot --outRawCounts $fingerprintData
		"""
  }

//Run sam stats
process samStats {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300'
	echo true
	cpus 2
	memory 16

	// KB 2018-06-04: Biowulf Centos7 update (June 04 2018)
	//module 'samtools/1.5'
	module 'samtools/1.8'
	//END 180604

	time '6h'

	tag { ibam }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  //file(ibam) from mergeITRBAMss
	  set file(ibam),file(ibamidx) from mergeBAMwithIDXss

	output:
	  file '*stats.tab' into dtSamStat
	  val 'OK' into ssDone

	script:
		iStat = ibam.name.replaceAll(/.bam/,".idxstats.tab")
		sStat = ibam.name.replaceAll(/.bam/,".samstats.tab")
		//samtools index $ibam
		"""
		samtools idxstats $ibam >$iStat
		samtools stats $ibam >$sStat
		"""
  }

//Make FWD/REV bigwig
process toFRBigWig {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:300 '

	cpus 2

  time { 6.hour }
  memory { 8.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 2

	module 'samtools/1.8'
	module 'ucsc/365'
	module 'bedtools/2.25.0'
	module 'picard/2.9.2'

	tag { ibam }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  //file(ibam) from mergeITRBAM
	  set file(ibam),file(ibamidx) from mergeBAMwithIDXfr

	output:
	  file '*.bigwig' into frBW
	  val 'OK' into frDone

	script:
		iName = ibam.name.replaceAll(/.bam/,".out")
		"""
		perl \$NXF_PIPEDIR/accessoryFiles/SSDS/scripts/ssDNA_to_bigwigs_FASTER_LOMEM.pl --bam $ibam --g ${params.genome} --o $iName --s 100 --w 1000	-v
		"""
  }

//Fnal process will make a multiQC report
process makeSSreport {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:30 --partition=quick,norm'

	cpus 2
	memory '4g'

	module 'python/2.7'
	module 'bedtools/2.25.0'
	// KB 2018-06-04: Biowulf Centos7 update (June 04 2018)
	//module 'samtools/1.5'
	module 'samtools/1.8'
	//END 180604
	time '2h'

	tag { bam }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	file ssdsBEDs from mergeITRBED.collect()
	file bam      from bamAligned.collect()

	output:
	file '*SSDSreport*' into SSDSreport
	val 'OK' into ssRepDone

	script:
	"""
	export TMPDIR=/lscratch/$SLURM_JOBID
	export HOTSPOTS=\$NXF_PIPEDIR/hotspots
	perl \$NXF_PIPEDIR/accessoryFiles/SSDS/scripts/makeSSMultiQCReport_nextFlow.pl $bam $ssdsBEDs
	"""
  }

//Fnal process will make a multiQC report
process multiQC {
	scratch '/lscratch/$SLURM_JOBID'
	clusterOptions ' --gres=lscratch:30 --partition=quick,norm'

	cpus 2
	memory '4g'
	module 'python/2.7'

	time '30m'

	tag { params.outName }
	publishDir params.outdir,  mode: 'copy', overwrite: false

	input:
	  file SSDSreport

	output:
	  file '*ultiQC*' into multiqcOut
	  val 'OK'        into mqDone

	script:
		"""
			export PYTHONPATH=\$NXF_PIPEDIR/accessoryFiles/SSDS/MultiQC_SSDS_Rev1/lib/python2.7/site-packages/
      \$NXF_PIPEDIR/accessoryFiles/SSDS/MultiQC_SSDS_Rev1/bin/multiqc -m ssds -f -n ${params.outName}.multiQC .
		"""
  }

//Print log message on completion
workflow.onComplete {
	if (workflow.success){
		def newFile = new File("${params.outdir}/${params.outName}.ssds_ok.done")
		newFile.createNewFile()
	}


	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
