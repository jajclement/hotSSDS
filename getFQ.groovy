#!/usr/bin/env nextflow

// Kevin Brick : Rev 1.0 : 06-03-18
// Get FASTQ script
// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "getFQ PIPELINE (Version 1.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/pipelines/pipeIt --pipe getfq \\"
  log.info " "
  log.info "/data/RDCO/code/pipelines/pipeIt --h for help \\ "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run $NXF_PIPEDIR/getFQ.groovy \\"
  log.info " --bam               <bam file>"
  log.info " --sra               <sra names (or comma separated list: align as one file)>"
  log.info " --fq1               <Read1 fastq file>"
  log.info " --fq2               <Read2 fastq file>"
  log.info " --obj               <Sample name regex for RDCO object store>"
  log.info " --genome            <string> \\"
  log.info " --threads           <number of threads: default = 6> \\"
  log.info " --outName           <string default=bam file stem> \\"
  log.info " --sample_name       <string: default = outName> \\"
  log.info " --outdir            <string: default = outName> \\"
  log.info " --outdir_tmp        <string>/tmp \\"
  log.info " --merge             Merge fastqs if multiple FQs are passed (i.e. SRA / obj store) <logical: default = true> \\"
  log.info " --withFQC           Generate fastQC & fastqScreen stats                            <logical: default = true> \\"
  log.info " --allowMultiSample  Allow object store regex to \"Get\" multiple samples           <logical: default = true> \\"
  log.info " -with-trace -with-timeline"
  log.info " "
  log.info "HELP: nextflow run /data/RDCO/code/pipelines/getFQ.groovy --help"
  log.info " "
  log.info "================================================================================================================="
  log.info "Required Arguments:"
  log.info " "
  log.info "          --bam        STRING     BAM FILE" OR
  log.info "          --sra        STRING     SRA ACCESSION(s) ; for multiple, separate with commas" OR
  log.info "          --obj        STRING     Sample name regex to search for in object storage" OR
  log.info "          --fq1        STRING     READ1 FASTQ FILE"
  log.info "          --fq2        STRING     READ2 FASTQ FILE"
  log.info "          --genome          STRING     FASTA FILE (GENOME REFERENCE)"
  log.info " "
  log.info "Output and Tempory directory Arguments"
  log.info " "
  log.info "          --outName      STRING        Output file prefix"
  log.info "          --outdir          DIR        Output directory"
  log.info "          --outdir_tmp      DIR        Tempory directory"
  log.info " "
  log.info "Options"
  log.info " "
  log.info "          --merge              true/false    Merge fastqs if multiple FQs are passed (i.e. SRA / obj store) "
  log.info "          --withFQC            true/false    Run fastQC & fastqScreen  "
  log.info "          --allowMultiSample   true/false    Allow object store regex to \"Get\" multiple samples "
  log.info "================================================================================================================"
  exit 1

}

// Define arg defaults

//number of threads
params.threads = 6

//sample information
params.bam = ""
params.fq1 = ""
params.fq2 = ""
params.sra = ""
params.type = "sr"

params.name = ""
params.gzipoutput = "false"

params.merge            = true
params.withFQC          = true
params.allowMultiSample = true //KB Nov 3 2019: allow multiple retrieve by default

params.splitSz = 10000000

if (params.name){
  outName = "${params.name}.${params.genome}"
}else{
  if (params.outName){
    outName     = params.outName
  }else{
    outName     = ""
  }
}

params.tmpName = "${params.name}.tmpFile.${params.genome}"
params.obj=""
params.sample_name=""
params.mem="16G"
params.debugmode=""

//output and tempory directories
params.outdir = "./${outName}"
params.outdir_tmp = "/tmp"

if (params.bam){
  params.bai = "${params.bam}.bai"
}

//genomeAndindex
params.genome_fasta = "$GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"
params.genome_faIdx = "$GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa.fai"

def inputType
if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq1){inputType = 'fastq'}

//log.info
log.info "===================================================================="
log.info "GET FQ PIPELINE : Get initial fastq for pipes"
log.info "===================================================================="
log.info "ref genome         : ${params.genome}"
log.info "genome fasta       : ${params.genome_fasta}"
log.info "bam                : ${params.bam}"
log.info "sra                : ${params.sra}"
log.info "obj                : ${params.obj}"
log.info "fq1                : ${params.fq1}"
log.info "fq2                : ${params.fq2}"
log.info "name               : ${params.outName}"
log.info "source type        : ${inputType}"
log.info "outdir             : ${params.outdir}"
log.info "temp_dir           : ${params.outdir_tmp}"
log.info "threads            : ${params.threads}"
log.info "mem                : ${params.mem}"
log.info "merge              : ${params.merge}"
log.info "withFQC            : ${params.withFQC}"
log.info "allowMultiSample   : ${params.allowMultiSample}"

// This switch deals with multiple input type to generate a fastq from
// either sra records, fq.gz files from the object store, a bam file,
// or directly from fastq/fastq.gz files
//
// PE/SR need not be specified ... it will be autodetected (or @ least it should be !!)
switch (inputType) {
  case 'sra':
    process getSRAfiles {
      scratch '/lscratch/$SLURM_JOBID'
      clusterOptions ' --gres=lscratch:400 --partition=norm'
      echo true
      cpus 2
      memory '4g'
      //time '24h'

      time { 24.hour }
      errorStrategy { 'retry' }
      maxRetries 4

      module 'samtools/1.8'
      module 'picard/2.9.2'
      module 'sratoolkit/2.9.2'

      tag { params.sra }

      publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

      input:
        //val sraNum from sraChannel

      output:
        file '*.R1.fastq' optional true into initFQ1
        file '*.R2.fastq' optional true into initFQ2
        file '*.fastq'      into allFQ, allFQ2
        file 'is*.expttype' into exptType

      script:
        //def randName = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
        // Sept 6 2019: KB
        // Modified to use fasterq-dump instead (fastq-dump is being deprecated)
        // Also modified to parse SRA line and download all SRA files (must be comma delimited)
        //                #fastq-dump ${params.sra} --split-files
        //                 #mv ${params.sra}_1.fastq ${params.sra}.R1.fastq
        //                #mv ${params.sra}_2.fastq ${params.sra}.R2.fastq
        """
        echo 1 >isSR.expttype

        sra=\$(echo "${params.sra}" | tr "," "\\n")

        for s in \$sra
        do
          fasterq-dump \$s --split-files -o \$s

          if [ -f \$s"_2.fastq" ] ; then
            mv \$s"_1.fastq" \$s".R1.fastq"
            mv \$s"_2.fastq" \$s".R2.fastq"

            rm -f isSR.expttype
            echo 2 >isPE.expttype
          else
            mv \$s \$s".R1.fastq"
            echo 1 >isSR.expttype
          fi
        done
        """
      }
    println('Otherside')
    break

  case 'obj':

    println('Retreiving from object storage ... ')

    process getFromObjectStore {
      scratch '/lscratch/$SLURM_JOBID'
      clusterOptions ' --gres=lscratch:600 --partition=norm'
      echo true
      cpus 2
      memory '4g'
      //time '24h'

      time { 24.hour }
      errorStrategy { 'retry' }
      maxRetries 4

      module 'samtools/1.8'
      module 'picard/2.9.2'
      module 'sratoolkit/2.9.2'

      publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

      input:

      output:
        file '*.R1.fastq' optional true into initFQ1
        file '*.R2.fastq' optional true into initFQ2
        file '*.fastq'      into allFQ, allFQ2
        file 'is*.expttype' into exptType

      script:
        def randName = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}

        if (params.allowMultiSample){
          """
          perl \$NXF_PIPEDIR/lsObj --v RDCO --p ${params.obj} --get --f

          gunzip *fastq.gz

          if [ `ls *.R2.*` ]; then
            echo 2 >isPE.expttype
          else
            echo 1 >isSR.expttype
          fi
          """
        }else{
          """
          perl \$NXF_PIPEDIR/lsObj --v RDCO --p ${params.obj} --get

          gunzip *fastq.gz

          if [ `ls *.R2.*` ]; then
            echo 2 >isPE.expttype
          else
            echo 1 >isSR.expttype
          fi
          """
        }
      }
    break

  case 'bam':
    println('BAM input detected ... ')

    inBAM = file(params.bam)

    process bamToFastq {
      scratch '/lscratch/$SLURM_JOBID'
      clusterOptions ' --gres=lscratch:600 --partition=norm'
      echo true
      cpus 2
      memory '4g'
      //time '24h'

      time { 24.hour }
      errorStrategy { 'retry' }
      maxRetries 4

      module 'samtools/1.8'
      module 'picard/2.9.2'

      tag { bam }

      publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

      input:
        file bam from inBAM

      output:
        file '*.R1.fastq' into initFQ1
        file '*.R2.fastq' optional true into initFQ2
        file '*.fastq'      into allFQ, allFQ2
        file 'is*.expttype' into exptType

      script:
        """
        n=`samtools view -h ${bam} |head -n 100000 |samtools view -f 1 -S /dev/stdin |wc -l`

        if [ \$n -eq 0 ]; then
          echo 1 >isSR.expttype

          java -jar \$PICARDJAR SortSam \
                         I=${bam}  \
                         O=querynameSort.bam \
                         SO=queryname \
                         TMP_DIR=/lscratch/\$SLURM_JOBID \
                         VALIDATION_STRINGENCY=LENIENT 2>b2fq.err

          java -jar \$PICARDJAR SamToFastq I=querynameSort.bam \
                         F=nxf.R1.fastq \
                         TMP_DIR=/lscratch/\$SLURM_JOBID \
                         VALIDATION_STRINGENCY=LENIENT
        else
          echo 2 >isPE.expttype
          java -jar \$PICARDJAR FixMateInformation \
                         I=${bam} \
                         O=fixMate.bam \
                         SORT_ORDER=queryname \
                         TMP_DIR=/lscratch/\$SLURM_JOBID \
                         VALIDATION_STRINGENCY=LENIENT

          java -jar \$PICARDJAR SortSam \
                         I=fixMate.bam  \
                         O=querynameSort.bam \
                         SO=queryname \
                         TMP_DIR=/lscratch/\$SLURM_JOBID \
                         VALIDATION_STRINGENCY=LENIENT 2>b2fq.err

          java -jar \$PICARDJAR SamToFastq I=querynameSort.bam \
                         F=nxf.R1.fastq F2=nxf.R2.fastq \
                         TMP_DIR=/lscratch/\$SLURM_JOBID \
                         VALIDATION_STRINGENCY=LENIENT
          fi
        """
        }
    break

  case 'fastq':

    println('FASTQ input(s) detected ... ')

    if (params.fq2){

      inFQ1 = file(params.fq1)
      inFQ2 = file(params.fq2)

      process initFQtoFQ_PE {
        scratch '/lscratch/$SLURM_JOBID'
        clusterOptions ' --gres=lscratch:600 --partition=norm'
        echo true
        cpus 2
        memory '4g'
        //time '24h'

        time { 24.hour }
        errorStrategy { 'retry' }
        maxRetries 4

        module 'samtools/1.8'
        module 'picard/2.9.2'

        tag { inFQ1 }

        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

        input:
          file inFQ1
          file inFQ2

        output:
          file 'nxf.R1.fastq' into initFQ1
          file 'nxf.R2.fastq' into initFQ2
          file 'nxf.R*.fastq' into allFQ, allFQ2
          file 'is*.expttype' into exptType

        script:
          if (inFQ1 =~ /.gz$/)
            """
            gunzip --stdout $inFQ1 >nxf.R1.fastq
            gunzip --stdout $inFQ2 >nxf.R2.fastq

            echo 2 >isPE.expttype
            """
          else
            """
            ln -s $inFQ1 nxf.R1.fastq
            ln -s $inFQ2 nxf.R2.fastq

            echo 2 >isPE.expttype
            """
        }
    }else{
      inFQ1 = file(params.fq1)

      process initFQtoFQ_SR{
        scratch '/lscratch/$SLURM_JOBID'
        clusterOptions ' --gres=lscratch:600 --partition=norm'
        echo true
        cpus 2
        memory '4g'
        //time '24h'

        time { 24.hour }
        errorStrategy { 'retry' }
        maxRetries 4

        module 'samtools/1.8'
        module 'picard/2.9.2'

        tag { inFQ1 }

        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

        input:
          file inFQ1

        output:
          file '*.R1.fastq'   into initFQ1
          file '*.fastq'      into allFQ, allFQ2
          file 'is*.expttype' into exptType

        script:
          if (inFQ1 =~ /.gz$/)
            """
            gunzip --stdout $inFQ1 >nxf.R1.fastq

            echo 1 >isSR.expttype
            """
          else
            """
            cp $inFQ1 nxf.R1.fastq

            echo 1 >isSR.expttype
            """
      }
    }
    break
}

if (params.withFQC){
  process fastQCall {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:400 --partition=norm'
    echo true
    cpus params.threads
    memory params.mem

    module 'fastqtools/0.8'
    module 'fastqc/0.11.8'
    module 'bwa/0.7.17'

    time { 24.hour }
    errorStrategy { 'retry' }
    maxRetries 4

    tag{$fq}

    //publishDir params.outdir, mode: 'copy', overwrite: false
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*zip"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*html"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*png"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*txt"

    input:
      each file(fq) from allFQ

    output:
      file '*zip'  into fqZip
      file '*html' into fqHtml
      file '*png'  into fqPng
      file '*txt'  into fqTxt
        file '*fastqc*' into fastQC
        file '*_screen*' into fastQScreen

    script:
      """
      head -n 10000000 ${fq} >${fq}.subset.fastq

      fastqc -t ${params.threads} ${fq}.subset.fastq

      mv ${fq}.subset_fastqc.html ${fq}c_report.html
      mv ${fq}.subset_fastqc.zip  ${fq}c_report.zip

        \$NXF_PIPEDIR/fastq_screen_v0.11.4/fastq_screen --threads ${params.threads} --force \
                                                       --aligner bwa ${fq} \
                                                       --conf \$NXF_GENOMES/fastq_screen.conf
      """
  }
}

if (!outName){outName='merge'}

println(params.merge)
if (params.merge){
  process mergeFQ {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:800'
    echo true
    cpus params.threads
    memory params.mem

    module 'fastqtools/0.8'

    time { 24.hour }
    errorStrategy { 'retry' }
    maxRetries 1

    publishDir params.outdir, mode: 'move', overwrite: false

    input:
      file fq from allFQ2.collect()

    output:
      file '*.R?.fastq*' into fastqFinal

    script:
    // Modified to work better with huge files
    // Now, we make soft links if at all possible instead of copying stuff
    // This also speeds things up quite a bit
      """
      R1files=\$(ls *R1*fastq 2> /dev/null | wc -l)
      R2files=\$(ls *R2*fastq 2> /dev/null | wc -l)

      if [ "\$R2files" != "0" ]
      then
        if [ "\$R2files" != "1" ]
        then
          cat *R2*fastq >mergeR2.fastq
          rm *.R2.fastq
        else
          fqR2=`ls *R2*fastq`
          ln -s \$fqR2 mergeR2.fastq
        fi
        fastq-sort --id mergeR2.fastq >${outName}.R2.fastq
        rm mergeR2.fastq
      fi

      if [ "\$R1files" != "1" ]
      then
        cat *.R1.fastq >mergeR1.fastq
        rm  *.R1.fastq
      else
        fqR1=`ls *R1*fastq`
        ln -s \$fqR1 mergeR1.fastq
      fi

      fastq-sort --id mergeR1.fastq >${outName}.R1.fastq
      rm mergeR1.fastq

      if [ "${params.gzipoutput}" == "true" ]
      then
          gzip ${outName}.R1.fastq
          gzip ${outName}.R2.fastq || true
      fi
      """
  }
}else{
  process omitFQ {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:60 --partition=norm'
    echo true
    cpus params.threads
    memory params.mem

    time { 4.hour }
    errorStrategy { 'retry' }
    maxRetries 4

    publishDir params.outdir, mode: 'move', overwrite: false

    input:
      file fq from allFQ2.collect()

    output:
      file '*.R?.fastq*' into fastqFinal

    script:
      """
      echo "AAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaagggh FAIL!!!! You forgot to merge !!!)"
      """
  }
}
