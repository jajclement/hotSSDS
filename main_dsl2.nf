#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
                        hotSSDS pipeline version 2.1 (DSL2)
                        Adapted from Kevin Brick (version 1.8_NF)
                        Pauline Auffret, 2020-2021
                        DSL2 migration, 2026
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
Single-Stranded-DNA-Sequencing (SSDS) Pipeline : Align, Parse and call Peaks in ssDNA
----------------------------------------------------------------------------------------
*/

// Construct help message (option --help)
def helpMessage() {
    log.info"""
=============================================================================
  hotSSDS pipeline version 2.1 : Align, parse and call hotspots from SSDNA
=============================================================================
    Usage:

    nextflow run main.nf -c conf/cluster.config -params-file conf/mm10.json \\
        --inputcsv tests/fastq/input.csv --name "runtest" \\
        --trim_cropR1 36 --trim_cropR2 40 --with_trimgalore \\
        -profile singularity -resume

    Runs with Nextflow DSL2 (requires Nextflow >=22.10.0)

    For full parameter list, see README.md
=================================================================================================
""".stripIndent()
}

// Pipeline parameters information
def paramsSection() {
log.info """
==========================================================================
   hotSSDS pipeline version 2.1 : Align, parse and call hotspots from SSDNA
==========================================================================
Run name                       : ${params.name}
Input file                     : ${params.inputcsv}
Reference genome               : ${params.genome}
Blacklist file                 : ${params.blacklist}
Hotspots file                  : ${params.hotspots}
Sequencer platform             : ${params.sequencer}
Output directory               : ${params.outdir}
Use input control files        : ${params.with_control}
Number of replicates           : ${params.nb_replicates}
Bigwig profile                 : ${params.bigwig_profile}
Run IDR analysis               : ${params.with_idr}
Plot saturation curve          : ${params.satcurve}
Use trimgalore                 : ${params.with_trimgalore}
""".stripIndent()
}

//***************************************************************************//
//                         PROCESS DEFINITIONS                               //
//***************************************************************************//

process get_sif {
    tag "${name}"
    label 'internet_access'
    publishDir "${params.outdir}/01_logs/get_sif", mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(dir), val(name)
    output:
        path '*.log'
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
    if grep -q "ERROR" wget_*.log || grep -q "failed" wget_*.log ;
    then
        echo "${name} file download failed. Check sif_url (${params.sif_url}) parameter."
        exit 1
    fi
    """
}

process check_design {
    tag "${design}"
    label 'process_basic'
    label "python"
    publishDir "${params.outdir}/02_results/qc/design/pipeline_info", mode: params.publishdir_mode, pattern: '*.csv'
    publishDir "${params.outdir}/01_logs/check_design",               mode: params.publishdir_mode, pattern: '*.log'
    input:
        path(design)
        path(check_design_script)
    output:
        path 'design_reads.csv',    emit: reads_csv
        path 'design_controls.csv', emit: controls_csv
        path '*.log'
    script:
    """
    python ${check_design_script} ${design} design_reads.csv design_controls.csv >& ${design}_check_design.log
    nb_line_ctrl_file=`cat design_controls.csv | wc -l`
    if [[ \$nb_line_ctrl_file == "1" && "${params.with_control}" == "true" ]]
    then
        echo "The option --with_control is set to true but there is no control files. Check the input csv file."
        exit 1
    fi
    """
}

process makeScreenConfigFile {
    tag "${params.name}"
    label 'process_basic'
    output:
        path "conf.fqscreen", emit: conf
        path "checkfile.ok",  emit: ok
    script:
    // Build the config file content as a Groovy string, then write it inside bash
    def glist = params.genome2screen
    def confLines = []
    confLines.add("This is a config file for FastQ Screen")
    confLines.add("")
    confLines.add("THREADS ${task.cpus}")
    confLines.add("")
    def present = false
    glist.each { item ->
        if (params.genomes.containsKey(item)) {
            confLines.add("DATABASE  ${params.genomes[item].genome_name}    ${params.genomes[item].genome_fasta}")
        }
        if (item == params.genome_name) {
            present = true
        }
    }
    if (!present) {
        confLines.add("DATABASE  ${params.genome_name}    ${params.genome_fasta}")
    }
    def confContent = confLines.join('\n')
    """
    cat <<'FQSCREEN_CONF' > conf.fqscreen
${confContent}
FQSCREEN_CONF
    if [ -f "conf.fqscreen" ]; then
        echo "conf.fqscreen created successfully." > checkfile.ok
    else
        echo "The configuration file for fastqscreen could not be generated. Please check genome2screen parameter."
        exit 1
    fi
    """
}

process crop {
    tag "${sampleId}"
    label 'process_low'
    label 'trimming'
    publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.html"
    publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.zip"
    publishDir "${params.outdir}/02_results/qc/fastqscreen", mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/02_results/qc/raw_fastqc",  mode: params.publishdir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/01_logs/trimming",           mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(sampleId), path(read1), path(read2)
        path(conf_fqscreen)
        path(ok)
    output:
        tuple val(sampleId), path("*_crop_R1.fastq.gz"), path("*_crop_R2.fastq.gz"), emit: cropped
        path('*'), emit: raw_fastqc
    script:
    ext=("${read1.getExtension()}")
    """
    if [[ ${ext} == "gz" ]]
    then
        ext="fastq.gz"
        mv ${read1} ${sampleId}_raw_R1.\${ext}
        mv ${read2} ${sampleId}_raw_R2.\${ext}
    else
        ext="fastq.gz"
        gzip -c ${read1} > ${sampleId}_raw_R1.\${ext}
        gzip -c ${read2} > ${sampleId}_raw_R2.\${ext}
    fi
    zcat ${sampleId}_raw_R1.\${ext} | fastx_trimmer -z -f 1 -l ${params.trim_cropR1} -o ${sampleId}_crop_R1.fastq.gz
    zcat ${sampleId}_raw_R2.\${ext} | fastx_trimmer -z -f 1 -l ${params.trim_cropR2} -o ${sampleId}_crop_R2.fastq.gz
    fastqc -t ${task.cpus} ${sampleId}_crop_R1.fastq.gz ${sampleId}_crop_R2.fastq.gz >& ${sampleId}_raw-fastqc.log 2>&1
    """
}

process trimming {
    tag "${sampleId}"
    label 'process_low'
    label 'trimming'
    publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*_report.txt"
    publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*trim*.html"
    publishDir "${params.outdir}/02_results/qc/trim_fastqc",       mode: params.publishdir_mode, pattern: "*trim*.zip"
    publishDir "${params.outdir}/02_results/trimming/trim_fastq",  mode: params.publishdir_mode, pattern: "*crop_trim_R1.fastq.gz"
    publishDir "${params.outdir}/02_results/trimming/trim_fastq",  mode: params.publishdir_mode, pattern: "*crop_trim_R2.fastq.gz"
    publishDir "${params.outdir}/01_logs/trimming",                mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(sampleId), path(cropread1), path(cropread2)
        path(trimgalore_adapters)
        path(trimmomatic_adapters)
    output:
        tuple val(sampleId), path('*crop_trim_R1.fastq.gz'), path('*crop_trim_R2.fastq.gz'), emit: trimmed
        path('*'), emit: trim_fastqc
    script:
    def global_name = "${params.name}"
    // NovaSeq 2-channel chemistry miscalls dark cycles as high-quality G runs; --nextseq tells
    // cutadapt (via trim_galore) to ignore quality at trailing G runs instead of --quality
    // (see VALIDATION_bwa_mem_vs_original.md, confirmed on real NovaSeq FastQC output)
    def quality_opt = (params.sequencer == 'novaseq') ? "--nextseq ${params.trimg_quality}" : "--quality ${params.trimg_quality}"
    if (params.with_trimgalore && params.trimgalore_adapters)
    """
    trim_galore ${quality_opt} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
                --cores ${task.cpus} --adapter "file:${trimgalore_adapters}" --gzip --paired --basename ${sampleId} ${cropread1} ${cropread2} >& ${sampleId}_trimgalore.log 2>&1
    mv ${sampleId}_val_1.fq.gz ${sampleId}_crop_trim_R1.fastq.gz
    mv ${sampleId}_val_2.fq.gz ${sampleId}_crop_trim_R2.fastq.gz
    fastqc -t ${task.cpus} ${sampleId}_crop_trim_R1.fastq.gz ${sampleId}_crop_trim_R2.fastq.gz >& ${sampleId}_trim-fastqc.log 2>&1
    """
    else if (params.with_trimgalore && !params.trimgalore_adapters)
    """
    trim_galore ${quality_opt} --stringency ${params.trimg_stringency} --length ${params.trim_minlen} \
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

process bwaAlign {
    tag "${sampleId}"
    label 'bam'
    publishDir "${params.outdir}/02_results/bwa/bam",      mode: params.publishdir_mode, pattern: "*.sorted.bam*"
    publishDir "${params.outdir}/02_results/qc/flagstat",  mode: params.publishdir_mode, pattern: "*.flagstat"
    publishDir "${params.outdir}/01_logs/bwa",             mode: params.publishdir_mode, pattern: "*.log"
    publishDir "${params.outdir}/01_logs/bwa",             mode: params.publishdir_mode, pattern: "*.out"
    input:
        tuple val(sampleId), path(fqR1), path(fqR2)
        path(bwa)
        path(fasta)
        path(amb)
        path(ann)
        path(bwt)
        path(pac)
        path(sa)
        path(fai)
    output:
        tuple val(sampleId), path('*.sorted.bam'), emit: sorted_bam
        path('*.flagstat'), emit: bwa_flagstat
        path('*.log')
        path('*.out')
    script:
    def tmp_name = "${params.name}.tmpFile"
    """
    # Align R1+R2 pair with bwa mem (soft-clips the non-homologous ITR part at the 5' of R2,
    # replacing bwa aln / bwa-ra aln / bwa sampe -- see bin/ITR_id_bwamem.pl for ITR parsing)
    ./${bwa} mem -M -t ${task.cpus} ${fasta} ${fqR1} ${fqR2} \
            >${tmp_name}.unsorted.sam 2> ${tmp_name}.unsorted.sam.log
    picard SamFormatConverter I=${tmp_name}.unsorted.sam O=${tmp_name}.unsorted.tmpbam \
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${tmp_name}.unsorted.picardSFC.out 2>&1
    picard SortSam I=${tmp_name}.unsorted.tmpbam O=${sampleId}.sorted.bam SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${sampleId}.sorted.picardSS.out 2>&1
    samtools index ${sampleId}.sorted.bam
    samtools flagstat ${sampleId}.sorted.bam > ${sampleId}.sorted.bam.flagstat
    if grep -q "^0 + 0 mapped" ${sampleId}.sorted.bam.flagstat
    then
        echo "Bam file ${sampleId}.sorted.bam is empty. Check genome parameter."
        exit 1
    fi
    picard -Xms8g -Xmx8g MarkDuplicatesWithMateCigar I=${sampleId}.sorted.bam O=${sampleId}.sorted.md.bam \
        PG=Picard2.9.2_MarkDuplicates M=${sampleId}.sorted.md.txt MINIMUM_DISTANCE=${params.picard_min_distance} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.picard_optdup_distance} REMOVE_DUPLICATES=false \
        CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${sampleId}.sorted.md.out 2>&1
    samtools flagstat ${sampleId}.sorted.md.bam > ${sampleId}.sorted.md.bam.flagstat
    if [[ ${params.with_multimap} == "false" ]]
    then
        # bwa mem reports MAPQ 0 for multi-mapped reads (no XT:A:R tag, that was specific to bwa aln)
        samtools view -h -q 1 ${sampleId}.sorted.bam | samtools view -b | samtools sort -o ${sampleId}.final.bam
        samtools index ${sampleId}.final.bam
        rm ${sampleId}.sorted.bam ; mv ${sampleId}.final.bam ${sampleId}.sorted.bam
        rm ${sampleId}.sorted.bam.bai ; mv ${sampleId}.final.bam.bai ${sampleId}.sorted.bam.bai
    fi
    """
}

process filterBam {
    tag "${sampleId}"
    label 'process_medium'
    label 'bam'
    publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/bed",  mode: params.publishdir_mode, pattern: "*.bed*"
    publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/bam",  mode: params.publishdir_mode, pattern: "*.bam*"
    publishDir "${params.outdir}/01_logs/filterbam/flag_${params.filtering_flag}/log",         mode: params.publishdir_mode, pattern: "*.out"
    input:
        tuple val(sampleId), path(bam)
    output:
        tuple val(sampleId), path('*.unparsed.bam'),                      emit: filtered_bam
        path('*.unparsed.bam'),                                           emit: unparsed_bam
        path('*.unparsed.bam.bai'),                                       emit: unparsed_bai
        tuple val(sampleId), path('*.unparsed.bed'),                      emit: input_bed
        tuple val(sampleId), path('*.unparsed.suppAlignments.bam'),       emit: supp_bam, optional: true
        tuple val(sampleId), path('*.unparsed.suppAlignments.bam.bai'),   emit: supp_bai, optional: true
        path('*.txt')
        path('*.out')
    script:
    """
    samtools view -F ${params.filtering_flag} -f 2 -q ${params.bed_trimqual} -hb ${bam} > ${bam.baseName}.ok.bam
    picard -Xms8g -Xmx8g MarkDuplicatesWithMateCigar I=${bam.baseName}.ok.bam O=${bam.baseName}.unparsed.bam \
        PG=Picard2.9.2_MarkDuplicates M=${bam.baseName}.MDmetrics.txt MINIMUM_DISTANCE=${params.picard_min_distance} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.picard_optdup_distance} REMOVE_DUPLICATES=true \
        CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=. >& ${bam.baseName}.unparsed.picardMD.out 2>&1
    samtools index ${bam.baseName}.unparsed.bam
    bedtools bamtobed -i ${bam.baseName}.unparsed.bam > ${bam.baseName}.unparsed.bed
    if ${params.get_supp} :
    then
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
    publishDir "${params.outdir}/02_results/qc/flagstat",                                                           mode: params.publishdir_mode, pattern: "*.flagstat"
    publishDir "${params.outdir}/01_logs/parseITRs/flag_${params.filtering_flag}",                                  mode: params.publishdir_mode, pattern: "*.out"
    publishDir "${params.outdir}/02_results/bwa/filterbam/flag_${params.filtering_flag}/parse_itr/norm_factors",    mode: params.publishdir_mode, pattern: "*norm_factors.txt"
    input:
        tuple val(sampleId), path(bam)
        path(fasta)
        path(ITR_id_bwamem_script)
    output:
        tuple val(sampleId), path("${bam}"), path('*.ssDNA_type1.bed'), path('*.ssDNA_type2.bed'), \
            path('*.dsDNA.bed'), path('*.dsDNA_strict.bed'), path('*.unclassified.bed'),           emit: itr_bed
        tuple val(sampleId), path("${bam}"), path("${bam}.bai"), path('*.f.ssDNA_type1.bam'), path('*.f.ssDNA_type1.bam.bai'), path('*.f.ssDNA_type2.bam'), \
            path('*.f.ssDNA_type2.bam.bai'), path('*.f.dsDNA.bam'), path('*.f.dsDNA.bam.bai'), path('*.f.dsDNA_strict.bam'), \
            path('*.f.dsDNA_strict.bam.bai'), path('*.f.unclassified.bam'), path('*.f.unclassified.bam.bai'),  emit: itr_bam
        tuple val(sampleId), path('*ssDNA_type1.bed'),                                            emit: t1_bed
        tuple val(sampleId), path('*.f.*.bam'), path('*.f.*.bam.bai'),                            emit: bam_with_idx
        tuple val(sampleId), path('*.f.ssDNA_type1.bam'), path('*.f.ssDNA_type1.bam.bai'),        emit: t1_bam_with_idx
        tuple val(sampleId), path('*norm_factors.txt'), path('*.ssDNA_type1.bed'), \
                path('*.ssDNA_type12.bed'),                                                       emit: bed_to_bw
        path('*.flagstat'),                                                                       emit: itr_flagstat
        path('*.out')
    script:
    """
    perl ${ITR_id_bwamem_script} ${bam} ${fasta} >& ${sampleId}_parseITR.log 2>&1
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
    libsizeT1=`wc -l  *.ssDNA_type1.bed | grep -v "_random"| awk '{print \$1}'`
    libsizeT2=`wc -l  *.ssDNA_type2.bed | grep -v "_random"| awk '{print \$1}'`
    libsizeDS=`wc -l  *.dsDNA.bed | grep -v "_random"| awk '{print \$1}'`
    libsizeDSstrict=`wc -l  *.dsDNA_strict.bed | grep -v "_random"| awk '{print \$1}'`
    libsizeUN=`wc -l  *.unclassified.bed | grep -v "_random"| awk '{print \$1}'`
    libsizeTotal=`expr \$libsizeT1 + \$libsizeT2 + \$libsizeDS + \$libsizeUN`
    echo \$libsizeTotal > ${sampleId}_norm_factors.txt
    cat ${bam}.ssDNA_type1.bed ${bam}.ssDNA_type2.bed | sort -k1,1 -k2,2n > ${bam}.ssDNA_type12.bed
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
    samtools index ${bam}
    samtools index ${bam}.f.ssDNA_type1.bam
    samtools index ${bam}.f.ssDNA_type2.bam
    samtools index ${bam}.f.dsDNA.bam
    samtools index ${bam}.f.dsDNA_strict.bam
    samtools index ${bam}.f.unclassified.bam
    rm -rf *.US.bam
    rm -rf *.sam
    samtools flagstat ${bam}.f.ssDNA_type1.bam > ${bam}.f.ssDNA_type1.bam.flagstat
    if grep -q "^0 + 0 mapped" ${bam}.f.ssDNA_type1.bam.flagstat
    then
        echo "Type 1 bam file ${bam}.f.ssDNA_type1.bam is empty. Check genome parameter."
        exit 1
    fi
    """
}

process makeBigwig {
    tag "${sampleId}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",  mode: params.publishdir_mode, pattern: "*.bw"
    publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",  mode: params.publishdir_mode, pattern: "*.bedgraph"
    publishDir "${params.outdir}/01_logs/bigwig/${params.bigwig_profile}",     mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(sampleId), path(libsizeTotal), path(ssDNA_type1_bed), path(ssDNA_type12_bed)
        path(chrsize)
    output:
        path('*.bw')
        path('*.bedgraph')
        path('*.log')
    script:
    """
    libsizeTot=`cat ${libsizeTotal}`
    libsizeT1=`wc -l ${ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
    normfactT1=`python -c "print(round((1000000/(\$libsizeT1)),5))"`
    normfactTot=`python -c "print(round((1000000/(\$libsizeTot)),5))"`
    genomeCoverageBed -i ${ssDNA_type1_bed} -g ${chrsize} -scale \$normfactT1 -bga > ${sampleId}_ssDNA_type1_RPM-T1.bedgraph
    sort -k1,1 -k2,2n ${sampleId}_ssDNA_type1_RPM-T1.bedgraph > ${sampleId}_ssDNA_type1_RPM-T1_sorted.bedgraph
    bedGraphToBigWig ${sampleId}_ssDNA_type1_RPM-T1_sorted.bedgraph ${chrsize} \
        ${sampleId}_ssDNA_type1_RPM-T1.bw >& ${sampleId}_ssDNA_type1.bedGraphToBigWig.normfactT1.log 2>&1
    genomeCoverageBed -i ${ssDNA_type1_bed} -g ${chrsize} -scale \$normfactTot -bga > ${sampleId}_ssDNA_type1_RPM-Tot.bedgraph
    sort -k1,1 -k2,2n ${sampleId}_ssDNA_type1_RPM-Tot.bedgraph > ${sampleId}_ssDNA_type1_RPM-Tot_sorted.bedgraph
    bedGraphToBigWig ${sampleId}_ssDNA_type1_RPM-Tot_sorted.bedgraph ${chrsize} \
        ${sampleId}_ssDNA_type1_RPM-Tot.bw >& ${sampleId}_ssDNA_type1.bedGraphToBigWig.normfactTot.log 2>&1
    if [[ ${params.bigwig_profile} == "T12" || ${params.bigwig_profile} == "T12rep" ]]
    then
        libsizeT12=`wc -l ${ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
        normfactT12=`python -c "print(round((1000000/(\$libsizeT12)),5))"`
        genomeCoverageBed -i ${ssDNA_type12_bed} -g ${chrsize} -scale \$normfactT12 -bga > ${sampleId}_ssDNA_type12_RPM-T12.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_ssDNA_type12_RPM-T12.bedgraph > ${sampleId}_ssDNA_type12_RPM-T12_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_ssDNA_type12_RPM-T12_sorted.bedgraph ${chrsize} \
            ${sampleId}_ssDNA_type12_RPM-T12.bw >& ${sampleId}_ssDNA_type12.bedGraphToBigWig.normfactT12.log 2>&1
        genomeCoverageBed -i ${ssDNA_type12_bed} -g ${chrsize} -scale \$normfactTot -bga > ${sampleId}_ssDNA_type12_RPM-Tot.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_ssDNA_type12_RPM-Tot.bedgraph > ${sampleId}_ssDNA_type12_RPM-Tot_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_ssDNA_type12_RPM-Tot_sorted.bedgraph ${chrsize} \
            ${sampleId}_ssDNA_type12_RPM-Tot.bw >& ${sampleId}_ssDNA_type12.bedGraphToBigWig.normfactTot.log 2>&1
    fi
    """
}

process makeBigwigReplicates {
    tag "${sampleId}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",  mode: params.publishdir_mode, pattern: "*.bw"
    publishDir "${params.outdir}/02_results/bigwig/${params.bigwig_profile}",  mode: params.publishdir_mode, pattern: "*.bedgraph"
    publishDir "${params.outdir}/01_logs/bigwig/${params.bigwig_profile}",     mode: params.publishdir_mode, pattern: "*.log"
    input:
        tuple val(sampleId), path(R1_libsizeTotal), path(R2_libsizeTotal), path(R1_ssDNA_type1_bed), path(R2_ssDNA_type1_bed), \
            path(R1_ssDNA_type12_bed), path(R2_ssDNA_type12_bed)
    output:
        path('*.bw')
        path('*.bedgraph')
        path('*.log')
    script:
    """
    R1libsizeTot=`cat ${R1_libsizeTotal}`
    R2libsizeTot=`cat ${R2_libsizeTotal}`
    R1libsizeT1=`wc -l ${R1_ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
    R2libsizeT1=`wc -l ${R2_ssDNA_type1_bed} | grep -v "_random"| awk '{print \$1}'`
    R1R2normfactT1=`python -c "print(round(1000000/(\$R1libsizeT1+\$R2libsizeT1),5))"`
    R1R2normfactTot=`python -c "print(round(1000000/(\$R1libsizeTot+\$R2libsizeTot),5))"`
    cat ${R1_ssDNA_type1_bed} ${R2_ssDNA_type1_bed} | grep -v "_random"| sort -k1,1 -k2,2n > ${sampleId}_R1R2_ssDNA_type1.bed
    genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type1.bed -g ${params.chrsize} -scale \$R1R2normfactT1 \
        -bga > ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bedgraph
    sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bedgraph > ${sampleId}_R1R2_ssDNA_type1_RPM-T1_sorted.bedgraph
    bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type1_RPM-T1_sorted.bedgraph \
        ${params.chrsize} ${sampleId}_R1R2_ssDNA_type1_RPM-T1.bw >& ${sampleId}_R1R2_ssDNA_type1.bedGraphToBigWig.normfactT1.log 2>&1
    genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type1.bed -g ${params.chrsize} -scale \$R1R2normfactTot \
        -bga > ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bedgraph
    sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bedgraph > ${sampleId}_R1R2_ssDNA_type1_RPM-Tot_sorted.bedgraph
    bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type1_RPM-Tot_sorted.bedgraph \
        ${params.chrsize} ${sampleId}_R1R2_ssDNA_type1_RPM-Tot.bw >& ${sampleId}_R1R2_ssDNA_type1.bedGraphToBigWig.normfactTot.log 2>&1
    if [[ ${params.bigwig_profile} == "T12rep" ]]
    then
        R1libsizeT12=`wc -l ${R1_ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
        R2libsizeT12=`wc -l ${R2_ssDNA_type12_bed} | grep -v "_random"| awk '{print \$1}'`
        R1R2normfactT12=`python -c "print(round(1000000/(\$R1libsizeT12+\$R2libsizeT12),5))"`
        cat ${R1_ssDNA_type12_bed} ${R2_ssDNA_type12_bed} | grep -v "_random" | sort -k1,1 -k2,2n > ${sampleId}_R1R2_ssDNA_type12.bed
        genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type12.bed -g ${params.chrsize} -scale \$R1R2normfactT12 \
            -bga > ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bedgraph > ${sampleId}_R1R2_ssDNA_type12_RPM-T12_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type12_RPM-T12_sorted.bedgraph ${params.chrsize} \
            ${sampleId}_R1R2_ssDNA_type12_RPM-T12.bw >& ${sampleId}_R1R2_ssDNA_type12.bedGraphToBigWig.normfactT12.log 2>&1
        genomeCoverageBed -i ${sampleId}_R1R2_ssDNA_type12.bed -g ${params.chrsize} -scale \$R1R2normfactTot \
            -bga > ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bedgraph
        sort -k1,1 -k2,2n ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bedgraph > ${sampleId}_R1R2_ssDNA_type12_RPM-Tot_sorted.bedgraph
        bedGraphToBigWig ${sampleId}_R1R2_ssDNA_type12_RPM-Tot_sorted.bedgraph ${params.chrsize} \
            ${sampleId}_R1R2_ssDNA_type12_RPM-Tot.bw >& ${sampleId}_R1R2_ssDNA_type12.bedGraphToBigWig.normfactTot.log 2>&1
    fi
    """
}

process makeDeeptoolsBigWig {
    tag "${sampleId}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}/plot",  mode: params.publishdir_mode, pattern: "*.png"
    publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}",       mode: params.publishdir_mode, pattern: "*.bigwig"
    publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}/tab",   mode: params.publishdir_mode, pattern: "*.tab"
    publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/deeptools/binsize${params.binsize}",          mode: params.publishdir_mode, pattern: "*log*"
    input:
        tuple val(sampleId), path(bam), path(bamidx)
    output:
        path('*')
        path('*.log')
    script:
    """
    bamCoverage --bam ${bam} --normalizeUsing RPKM --binSize ${params.binsize} \
        --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.RPKM.bigwig >& ${sampleId}_makeDeeptoolsBigWig_bamCoverage.log 2>&1
    plotCoverage --bamfiles ${bam} --numberOfProcessors ${task.cpus} -o ${bam.baseName}.deeptools.coveragePlot.png >& ${sampleId}_makeDeeptoolsBigWig_plotCoverage.log 2>&1
    plotFingerprint --bamfiles ${bam} --labels ${bam} --numberOfProcessors ${task.cpus} \
        --minMappingQuality 30 --skipZeros --plotFile ${bam.baseName}.deeptools.fingerprints.png \
        --outRawCounts ${bam.baseName}.deeptools.fingerprints.tab >& ${sampleId}_makeDeeptoolsBigWig_plotFingerprint.log 2>&1
    """
}

process toFRBigWig {
    tag "${sampleId}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/FRbigwig",    mode: params.publishdir_mode, pattern: '*.out'
    publishDir "${params.outdir}/01_logs/bigwig/kbrick_bigwig/FRbigwig",    mode: params.publishdir_mode, pattern: '*.log'
    publishDir "${params.outdir}/02_results/bigwig/kbrick_bigwig/FRbigwig", mode: params.publishdir_mode, pattern: '*.bigwig'
    input:
        tuple val(sampleId), path(bam), path(bamidx)
        path(ssDNA_to_bigwigs_FASTER_LOMEM_script)
        path(fasta)
        path(fai)
    output:
        path('*')
        path('*.log')
    script:
    """
    perl ${ssDNA_to_bigwigs_FASTER_LOMEM_script} --bam ${bam} --g ${fasta} --o ${bam.baseName}.out \
        --s 100 --w 1000 --sc "." --gIdx ${fai} -v >& ${sampleId}_toFRBigwig.log 2>&1
    """
}

// shufBEDs_ct: with-control version
process shufBEDs_ct {
    tag "${id_ip}"
    label 'process_medium'
    label 'peakcalling'
    publishDir "${params.outdir}/02_results/peaks/with-input/bed_shuffle/trim_q${params.bed_trimqual}", mode: params.publishdir_mode
    input:
        tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed)
    output:
        tuple val(id_ip), path("*.IP.sq30.bed"), path("*.CT.sq30.bed"), emit: shuf_beds
    script:
    def ip_Q30_bed      = ip_bed.name.replaceFirst(".bed",".IP.q30.bed")
    def ip_Q30_shuf_bed = ip_bed.name.replaceFirst(".bed",".IP.sq30.bed")
    def ct_Q30_bed      = ct_bed.name.replaceFirst(".bed",".CT.q30.bed")
    def ct_Q30_shuf_bed = ct_bed.name.replaceFirst(".bed",".CT.sq30.bed")
    """
    # Quality trimming IP
    perl -lane '@F = split(/\\t/,\$_); @Q = split(/_/,\$F[3]); print join("\\t",@F) \
        if (\$Q[0] >= ${params.bed_trimqual} && \$Q[1] >= ${params.bed_trimqual})' ${ip_bed} >${ip_Q30_bed}
    # Quality trimming CT
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

// shufBEDs: without-control version
process shufBEDs {
    tag "${id_ip}"
    label 'process_medium'
    label 'peakcalling'
    publishDir "${params.outdir}/02_results/peaks/without-input/bed_shuffle/trim_q${params.bed_trimqual}", mode: params.publishdir_mode
    input:
        tuple val(id_ip), path(ip_bed)
    output:
        tuple val(id_ip), path("*.IP.sq30.bed"), emit: ip_shuf
    script:
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

// Unified callPeaks process: handles both with-control and without-control cases
process callPeaks {
    tag "${id_ip}"
    label 'process_medium'
    label 'peakcalling'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/saturation_curve/${params.sctype}", mode: params.publishdir_mode, pattern: "*peaks_sc.bed"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/macs2/pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}/bed", mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks_sc.bed"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/macs2/pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}/bed", mode: params.publishdir_mode, pattern: "*summits.bed"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/macs2/pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}/xls", mode: params.publishdir_mode, pattern: "*peaks*.xls"
    publishDir "${params.outdir}/01_logs/peaks/${params.with_control ? 'with-input' : 'without-input'}/macs2/pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}", mode: params.publishdir_mode, pattern: "*.macs2.log"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/finalpeaks", mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.bed"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/finalpeaks", mode: params.publishdir_mode, pattern: "*finalPeaks_noIDR.xls"
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/macs2/pv${params.macs_pv}_qv${params.macs_qv}_bw${params.macs_bw}_sloc${params.macs_slocal}_extsize${params.macs_extsize}/narrowPeak", mode: params.publishdir_mode, pattern: "*1.00pc.0_peaks.narrowPeak"
    input:
        tuple val(id_ip), path(ip_bed), val(shuffle_percent)
        path(ct_bed)
        path(blacklist)
        path(fai)
        path(pickNlines_script)
        val(satCurveReps)
    output:
        path("*peaks*.bed"),                                              emit: allbed
        path("*peaks_sc.bed"),                                            emit: satcurvebed, optional: true
        path("*peaks*.xls"),                                              optional: true
        path("*peaks*.bedgraph"),                                         optional: true
        path("*.macs2.log")
        tuple val(id_ip), path(ip_bed), path("*peaks.bed"),              emit: all_peaks, optional: true
        path("genome_size.txt"),                                          emit: genome_size
        tuple val(id_ip), path("*finalPeaks_noIDR.bed"),                  emit: final_peaks_noidr, optional: true
        path("*")
    script:
    def ct_flag = params.with_control ? "-c ${ct_bed}" : ""
    """
    ## SELECT N LINES FROM IP BED FILE ACCORDING TO satCurve parameter
    nT=`cat ${ip_bed} |wc -l`
    nPC=`perl -e 'print int('\$nT'*${shuffle_percent})'`
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

    # Write genome_size to file for downstream processes (resume-safe)
    echo -n \$genome_size > genome_size.txt

    ## CALL PEAKS WITH MACS2 N TIMES
    for i in {0..${satCurveReps}}; do
        peakFileName=${id_ip}'.N'\$nPC'_${shuffle_percent}pc.'\$i
        if [ ${params.macs_pv} != -1 ]; then
            macs2 callpeak -g \$genome_size -t \$nPC.IP.bed ${ct_flag} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$peakFileName --nomodel --extsize ${params.macs_extsize} --pvalue ${params.macs_pv} > IP-${id_ip}.macs2.log 2>&1
        else
            macs2 callpeak -g \$genome_size -t \$nPC.IP.bed ${ct_flag} \
                    --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} \
                    --name \$peakFileName --nomodel --extsize ${params.macs_extsize} --qvalue ${params.macs_qv} > IP-${id_ip}.macs2.log 2>&1
        fi

        # Filter out peaks from blacklist bed file
        if [ ${params.blacklist} != "None" ]
        then
            intersectBed -a \$peakFileName'_peaks.narrowPeak' -b ${blacklist} -v >\$peakFileName'.peaks_sc.noBL'
            cut -f1-3 \$peakFileName'.peaks_sc.noBL' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
        else
            cut -f1-3 \$peakFileName'_peaks.narrowPeak' |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$peakFileName'_peaks_sc.noM'
        fi
        if [ ${params.no_chrY} ]
        then
            grep -v 'chrY' \$peakFileName'_peaks_sc.noM' |sort -k1,1 -k2n,2n > \$peakFileName'_peaks_sc.bed'
        else
            mv \$peakFileName'_peaks_sc.noM' \$peakFileName'_peaks_sc.bed'
        fi
        mv \$peakFileName'_peaks.xls' \$peakFileName'_peaks_sc.xls'
    done

    ##MERGE BED FILES FROM THE LOOP
    sort -k1,1 -k2n,2n -k3n,3n ${id_ip}*peaks_sc.bed |mergeBed -i - >${id_ip}.${shuffle_percent}.peaks_sc.bed

    ##POSTPROCESS IF shuffle_percent is 100%
    if [ ${shuffle_percent} == 1.00 ]; then
        mv ${id_ip}.${shuffle_percent}.peaks_sc.bed ${id_ip}.peaks.bed
        cp ${id_ip}.peaks.bed ${id_ip}.finalPeaks_noIDR.bed
        cat *1.00pc.0_peaks_sc.xls >${id_ip}.peaks.xls
        cp ${id_ip}.peaks.xls ${id_ip}.finalPeaks_noIDR.xls
    fi
    rm *.tmp
    """
}

process samStats {
    tag "${sampleId}"
    label 'process_basic'
    label 'ssds'
    publishDir "${params.outdir}/02_results/qc/samstats", mode: params.publishdir_mode, pattern: "*.tab"
    input:
        tuple val(sampleId), path(bam), path(bamidx)
    output:
        path '*stats.tab', emit: samstats
    script:
    """
    samtools idxstats ${bam} > ${bam.baseName}.idxstats.tab
    samtools stats ${bam} > ${bam.baseName}.samstats.tab
    """
}

process makeSSreport {
    tag "${sampleId}"
    label 'process_basic'
    label 'ssds'
    publishDir "${params.outdir}/02_results/qc/multiqc",  mode: params.publishdir_mode, pattern: '*SSDSreport*'
    publishDir "${params.outdir}/01_logs/makeSSreport",    mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(sampleId), path(bam), path(T1), path(T2), path(ds), path(dss), path(unclassified)
        path(hotspots)
        path(makeSSMultiQCReport_nextFlow_script)
    output:
        tuple val(sampleId), path('*SSDSreport*'), emit: ssds_report
        path('*SSDSreport*'),                      emit: ssds_report_flat
        path('*.log')
    script:
    """
    perl ${makeSSMultiQCReport_nextFlow_script} ${bam} ${T1} ${T2} ${ds} ${dss} ${unclassified} \
        --g ${params.genome} --h ${hotspots} >& ${sampleId}_makeSSreport.log 2>&1
    """
}

process makeFingerPrint {
    tag "${params.name}"
    label 'process_low'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/qc/fingerprint", mode: params.publishdir_mode, pattern: '*.png'
    publishDir "${params.outdir}/02_results/qc/fingerprint", mode: params.publishdir_mode, pattern: '*.tab'
    publishDir "${params.outdir}/01_logs/fingerprint",       mode: params.publishdir_mode, pattern: '*.log'
    input:
        path(all_bam)
        path(all_bai)
    output:
        path('*.png'), emit: png
        path('*.tab'), emit: tab
        path('*.log')
    script:
    def global_name = "${params.name}"
    """
    plotFingerprint --bamfiles ${all_bam} --smartLabels --numberOfProcessors ${task.cpus} \
        --minMappingQuality 30 --skipZeros --plotFile ${global_name}.deeptools.fingerprints.png \
        --outRawCounts ${global_name}.deeptools.fingerprints.rawcounts.tab --outQualityMetrics ${global_name}.deeptools.fingerprints.qual.tab \
        >& ${global_name}_makeFingerPrint.log
    """
}

process computeFRIP {
    tag "${sampleId}"
    label 'process_low'
    label 'frip'
    publishDir "${params.outdir}/02_results/qc/multiqc", mode: params.publishdir_mode, pattern: '*SSDSreport*'
    publishDir "${params.outdir}/02_results/qc/frip",    mode: params.publishdir_mode, pattern: '*'
    publishDir "${params.outdir}/01_logs/frip",           mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(sampleId), path(bam), path(bai), path(T1), path(T1bai), path(T2), path(T2bai), path(ds), path(dsbai), \
            path(dss), path(dssbai), path(unc), path(uncbai)
        path(peaks)
        path(report)
        path(get_frip_script)
    output:
        tuple val(sampleId), path('*SSDSreport*'), emit: ssds_report_denovo
        path('*SSDSreport*'),                      emit: ssds_report_denovo_flat
        path("*")
    script:
    """
    python ${get_frip_script} ${peaks} ${bam} ${task.cpus} total ${params.genome_name} ${sampleId}_total.frip
    python ${get_frip_script} ${peaks} ${T1} ${task.cpus} ssType1 ${params.genome_name} ${sampleId}_type1.frip >& ${sampleId}_computeFRIP_type1.log 2>&1
    python ${get_frip_script} ${peaks} ${T2} ${task.cpus} ssType2 ${params.genome_name} ${sampleId}_type2.frip >& ${sampleId}_computeFRIP_type2.log 2>&1
    python ${get_frip_script} ${peaks} ${ds} ${task.cpus} dsDNA ${params.genome_name} ${sampleId}_ds.frip >& ${sampleId}_computeFRIP_ds.log 2>&1
    python ${get_frip_script} ${peaks} ${dss} ${task.cpus} dsDNA ${params.genome_name} ${sampleId}_dss.frip >& ${sampleId}_computeFRIP_dss.log 2>&1
    python ${get_frip_script} ${peaks} ${unc} ${task.cpus} unclassified ${params.genome_name} ${sampleId}_unc.frip >& ${sampleId}_computeFRIP_unc.log 2>&1
    cat ${report} ${sampleId}_type1.frip ${sampleId}_type2.frip ${sampleId}_ds.frip ${sampleId}_dss.frip ${sampleId}_unc.frip > ${sampleId}_denovo.SSDSreport.tab
    """
}

process ssds_multiqc_report {
    tag "${sampleId}"
    label 'process_basic'
    label 'plot'
    publishDir "${params.outdir}/02_results/qc/multiqc", mode: params.publishdir_mode
    publishDir "${params.outdir}/01_logs/multiqc",        mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(sampleId), path(report)
        path(runPlotSSDSqc_script)
        path(plot_ssds_stat_script)
    output:
        path('ssds/*.png'), emit: multiqc_plots
        path('ssds/*.pdf')
    script:
    """
    mkdir -p ssds
    Rscript ${runPlotSSDSqc_script} ${plot_ssds_stat_script} ${report} ${sampleId} "./ssds" >& ${sampleId}_ssds_multiqc.log
    cd ssds ; for file in \$(ls *.png); do mv \$file \${file%.*}_mqc.png; done ; cd ..
    """
}

process renameInputBed {
    tag "${id_ip}"
    label 'process_basic'
    input:
        tuple val(id_ip), val(id_ct), path(ip_bed), path(ct_bed)
    output:
        tuple val(id_ip), val(id_ct), path(ip_bed), path('*renamed*'), emit: renamed
    script:
    """
    random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
    cat ${ct_bed} > ${ct_bed}_\${random_id}_renamed.bed
    """
}

process createPseudoReplicates {
    tag "${id_ip}"
    label 'process_basic'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/pseudo_replicates", mode: params.publishdir_mode, pattern: '*.bed'
    input:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2)
    output:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path('*r1_pseudorep_r1.bed'), path('*r1_pseudorep_r2.bed'), \
            path('*r2_pseudorep_r1.bed'), path('*r2_pseudorep_r2.bed'), path('*pool_r1.bed'), path('*pool_r2.bed'), \
            path('*poolT.bed'), emit: all_rep
    script:
    """
    nlines_r1=\$((`cat ${t1_rep1} | wc -l`/2))
    nlines_r2=\$((`cat ${t1_rep2} | wc -l`/2))
    shuf ${t1_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
    shuf ${t1_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
    mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
    mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
    mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
    mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed
    nlines_pool=\$((`cat ${t1_rep1} ${t1_rep2} | wc -l`/2))
    cat ${t1_rep1} ${t1_rep2} > ${id_ip}_poolT.bed
    shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
    mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
    mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
    """
}

process createPseudoReplicates_ct {
    tag "${id_ip}"
    label 'process_medium'
    publishDir "${params.outdir}/02_results/idr/with-input/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/pseudo_replicates", mode: params.publishdir_mode, pattern: '*.bed'
    input:
        tuple val(id_ip), val(id_ct), path(t1_rep1), path(t1_rep2), path(ct_r1), path(ct_r2)
    output:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path('*r1_pseudorep_r1.bed'), path('*r1_pseudorep_r2.bed'), \
            path('*r2_pseudorep_r1.bed'), path('*r2_pseudorep_r2.bed'), path('*pool_r1.bed'), path('*pool_r2.bed'), \
            path('*poolT.bed'), path('*ct_pool.bed'), emit: all_rep
    script:
    """
    nlines_r1=\$((`cat ${t1_rep1} | wc -l`/2))
    nlines_r2=\$((`cat ${t1_rep2} | wc -l`/2))
    shuf ${t1_rep1} | split -d -l \$nlines_r1 - ${id_ip}_r1_pseudorep_
    shuf ${t1_rep2} | split -d -l \$nlines_r2 - ${id_ip}_r2_pseudorep_
    mv ${id_ip}_r1_pseudorep_00 ${id_ip}_r1_pseudorep_r1.bed
    mv ${id_ip}_r1_pseudorep_01 ${id_ip}_r1_pseudorep_r2.bed
    mv ${id_ip}_r2_pseudorep_00 ${id_ip}_r2_pseudorep_r1.bed
    mv ${id_ip}_r2_pseudorep_01 ${id_ip}_r2_pseudorep_r2.bed
    nlines_pool=\$((`cat ${t1_rep1} ${t1_rep2} | wc -l`/2))
    cat ${t1_rep1} ${t1_rep2} > ${id_ip}_poolT.bed
    shuf ${id_ip}_poolT.bed | split -d -l \$nlines_pool - ${id_ip}_pool_
    mv ${id_ip}_pool_00 ${id_ip}_pool_r1.bed
    mv ${id_ip}_pool_01 ${id_ip}_pool_r2.bed
    if cmp -s ${ct_r1} ${ct_r2}
    then
        cat ${ct_r1} > ${id_ip}_ct_pool.bed
    else
        cat ${ct_r1} ${ct_r2} | sort -n | uniq > ${id_ip}_ct_pool.bed
    fi
    """
}

process callPeaksForIDR {
    tag "${id_ip}"
    label 'process_medium'
    label 'peakcalling'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/macs2", mode: params.publishdir_mode, pattern: '*narrowPeak*'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/macs2", mode: params.publishdir_mode, pattern: '*regionPeak*'
    publishDir "${params.outdir}/01_logs/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/macs2", mode: params.publishdir_mode, pattern: '*.macs2.log'
    input:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(r1_pseudorep_r1), path(r1_pseudorep_r2), \
            path(r2_pseudorep_r1), path(r2_pseudorep_r2), path(pool_r1), path(pool_r2), path(poolT)
        val(genome_size)
        path(ctpool)
    output:
        path('*narrowPeak*')
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path('*_R1*type1*.regionPeak'), path('*_R2*type1*.regionPeak'), \
            path('*r1_pseudorep_r1*.regionPeak'), path('*r1_pseudorep_r2*.regionPeak'), \
            path('*r2_pseudorep_r1*.regionPeak'), path('*r2_pseudorep_r2*.regionPeak'), \
            path('*pool_r1*.regionPeak'), path('*pool_r2*.regionPeak'), path('*poolT*.regionPeak'), emit: all_peaks_rep
        path('*.macs2.log')
    script:
    def ct_flag = params.with_control ? "-c ${ctpool}" : ""
    def ct_skip = params.with_control ? "ctname=`basename -- ${ctpool} .bed`" : "ctname='__NO_CONTROL__'"
    """
    ${ct_skip}
    for file in \$(ls *.bed) ;
    do
        name=`basename -- \$file .bed`
        if [ \$ctname != \$name ]
        then
            random_id=`shuf -zer -n20  {A..Z} {a..z} {0..9}`
            mkdir \${random_id}
            if [[ ${params.idr_macs_pv} == -1 && ${params.idr_macs_qv} == -1 ]];
            then
                macs2 callpeak -g ${genome_size} -t \$file ${ct_flag} \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                --name \$name --nomodel --extsize ${params.macs_extsize} > \$name.macs2.log 2>&1
            elif [ ${params.idr_macs_pv} != -1 ];
            then
                macs2 callpeak -g ${genome_size} -t \$file ${ct_flag} \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                --name \$name --nomodel --extsize ${params.macs_extsize} --pvalue ${params.idr_macs_pv} > \$name.macs2.log 2>&1
            elif [ ${params.idr_macs_qv} != -1 ];
            then
                macs2 callpeak -g ${genome_size} -t \$file ${ct_flag} \
                --bw ${params.macs_bw} --keep-dup all --slocal ${params.macs_slocal} --tempdir \${random_id} \
                --name \$name --nomodel --extsize ${params.macs_extsize} --qvalue ${params.idr_macs_qv} > \$name.macs2.log 2>&1
            fi
            npeaks=`cat \${name}_peaks.narrowPeak | wc -l`
            if [ \$npeaks -lt 20 ]
            then
                echo "Not enough peaks in \${name}_peaks.narrowPeak file to perform IDR analysis."
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

process IDRanalysis {
    tag "${id_ip}"
    label 'process_medium'
    label 'idr'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/plot",                mode: params.publishdir_mode, pattern: '*.png'
    publishDir "${params.outdir}/01_logs/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}",                         mode: params.publishdir_mode, pattern: '*.log'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/unthresholded-peaks",   mode: params.publishdir_mode, pattern: '*unthresholded-peaks*'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/bfilt",                 mode: params.publishdir_mode, pattern: '*.bfilt.gz'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/peaks",                 mode: params.publishdir_mode, pattern: '*Peak.gz'
    input:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(ip_rep1), path(ip_rep2), path(r1_pseudorep_r1), path(r1_pseudorep_r2), \
            path(r2_pseudorep_r1), path(r2_pseudorep_r2), path(pool_r1), path(pool_r2), path(poolT)
        path(blacklist)
        path(chrsize)
        path(encode_idr_script)
        path(fai)
    output:
        path('*unthresholded-peaks*')
        path('*.log')
        path('*.png')
        path('*.bfilt.gz')
        path('*Peak.gz')
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path('*r1_pseudorep*.*Peak.gz'), path('*r2_pseudorep*.*Peak.gz'), \
            path('*truerep*.*Peak.gz'), path('*poolrep*.*Peak.gz'), emit: idr_peaks
    script:
    """
    if [ ${params.blacklist} != "None" ]
    then
        blacklist_params="--blacklist ${blacklist}"
    else
        blacklist_params=""
    fi
    genome_size=`cut -f3 ${fai} |tail -n1`
    npeaksR1=`cat ${ip_rep1} | wc -l`
    npeaksR2=`cat ${ip_rep2} | wc -l`
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
    find . -name "${id_ip}_r1_pseudorep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_r1_pseudorep.bfilt.gz \\;
    find . -name "${id_ip}_r2_pseudorep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_r2_pseudorep.bfilt.gz \\;
    find . -name "${id_ip}_truerep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {}  ${id_ip}_truerep.bfilt.gz \\;
    find . -name "${id_ip}_poolrep.idr*.bfilt.${params.idr_peaktype}.gz" -exec mv {} ${id_ip}_poolrep.bfilt.gz \\;
    """
}

process IDRpostprocess {
    tag "${id_ip}"
    label 'process_basic'
    label 'idr'
    publishDir "${params.outdir}/01_logs/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/qc",    mode: params.publishdir_mode, pattern:'*.log'
    publishDir "${params.outdir}/02_results/idr/${params.with_control ? 'with-input' : 'without-input'}/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/peaks", mode: params.publishdir_mode, pattern:'*.finalPeaks*'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/finalpeaks/idr/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}", mode: params.publishdir_mode, pattern:'*.finalPeaks*'
    input:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(r1_pseudorep_peaks), path(r2_pseudorep_peaks), \
            path(truerep_peaks), path(poolrep_peaks)
    output:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path('*.finalPeaks_IDR.bed'), emit: idr_final_peaks
        path('*_IDR_QC.log')
        path('*.finalPeaks*')
    script:
    """
    n1=`zcat ${r1_pseudorep_peaks} | wc -l`
    n2=`zcat ${r2_pseudorep_peaks} | wc -l`
    np=`zcat ${poolrep_peaks} | wc -l`
    nt=`zcat ${truerep_peaks} | wc -l`
    if [ \$n1 -gt \$n2 ]
    then
        scr=`python -c "print(\$n1/\$n2)"`
    else
        scr=`python -c "print(\$n2/\$n1)"`
    fi
    echo "Self-consistency Ratio = \$scr" >> ${id_ip}_IDR_QC.log
    if [ \$np -gt \$nt ]
    then
        rr=`python -c "print(\$np/\$nt)"`
        zcat ${poolrep_peaks} > ${poolrep_peaks}.finalPeaks_IDR.bed
    else
        rr=`python -c "print(\$nt/\$np)"`
        zcat ${truerep_peaks} > ${truerep_peaks}.finalPeaks_IDR.bed
    fi
    echo "Rescue Ratio = \$rr" >> ${id_ip}_IDR_QC.log
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

process normalizePeaks {
    tag "${id_ip}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/no-idr",      mode: params.publishdir_mode, pattern: '*.bedgraph'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/no-idr/tab",  mode: params.publishdir_mode, pattern: '*.tab'
    publishDir "${params.outdir}/01_logs/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/no-idr",         mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(id_ip), path(ip_bed), path(peaks_bed)
        path(norm_script)
        path(reverse_script)
    output:
        tuple val(id_ip), path("*normpeaks.bedgraph"), path("*normpeaks.tab")
        path('*.normpeaks.log')
    script:
    """
    perl ${norm_script} --bed ${peaks_bed} \
        --in ${ip_bed} --out ${id_ip}.no-idr.normpeaks.bedgraph \
        --rc --rev_src ${reverse_script} > ${id_ip}.no-idr.normpeaks.log 2>&1
    """
}

process normalizePeaks_idr {
    tag "${id_ip}"
    label 'process_basic'
    label 'bigwig'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/idr/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}",     mode: params.publishdir_mode, pattern: '*.bedgraph'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/idr/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}/tab", mode: params.publishdir_mode, pattern: '*.tab'
    publishDir "${params.outdir}/01_logs/peaks/${params.with_control ? 'with-input' : 'without-input'}/normalized/idr/${params.idr_peaktype}_macs2pv${params.idr_macs_pv}_macs2qv${params.idr_macs_qv}_setup-${params.idr_setup}",        mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(id_ip), path(t1_rep1), path(t1_rep2), path(peaks_bed)
        path(norm_script)
        path(reverse_script)
    output:
        tuple val(id_ip), path("*normpeaks*.bedgraph"), path("*normpeaks*.tab")
        path('*.normpeaks*.log')
    script:
    """
    cat ${t1_rep1} ${t1_rep2} | sort -k1,1 -k2,2n > ${id_ip}_rep1_rep2.bed
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

process mergeFinalPeaks {
    tag "${id_ip}"
    label 'peakcalling'
    label 'process_low'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/finalpeaks/no-idr", mode: params.publishdir_mode, pattern: '*.bed'
    publishDir "${params.outdir}/01_logs/peaks/${params.with_control ? 'with-input' : 'without-input'}/finalpeaks/no-idr",    mode: params.publishdir_mode, pattern: '*.log'
    input:
        tuple val(id_ip), path(R1peaks), path(R2peaks)
    output:
        path('*mergePeaks.bed')
        path('*.log')
    script:
    """
    cat ${R1peaks} ${R2peaks} | sort -k1,1 -k2,2n > ${id_ip}_R1R2_no-idr.mergePeaks.bed
    mergeBed -i ${id_ip}_R1R2_no-idr.mergePeaks.bed >& ${id_ip}_R1R2_no-idr.mergeBed.log 2>&1
    """
}

process makeSatCurve {
    tag "${params.name}"
    label 'process_basic'
    label 'peakcalling'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/saturation_curve/${params.sctype}", mode: params.publishdir_mode, pattern: '*.tab'
    publishDir "${params.outdir}/02_results/peaks/${params.with_control ? 'with-input' : 'without-input'}/saturation_curve/${params.sctype}", mode: params.publishdir_mode, pattern: '*.png'
    publishDir "${params.outdir}/01_logs/peaks/${params.with_control ? 'with-input' : 'without-input'}/saturation_curve/${params.sctype}",    mode: params.publishdir_mode, pattern: '*.log'
    input:
        path(saturation_curve_data)
        path(getPeaksBedFiles_script)
        path(runSatCurve_script)
        path(satCurveHS_script)
    output:
        path("*.tab"), emit: tab
        path("*.png"), emit: png
        path("*.log")
    script:
    def global_name = "${params.name}"
    """
    perl ${getPeaksBedFiles_script} -tf satCurve.tab -dir "." >& ${global_name}_makeSatCurve.log 2>&1
    Rscript ${runSatCurve_script} ${satCurveHS_script} satCurve.tab ${global_name} >& ${global_name}_makeSatCurve_Rscript.log 2>&1
    """
}

process general_multiqc {
    tag "${params.name}"
    label 'process_basic'
    label 'multiqc'
    label 'internet_access'
    publishDir "${params.outdir}/02_results/qc/multiqc", mode: 'copy'
    publishDir "${params.outdir}/01_logs/multiqc",       mode: 'copy', pattern: '*.log'
    input:
        path(logo)
        path(multiqc_configfile)
        path(png_fingerprint)
        path(tab_fingerprint)
        path(bwa_flagstat)
        path(itr_flagstat)
        path(ssdsreport)
        path(ssds_report_denovo)
        path(multiqc_denovo)
        path(multiqc)
        path(samstats)
        path(raw_fastqc)
        path(trim_fastqc)
        path(satcurve)
    output:
        path('*')
    script:
    def global_name = "${params.name}"
    """
    multiqc --export -c ${multiqc_configfile} \
        -n ${global_name}.multiQC.quality-control.report "." \
        >& ${global_name}.multiQC.quality-control.log 2>&1
    """
}

//***************************************************************************//
//                           WORKFLOW                                        //
//***************************************************************************//

// Default values for CLI-only params
params.version = false
params.help = false
params.get_sif = false

// Sequencing platform: 'hiseq' (default) or 'novaseq' -- adjusts optical duplicate marking
// distance and trim_galore's G-tailing handling (validated/adjusted inside workflow{}, since
// NB: params loaded from a JSON -params-file come through as Strings, e.g. "100")
params.sequencer = params.sequencer ?: 'hiseq'
params.picard_optdup_distance = params.picard_optdup_distance ?: 100

// Define genome variables (must be at top-level)
params.genome_fasta = params.genome ? params.genomes[ params.genome ].genome_fasta ?: false : false
params.genomedir = params.genome ? params.genomes[ params.genome ].genomedir ?: false : false
params.genome_name = params.genome ? params.genomes[ params.genome ].genome_name ?: false : false
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false

workflow {

    // Show help / version
    if (params.version){
        println("This is $workflow.manifest.name version $workflow.manifest.version")
        println("$workflow.manifest.description")
        return
    }
    if (params.help){
        helpMessage()
        return
    }

    // Validate parameters
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        exit 1, "The provided genome '${params.genome}' is not available. Available: ${params.genomes.keySet().join(", ")}"
    }
    if(params.bigwig_profile != "T1" && params.bigwig_profile != "T12" && params.bigwig_profile != "T1rep" && params.bigwig_profile != "T12rep") {
        exit 1, "Error : --bigwig_profile must be T1, T12, T1rep or T12rep."
    }
    if(params.idr_setup != "auto" && params.idr_setup != "custom") {
        exit 1, "Error : --idr_setup must be auto or custom."
    }
    if (!['copy','symlink','rellink','link','copyNoFollow','move'].contains(params.publishdir_mode)) {
        exit 1, "Error : --publishdir_mode must be symlink, rellink, link, copy, copyNoFollow, move, see https://www.nextflow.io/docs/latest/process.html"
    }
    if(params.satcurve) {
        if(params.sctype != "minimal" && params.sctype != "standard" && params.sctype != "expanded") {
            exit 1, "Error : --sctype must be minimal, standard or expanded."
        }
        if(params.reps < 3) {
            exit 1, "Error : when --satcurve is true, --reps must be >= 3"
        }
    }
    if(params.with_idr && params.nb_replicates != "2") {
        exit 1, "Error : when --with_idr is true, --nb_replicates must be 2"
    }

    // Determine correct picard optical duplicate distance based on sequencer
    if (params.sequencer == 'novaseq' && params.picard_optdup_distance.toString() == '100') {
        params.picard_optdup_distance = 2500
        log.info "INFO: --sequencer 'novaseq' detected and --picard_optdup_distance is at its default for unpatterned flowcells (100). Adjusting to 2500."
    } else if (params.sequencer == 'novaseq' && params.picard_optdup_distance.toString() != '2500') {
        log.warn "WARNING: --sequencer 'novaseq' was specified, but --picard_optdup_distance is set to ${params.picard_optdup_distance}. For NovaSeq, 2500 is recommended."
    } else if (params.sequencer && !['hiseq', 'novaseq'].contains(params.sequencer)) {
        log.warn "WARNING: Unknown sequencer '${params.sequencer}'. Using provided/default picard_optdup_distance: ${params.picard_optdup_distance}."
    }
    // trimmomatic has no equivalent to trim_galore's --nextseq mode for NovaSeq's 2-channel chemistry
    if (params.sequencer == 'novaseq' && !params.with_trimgalore) {
        log.warn "WARNING: --sequencer 'novaseq' was specified with trimmomatic (not --with_trimgalore). Trimmomatic has no equivalent to trim_galore's --nextseq mode for NovaSeq's 2-channel chemistry G-tailing artifact; consider --with_trimgalore for NovaSeq data."
    }

    // Check bwa index files exist
    file("${params.genome_fasta}.amb", checkIfExists: true)
    file("${params.genome_fasta}.ann", checkIfExists: true)
    file("${params.genome_fasta}.bwt", checkIfExists: true)
    file("${params.genome_fasta}.pac", checkIfExists: true)
    file("${params.genome_fasta}.sa",  checkIfExists: true)

    paramsSection()

    // Saturation curve settings
    def satCurveReps
    def satCurvePCsList
    if (params.satcurve){
        satCurveReps = params.reps - 1
        if (params.sctype == 'expanded'){
            satCurvePCsList = [0.025,0.05,0.075,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]
        } else if (params.sctype == 'standard'){
            satCurvePCsList = [0.20,0.40,0.60,0.80,1.00]
        } else {
            satCurvePCsList = [0.10,0.50,1.00]
        }
    } else {
        satCurvePCsList = [1.00]
        satCurveReps = 0
    }

    //=========================================================================
    // If get_sif mode, just download containers and stop
    //=========================================================================
    if (params.get_sif) {
        sif_ch = channel.from(
            ["bam-box-1.0","bam-box_1.0.sif"],
            ["bigwig-box-1.0","bigwig-box-1.0.sif"],
            ["frip-box-1.0","frip-box_1.0.sif"],
            ["idr-box-1.0","idr-box_1.0.sif"],
            ["multiqc-box-1.0","multiqc-box_1.0.sif"],
            ["peak-calling-box-1.0","peak-calling-box_1.0.sif"],
            ["plot-box-1.0","plot-box_1.0.sif"],
            ["python-3.8","python-3.8.sif"],
            ["ssds-qc-box-1.0","ssds-qc-box_1.0.sif"],
            ["trimming-box-1.0","trimming-box_1.0.sif"]
        )
        get_sif(sif_ch)
        return
    }

    //=========================================================================
    // Input files and scripts
    //=========================================================================

    // If test profile, create input csv
    if (workflow.profile ==~ /.*test.*/) {
        def input_file = new File("${baseDir}/tests/fastq/input.csv")
        input_file.write "group,replicate,fastq_1,fastq_2,antibody,control\n"
        input_file << "TEST_IP,1,${baseDir}/tests/fastq/ssdsLong.100k.R1.fastq,${baseDir}/tests/fastq/ssdsLong.100k.R2.fastq,antiDMC1,"
    }

    // Input csv
    if (!params.inputcsv) { exit 1, 'Samples design file not specified!' }
    def input_ch = file(params.inputcsv, checkIfExists: true)

    // Input files as channels (DSL2: can be consumed multiple times)
    def chrsize_ch    = params.chrsize ? file(params.chrsize, checkIfExists: true) : file("${workDir}/NO_CHRSIZE")
    def hotspots_ch   = (params.hotspots && params.hotspots != "None") ? file(params.hotspots, checkIfExists: true) : file("${workDir}/NO_HOTSPOTS")
    def blacklist_ch  = (params.blacklist && params.blacklist != "None") ? file(params.blacklist, checkIfExists: true) : file("${workDir}/NO_BLACKLIST")
    def multiqc_cfg   = params.multiqc_configfile ? file(params.multiqc_configfile, checkIfExists: true) : file("${workDir}/NO_MQCFG")
    def tg_adapters   = params.trimgalore_adapters ? file(params.trimgalore_adapters, checkIfExists: true) : file("${workDir}/NO_TGA")
    def tm_adapters   = params.trimmomatic_adapters ? file(params.trimmomatic_adapters, checkIfExists: true) : file("${workDir}/NO_TMA")
    def logo_ch       = file(params.logo, checkIfExists: true)

    // External scripts
    def check_design_script_ch              = file("${params.src}/check_design.py", checkIfExists: true)
    def ITR_id_bwamem_script_ch             = file("${params.src}/ITR_id_bwamem.pl", checkIfExists: true)
    def ssDNA_to_bigwigs_script_ch          = file("${params.src}/ssDNA_to_bigwigs_FASTER_LOMEM.pl", checkIfExists: true)
    def makeSSMultiQCReport_script_ch       = file("${params.src}/makeSSMultiQCReport_nextFlow.pl", checkIfExists: true)
    def pickNlines_script_ch                = file("${params.src}/pickNlines.pl", checkIfExists: true)
    def satCurveHS_script_ch                = file("${params.src}/satCurveHS.R", checkIfExists: true)
    def norm_script_ch                      = file("${params.src}/normalizeStrengthByAdjacentRegions.pl", checkIfExists: true)
    def reverse_script_ch                   = file("${params.src}/reverseStrandsForOriCalling.pl", checkIfExists: true)
    def getPeaksBedFiles_script_ch          = file("${params.src}/getPeaksBedFiles.pl", checkIfExists: true)
    def runSatCurve_script_ch               = file("${params.src}/runSatCurve.R", checkIfExists: true)
    def encode_idr_script_ch               = file("${params.src}/encode-dcc_chip-seq-pipeline2_src/encode_task_idr.py", checkIfExists: true)
    def get_frip_script_ch                  = file("${params.src}/compute_frip.py", checkIfExists: true)
    def plot_ssds_stat_script_ch            = file("${params.src}/plotSSDSqc.R", checkIfExists: true)
    def runPlotSSDSqc_script_ch             = file("${params.src}/runPlotSSDSqc.R", checkIfExists: true)

    // Genome files
    def bwa_ch      = file(params.custom_bwa, checkIfExists: true)
    def fasta_ch    = file(params.genome_fasta, checkIfExists: true)
    def index_amb   = file("${params.genomedir}/*.amb")
    def index_ann   = file("${params.genomedir}/*.ann")
    def index_bwt   = file("${params.genomedir}/*.bwt")
    def index_pac   = file("${params.genomedir}/*.pac")
    def index_sa    = file("${params.genomedir}/*.sa")
    def fai_ch      = file(params.fai, checkIfExists: true)

    //=========================================================================
    // SECTION 1 : INPUT SETTINGS
    //=========================================================================
    check_design(input_ch, check_design_script_ch)

    // Parse design CSVs into channels
    def fq_ch = check_design.out.reads_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }

    def ch_design_controls = check_design.out.controls_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(), row.multipleGroups.toBoolean() ] }

    makeScreenConfigFile()

    //=========================================================================
    // SECTION 2 : TRIMMING
    //=========================================================================
    def fq_tuple = fq_ch.map { it -> it[0,1].flatten() }

    crop(fq_tuple, makeScreenConfigFile.out.conf, makeScreenConfigFile.out.ok)

    trimming(crop.out.cropped, tg_adapters, tm_adapters)

    //=========================================================================
    // SECTION 3 : MAPPING AND PARSING
    //=========================================================================
    bwaAlign(
        trimming.out.trimmed,
        bwa_ch,
        fasta_ch,
        index_amb,
        index_ann,
        index_bwt,
        index_pac,
        index_sa,
        fai_ch
    )

    filterBam(bwaAlign.out.sorted_bam)

    // If with_control, filter out control samples from parseITRs input
    def parseITRs_input
    if (params.with_control) {
        parseITRs_input = ch_design_controls
            .combine(filterBam.out.filtered_bam)
            .map { it -> [ it[0], it[1], it[5], it[5].split('_')[0..-2].join('_'), it[6] ] }
            .filter { it -> it[0] == it[3] }
            .map { it -> it[2,4] }
    } else {
        parseITRs_input = filterBam.out.filtered_bam
    }

    parseITRs(parseITRs_input, fasta_ch, ITR_id_bwamem_script_ch)

    //=========================================================================
    // SECTION 4 : BIGWIG
    //=========================================================================
    makeBigwig(parseITRs.out.bed_to_bw, chrsize_ch)

    // Optional: makeBigwigReplicates
    if ( (params.bigwig_profile == "T1rep" || params.bigwig_profile == "T12rep") && params.nb_replicates == "2" ) {
        def bw_rep_ch = parseITRs.out.bed_to_bw
            .map { it -> [ it[0].split('_')[0..-3].join('_'), it[1], it[2], it[3] ] }
            .groupTuple(by: [0])
            .map { it -> it[0,1,2,3].flatten() }
        makeBigwigReplicates(bw_rep_ch)
    } else if (params.bigwig_profile == "T1rep" || params.bigwig_profile == "T12rep") {
        // DSL1 treated this combination as a fatal config error (previousmain.nf:1178-1181) rather
        // than silently skipping the process -- see DSL2_PORT_REVIEW.md point 2
        exit 1, "Error : selected bigwig_profile is ${params.bigwig_profile} but nb_replicates is ${params.nb_replicates}, check the parameters."
    }

    // Optional: kbrick bigwig
    if ( params.kbrick_bigwig ) {
        // bam_with_idx emits multiple bam types per sample; flatten for per-bam processing
        def bam_dt_ch = parseITRs.out.bam_with_idx.transpose()
        makeDeeptoolsBigWig(bam_dt_ch)

        def bam_fr_ch = parseITRs.out.bam_with_idx.transpose()
        toFRBigWig(bam_fr_ch, ssDNA_to_bigwigs_script_ch, fasta_ch, fai_ch)
    }

    //=========================================================================
    // SECTION 5 : PEAK CALLING
    //=========================================================================

    def satCurvePCs_ch = channel.from(satCurvePCsList)

    if (params.with_control) {
        // Group T1 beds by sample group (strip replicate suffix)
        def t1_grouped = parseITRs.out.t1_bed
            .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
            .groupTuple(by: [0])
            .map { it -> [ it[0], it[1].flatten() ] }

        def input_bed_grouped = filterBam.out.input_bed
            .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
            .groupTuple(by: [0])
            .map { it -> [ it[0], it[1].flatten() ] }

        // Join IP T1 beds with control beds via design -> tuple(id_ip, id_ct, ip_bed, ct_bed)
        def t1_with_ct = ch_design_controls
            .combine(t1_grouped)
            .combine(input_bed_grouped)
            .filter { it -> it[0] == it[5] && it[1] == it[7] }
            .map { it -> it[0,1,6,8].flatten() }

        shufBEDs_ct(t1_with_ct)

        // Fan out per saturation-curve percentage while keeping each IP group's OWN control
        // bed correctly paired throughout. shufBEDs_ct.out.shuf_beds already carries the right
        // [id_ip, ip_shuf_bed, ct_shuf_bed] triple per group (joined above); the previous code
        // threw that pairing away with `.first()`, broadcasting a single group's control to
        // every group in multi-group designs (silently wrong peak calls) -- see
        // DSL2_PORT_REVIEW.md point 1. multiMap keeps both derived channels emission-for-emission
        // in lockstep, so callPeaks always receives the matching control for each IP task.
        def peaks_ct_input = shufBEDs_ct.out.shuf_beds
            .combine(satCurvePCs_ch)
            .multiMap { id_ip, ip_bed, ct_bed, pct ->
                ip_tuple: tuple(id_ip, ip_bed, pct)
                ct_path:  ct_bed
            }

        callPeaks(
            peaks_ct_input.ip_tuple,
            peaks_ct_input.ct_path,
            blacklist_ch,
            fai_ch,
            pickNlines_script_ch,
            satCurveReps
        )
    } else {
        def t1_bed_grouped = parseITRs.out.t1_bed
            .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
            .groupTuple(by: [0])
            .map { it -> [ it[0], it[1].flatten() ] }

        shufBEDs(t1_bed_grouped)

        def shuf_with_pct = shufBEDs.out.ip_shuf.combine(satCurvePCs_ch)

        callPeaks(
            shuf_with_pct,
            file("${workDir}/NO_CONTROL"),
            blacklist_ch,
            fai_ch,
            pickNlines_script_ch,
            satCurveReps
        )
    }

    // Read genome_size from file (resume-safe replacement for stdout)
    def genome_size_ch = callPeaks.out.genome_size.first().map { it -> it.text.trim() }

    //=========================================================================
    // SECTION 6 : SSDS QC
    //=========================================================================
    // samStats on all parsed bam types
    def bam_ss_ch = parseITRs.out.bam_with_idx.transpose()
    samStats(bam_ss_ch)

    makeSSreport(parseITRs.out.itr_bed, hotspots_ch, makeSSMultiQCReport_script_ch)

    makeFingerPrint(
        filterBam.out.unparsed_bam.collect(),
        filterBam.out.unparsed_bai.collect()
    )

    // computeFRIP (only when hotspots == "None", i.e. de novo)
    def ssds_report_denovo_ch = channel.empty()
    def multiqc_denovo_ch     = channel.empty()
    def multiqc_ch            = channel.empty()

    if (params.hotspots == "None") {
        // Join each replicate's own bam-set with its own group's finalPeaks_noIDR.bed and its
        // own SSDSreport, instead of broadcasting the very first emission of each to every
        // sample via .first() -- same bug class as DSL2_PORT_REVIEW.md point 1 (and confirmed
        // against previousmain.nf:1762-1769, which does an explicit per-sample join here too).
        // parseITRs.out.itr_bam is keyed by sampleId (per replicate); callPeaks.out.final_peaks_noidr
        // is keyed by id_ip (per group, replicates sharing a peak call use the same one) -- same
        // '_'-strip used to build callPeaks' own input groups earlier in this workflow.
        def frip_input = parseITRs.out.itr_bam
            .join(makeSSreport.out.ssds_report, by: 0)
            .map { row -> tuple(row[0].split('_')[0..-2].join('_'), row) }
            .join(callPeaks.out.final_peaks_noidr, by: 0)
            .multiMap { id_ip, row, peaks ->
                bam_tuple: tuple(*row[0..-2])
                report:    row[-1]
                peaks_bed: peaks
            }

        computeFRIP(
            frip_input.bam_tuple,
            frip_input.peaks_bed,
            frip_input.report,
            get_frip_script_ch
        )

        if (params.with_ssds_multiqc) {
            ssds_multiqc_report(
                computeFRIP.out.ssds_report_denovo,
                runPlotSSDSqc_script_ch,
                plot_ssds_stat_script_ch
            )
            ssds_report_denovo_ch = computeFRIP.out.ssds_report_denovo_flat
            multiqc_denovo_ch     = ssds_multiqc_report.out.multiqc_plots
        }
    }

    if (params.with_ssds_multiqc && params.hotspots != "None") {
        ssds_multiqc_report(
            makeSSreport.out.ssds_report,
            runPlotSSDSqc_script_ch,
            plot_ssds_stat_script_ch
        )
        multiqc_ch = ssds_multiqc_report.out.multiqc_plots
    }


    //=========================================================================
    // SECTION 7 : OPTIONAL IDR ANALYSIS
    //=========================================================================
    def idr_final_peaks_ch = channel.empty()

    if (params.with_idr && params.nb_replicates == "2") {

        if (params.with_control) {
            // Prepare replicate channel with controls
            def t1_rep_ct = ch_design_controls
                .combine(parseITRs.out.t1_bed)
                .combine(filterBam.out.input_bed)
                .filter { it -> it[0] == it[5].split('_')[0..-2].join('_') && it[1] == it[7].split('_')[0..-2].join('_') }
                .map { it -> [it[0], it[1], it[6], it[8]] }

            renameInputBed(t1_rep_ct)

            def rep_renamed = renameInputBed.out.renamed
                .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1].split('_')[0..-2].join('_'), it[2], it[3] ] }
                .groupTuple(by: [0])
                .groupTuple(by: [1])
                .map { it -> [ it[0], it[1][0], it[2], it[3] ] }
                .map { it -> it[0,1,2,3].flatten() }

            createPseudoReplicates_ct(rep_renamed)

            // Same fix as callPeaks above: keep each IP group's own ct_pool.bed paired via
            // multiMap instead of broadcasting one group's control to all via .first()
            // (DSL2_PORT_REVIEW.md point 1). genome_size_ch.first() below is unrelated and
            // fine to broadcast: genome size is computed identically for every group (it only
            // depends on the shared fai/blacklist, not per-group data).
            def idr_ct_input = createPseudoReplicates_ct.out.all_rep
                .multiMap { id_ip, t1_rep1, t1_rep2, r1pr1, r1pr2, r2pr1, r2pr2, pool_r1, pool_r2, poolT, ct_pool ->
                    rep_tuple: tuple(id_ip, t1_rep1, t1_rep2, r1pr1, r1pr2, r2pr1, r2pr2, pool_r1, pool_r2, poolT)
                    ct_path:   ct_pool
                }

            callPeaksForIDR(
                idr_ct_input.rep_tuple,
                genome_size_ch,
                idr_ct_input.ct_path
            )
        } else {
            // Without control: group replicates
            def t1_rep = parseITRs.out.t1_bed
                .map { it -> [ it[0].split('_')[0..-3].join('_'), it[1] ] }
                .groupTuple(by: [0])
                .map { it -> it[0,1].flatten() }

            createPseudoReplicates(t1_rep)

            callPeaksForIDR(
                createPseudoReplicates.out.all_rep,
                genome_size_ch,
                file("${workDir}/NO_CONTROL")
            )
        }

        IDRanalysis(
            callPeaksForIDR.out.all_peaks_rep,
            blacklist_ch,
            chrsize_ch,
            encode_idr_script_ch,
            fai_ch
        )

        IDRpostprocess(IDRanalysis.out.idr_peaks)

        idr_final_peaks_ch = IDRpostprocess.out.idr_final_peaks
    }

    //=========================================================================
    // SECTION 8 : PEAK POST PROCESSING
    //=========================================================================
    if (!params.with_idr) {
        // normalizePeaks (no IDR)
        normalizePeaks(
            callPeaks.out.all_peaks,
            norm_script_ch,
            reverse_script_ch
        )

        // mergeFinalPeaks (only for 2 replicates)
        if (params.nb_replicates == "2") {
            def peaks_to_merge = callPeaks.out.all_peaks
                .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1], it[2] ] }
                .groupTuple(by: [0])
                .map { it -> it[0,2].flatten() }
            mergeFinalPeaks(peaks_to_merge)
        }
    }

    if (params.with_idr && params.nb_replicates == "2") {
        // normalizePeaks_idr
        normalizePeaks_idr(
            idr_final_peaks_ch,
            norm_script_ch,
            reverse_script_ch
        )
    }

    //=========================================================================
    // SECTION 8b : SATURATION CURVE
    //=========================================================================
    def satcurve_png_ch = channel.empty()
    if (params.satcurve) {
        makeSatCurve(
            callPeaks.out.allbed.collect(),
            getPeaksBedFiles_script_ch,
            runSatCurve_script_ch,
            satCurveHS_script_ch
        )
        satcurve_png_ch = makeSatCurve.out.png
    }

    //=========================================================================
    // SECTION 9 : GENERAL QC
    //=========================================================================
    // Collect all outputs for multiqc, using Channel.empty() for optional outputs
    general_multiqc(
        logo_ch,
        multiqc_cfg,
        makeFingerPrint.out.png.collect().ifEmpty([]),
        makeFingerPrint.out.tab.collect().ifEmpty([]),
        bwaAlign.out.bwa_flagstat.collect().ifEmpty([]),
        parseITRs.out.itr_flagstat.collect().ifEmpty([]),
        makeSSreport.out.ssds_report_flat.collect().ifEmpty([]),
        ssds_report_denovo_ch.collect().ifEmpty([]),
        multiqc_denovo_ch.collect().ifEmpty([]),
        multiqc_ch.collect().ifEmpty([]),
        samStats.out.samstats.collect().ifEmpty([]),
        crop.out.raw_fastqc.collect().ifEmpty([]),
        trimming.out.trim_fastqc.collect().ifEmpty([]),
        satcurve_png_ch.collect().ifEmpty([])
    )
}

//***************************************************************************//
//                          END OF PIPELINE                                  //
//***************************************************************************//

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Pipeline duration: $workflow.duration"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution profile: $workflow.profile"
    log.info "Script ID: $workflow.scriptId"
    log.info "Run name: $workflow.runName"
}

workflow.onError {
    log.error "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    log.error "Error report: ${workflow.errorReport}"
}