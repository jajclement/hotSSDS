module load nextflow/0.30.2

nextflow run -c $NXF_PIPEDIR/nextflow.local.config $NXF_PIPEDIR/SSDSPipeline_1.6.groovy \
    --fq1 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq \
    --fq2 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq \
    --r1Len 36 \
    --r2Len 40 \
    --genome mm10 \
    --outName testSSDS1.6 \
    --outdir SSDS1.6_test \
    -with-trace -with-timeline
