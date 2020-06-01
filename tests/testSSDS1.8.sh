module load nextflow/0.30.2

nextflow run -c $NXF_PIPEDIR/config.nf -profile local \
    $NXF_PIPEDIR/SSDSPipeline_1.8.nf \
    --fq1 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq \
    --fq2 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq \
    --r1Len 36 \
    --r2Len 40 \
    --genome mm10 \
    --name testSSDS1.8 \
    --outdir SSDS1.8_test 
