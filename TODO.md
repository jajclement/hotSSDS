# TODO list for ssds nextflow pipeline
## Little tasks
* Rename pipeline... or not ? ssds_chipseq
* Check if output files in processes match publishDir files
* **Comment code in main.nf** started 2020-12-01
* Connect the filter process or not ?

## Medium tasks
* Implement log.info with parameters or reorder parameters in report
* **Make a changelog file** started 2020-10-29
* Update README with examples & output description
* Add ``when`` tag for peak calling
* Test with nexflow latest version
* Consistency of variables calling in ``main.nf``

## Big tasks
* Implement a Global and pretty QC
* Test pipeline on SRA and Julie's data
* **Add (optional) IDR analysis process** started 2020-12-17 
* Migration to Singularity/Docker container
* Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation
* Fix the 'don't go beyond' resources function in config
* Test running pipeline on Genotoul cluster
* Pipeline starts from bam or another check point

## Bug to investigate
* bwa job cancelled time limit but only 48h
* Pipeline do not resume properly when job cancelled by SLURM
* satcurve process : one per bam or one per pipeline ?

## Features to consider
* Downstream analyses integration
* Handling n>2 replicates for IDR
* SRA inputs

## Done
* **Processing of control files -> not DSBED for peak calling but ouptut from bwa** ok 2021-01-08
* **Check if --version works in the main.nf** ok 2020-10-27
* **Remove build from conda environment yml files** ok 2020-10-27
* **Edit manifest section in nextflow.config** ok 2020-10-28
* **Genomes paths modified on cluster ?** start 2020-10-29 end 2020-10-30
* **filename trimming report** 
* **Correct bug with hard trimming**
* **Correct shuffle_percent input un callPeaks process (it's satCurvePCs channel)** ok 2020-10-28
* **Clean macs2 script in callPeaks process** ish 2020-11-05
* **Reorder multiQC output report** ok 2020-10-28
* **No error when pipeline is launched with option "with-control" without control** ok 2020-11-04
* **No error when empty bam** ok 2020-11-06 
* **Add multimapper handling process** ok 2020-11-04
* **git merge branch dev when peak calling is functionnal** ok 2020-12-09
* **Update  -h option in main.nf** ok 2020-12-02
* **Comment code in main.nf** ok 2020-12-04
* *Implement frombam option* no longer considered
* *Handling of concatenation of SRA files* no longer considered



