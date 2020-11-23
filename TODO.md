# TODO list for ssds nextflow pipeline
## Little tasks
* Rename pipeline... or not ? ssds_chipseq
* git merge branch dev when peak calling is functionnal
* Update  -h option in main.nf
* Check if output files in processes match publishDir files
* Comment code in main.nf

## Medium tasks
* Implement log.info with parameters or reorder parameters in report
* Make a changelog file start 2020-10-29
* Update README with examples & output description
* Handling of concatenation of SRA files

## Big tasks
* Implement a Global and pretty QC
* Test pipeline on SRA and Julie's data
* Add (optional) IDR analysis process
* Migration to Singularity/Docker container
* Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation
* Implement frombam option
* Fix the 'don't go beyond' resources function in config
* Test running pipeline on Genotoul cluster

## One day maybe
* Replicates handling
* SRA inputs
* controls when not in input (dsDNA ? pooled ssDNA ?)

## Done
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

