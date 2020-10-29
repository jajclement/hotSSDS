# TODO list for ssds nextflow pipeline
## Little tasks
* Rename pipeline... or not ?
* git merge branch dev when peak calling is functionnal
* Update  -h option in main.nf
* **Check if --version works in the main.nf** ok 2020-10-27
* **Remove build from conda environment yml files** ok 2020-10-27
* **Edit manifest section in nextflow.config** ok 2020-10-28
* Check if output files in processes match publishDir files
* Genomes paths modified on cluster ?
* Comment code in main.nf

## Medium tasks
* **Correct shuffle_percent input un callPeaks process (it's satCurvePCs channel)** ok 2020-10-28
* Implement log.info with parameters or reorder parameters in report
* Clean macs2 script in callPeaks process
* Make a changelog file
* **Reorder multiQC output report ** ok 2020-10-28
* Update README with examples & output description
* No error when pipeline is launched with option "with-control" without control
* Handling of concatenation of SRA files
* Add multimapper handling process

## Big tasks
* Implement a Global and pretty QC
* Test pipeline on SRA and Julie's data
* Add (optional) IDR analysis process
* Migration to Singularity/Docker container
* Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation
* Implement frombam option
* Fix the 'don't go beyond' resources function in config
* Test running pipeline on Genotoul cluster

## Things to think about
* Replicates handling
* SRA inputs
* controls when not in input (dsDNA ? pooled ssDNA ?)



