# TODO list for ssds nextflow pipeline
## Little tasks
* Rename pipeline... or not ?
* git merge branch dev when peak calling is functionnal
* Update  -h option in main.nf
* **Check if --version works in the main.nf** ok 2020-10-27
* **Remove build from conda environment yml files** ok 2020-10-27
* Edit manifest section in nextflow.config

## Medium tasks
* Correct shuffle_percent input un callPeaks process (it's satCurvePCs channel)
* Implement log.info with parameters or reorder parameters in report
* Clean macs2 script in callPeaks process
* Make a changelog file
* Reorder multiQC output report
* Update README with examples & output description

## Big tasks
* Implement a Global and pretty QC
* Test pipeline on SRA and Julie's data
* Add (optional) IDR analysis process
* Migration to Singularity/Docker container
* Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation

## Things to think about
* Replicates handling
* SRA inputs
* controls when not in input (dsDNA ? pooled ssDNA ?)



