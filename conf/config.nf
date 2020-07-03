
profiles {
  local {
    process.maxForks = 1
    process.executor='local'
    env{
      TMPDIR='./'
    }
  }

  slurm {
    process.executor='slurm'
//    process.scratch = '/lscratch/$SLURM_JOBID'
    process.scratch = '$SCRATCH/$SLURM_JOBID'
//    process.clusterOptions = ' --gres=lscratch:800 '
    env{
//      TMPDIR='/lscratch/$SLURM_JOBID'
        TMPDIR='$SCRATCH/$SLURM_JOBID'
    }
  }

  standard {
    process.executor='slurm'
    process.scratch = '/lscratch/$SLURM_JOBID'
    process.clusterOptions = ' --gres=lscratch:800 '
    env{
      TMPDIR='/lscratch/$SLURM_JOBID'
    }
  }
}

params.outdir = './nxfOut'

report {
  enabled = true
  file = "${params.outdir}/nxfReports/report.html"
}

timeline {
  enabled = true
  file = "${params.outdir}/nxfReports/timeline.html"
}

trace {
  enabled = true
  file = "${params.outdir}/nxfReports/trace.txt"
}

manifest {
  description = '2020: Kevin Brick'
}

process {
	//MODULE DETAILS
	$getFQs.module = ['nextflow/0.30.2']
	$runFASTQC.module = ['fastqtools/0.8','fastqc/0.11.8','bwa/0.7.17']
	$trimFASTQs.module = ['trimgalore/0.4.5','cutadapt/2.9']
	$bwaAlign.module = ['java/1.8.0_92','picard/2.9.2','fastxtoolkit/0.0.14','samtools/1.8']
  $mergeInitBAMs.module = ['java/1.8.0_92','picard/2.9.2','samtools/1.8']
  $parseITRs.module = ['picard/2.9.2','samtools/1.8']
  $gatherOutputs.module = ['picard/2.9.2','samtools/1.8']
  $makeDeeptoolsBigWig.module = ['samtools/1.8','deeptools/3.0.1','ucsc/365']
	$samStats.module = ['samtools/1.8']
	$toFRBigWig.module = ['samtools/1.8','ucsc/365','bedtools/2.25.0','picard/2.9.2']
	$makeSSreport.module = ['python/2.7','bedtools/2.25.0','samtools/1.8']
	$multiQC.module = ['python/2.7']

	project = 'SSDS_alignment_Pipeline'

	// Defaults maxs
	max_memory = 128.GB
	max_cpus = 16
	max_time = 240.h

	//DEFAULT PROCESS PROPS
	cpus = { 1 * task.attempt }
	memory = { 8.GB * task.attempt }
	time = { 2.h * task.attempt }

	errorStrategy = { task.exitStatus == 143 ? 'retry' : 'finish' }
	maxRetries = 2
	maxErrors = '-1'

	echo = true

	//PROCESS-SPECIFIC RESOURCES
	$getFQs {
		cpus = { 12 }
		memory = { 32.GB }
		time = { 24.hour * task.attempt }
		clusterOptions = ' --gres=lscratch:800'
	  }

	$runFASTQC {
		cpus = { 16 }
		memory = { 16.GB }
		time = { 8.hour * task.attempt }
	  }

	$trimFASTQs {
		cpus = { 16 }
		memory = { 4.GB * task.attempt }
		time = { 8.hour * task.attempt }
	  }

	$bwaAlign {
		cpus = { 16 }
		memory = { 32.GB }
		time = { 60.hour * task.attempt }
	  }

  $mergeInitBAMs {
		cpus = { 8 }
		memory = { 32.GB }
		time = { 6.hour * task.attempt }
	  }

  $parseITRs {
		cpus = { 16 }
		memory = { 32.GB }
		time = { 8.hour * task.attempt }
	  }

  $gatherOutputs {
		cpus = { 16 }
		memory = { 32.GB }
		time = { 8.hour * task.attempt }
	  }

  $makeDeeptoolsBigWig {
		cpus = { 16 }
		memory = { 32.GB }
		time = { 8.hour * task.attempt }
	  }

	$samStats {
		cpus = { 2 }
		memory = { 16.GB }
		time = { 6.hour * task.attempt }
	  }

	$toFRBigWig {
		cpus = { 2 }
		memory = { 8.GB * task.attempt }
		time = { 6.hour * task.attempt }
	  }

	$makeSSreport {
		cpus = { 2 }
		memory = { 8.GB * task.attempt }
		time = { 6.hour * task.attempt }
	  }

	$multiQC {
		cpus = { 2 }
		memory = { 8.GB * task.attempt }
		time = { 2.hour * task.attempt }
	  }
}
