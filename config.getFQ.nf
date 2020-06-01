
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
    process.scratch = '/lscratch/$SLURM_JOBID'
    process.clusterOptions = ' --gres=lscratch:800 '
    env{
      TMPDIR='/lscratch/$SLURM_JOBID'
    }
  }

  none {
    // Add custom configs here
  }
}

params.outdir = './nxfOut'

report {
  enabled = false
}

timeline {
  enabled = false
}

trace {
  enabled = false
}

manifest {
  description = '2020: Kevin Brick'
}

process {
	//MODULE DETAILS
	$getSRAfiles.module = ['sratoolkit/2.9.2','picard/2.9.2','samtools/1.8']
	$getFromObjectStore.module = ['sratoolkit/2.9.2','picard/2.9.2','samtools/1.8']
	$bamToFastq.module = ['picard/2.9.2','samtools/1.8']
	$initFQtoFQ_PE.module = ['picard/2.9.2','samtools/1.8']
  $initFQtoFQ_SR.module = ['picard/2.9.2','samtools/1.8']
  $fastQCall.module = ['fastqtools/0.8','fastqc/0.11.8','bwa/0.7.17']
  $mergeFQ.module = ['fastqtools/0.8']

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
	$getSRAfiles {
		cpus = { 2 }
		memory = { 4.GB }
		time = { 6.hour * task.attempt }
	  }

	$getFromObjectStore {
		cpus = { 4 }
		memory = { 8.GB }
		time = { 6.hour * task.attempt }
	  }

	$bamToFastq {
		cpus = { 4 }
		memory = { 16.GB * task.attempt }
		time = { 6.hour * task.attempt }
	  }

	$initFQtoFQ_PE {
		cpus = { 8 }
		memory = { 32.GB }
		time = { 6.hour * task.attempt }
	  }

  $initFQtoFQ_SR {
		cpus = { 8 }
		memory = { 32.GB }
		time = { 6.hour * task.attempt }
	  }

  $fastQCall {
		cpus = { 12 }
		memory = { 32.GB }
		time = { 8.hour * task.attempt }
	  }

  $mergeFQ {
		cpus = { 8 }
		memory = { 16.GB }
		time = { 12.hour * task.attempt }
	  }

  $omitFQ {
		cpus = { 2 }
		memory = { 2.GB }
		time = { 0.2.hour * task.attempt }
		clusterOptions = ' --gres=lscratch:60'
	  }

}
