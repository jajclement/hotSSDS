// --------------------------------------------------------
//  SSDSnextflowPipeline shenron cluster getFQ config file 
//  --------------------------------------------------------

//NEXTFLOW REPORT PARAMETERS
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


//MAIN PROCESS DEFINITION  
process {
 project = 'SSDS_alignment_Pipeline'
	
 //MODULE DETAILS
 $getSRAfiles.module = ['sratoolkit/2.9.2','picard/2.9.2','samtools/1.8']
 $getFromObjectStore.module = ['sratoolkit/2.9.2','picard/2.9.2','samtools/1.8']
 $bamToFastq.module = ['picard/2.9.2','samtools/1.8']
 $initFQtoFQ_PE.module = ['picard/2.9.2','samtools/1.8']
 $initFQtoFQ_SR.module = ['picard/2.9.2','samtools/1.8']
 $fastQCall.module = ['fastqtools/0.8','fastqc/0.11.8','bwa/0.7.17']
 $mergeFQ.module = ['fastqtools/0.8']

 //DEFAULT CONFIGURATION
 cpus = { 1 * task.attempt }
 memory = { 8.GB * task.attempt } 
 time = { 2.h * task.attempt }
 executor = 'slurm'
 scratch = '${SCRATCH}'
 errorStrategy = { task.exitStatus == 143 ? 'retry' : 'finish' }
 maxRetries = 2
 maxErrors = '-1'
 env{
    TMPDIR='${SCRATCH}'
 }

 //MAX CONFIG
 max_memory = { 128.GB }
 max_cpus = { 32 }
 max_time = { 96h } 

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
  }
}
