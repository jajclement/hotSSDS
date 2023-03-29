#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
#######################################################################
## Use Deeptools to compute FRIP scores
## Created on 2021-08-10
#######################################################################
#######################################################################
import deeptools.countReadsPerBin as crpb
import os
import sys
import argparse
import pysam

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = "Compute FRIP score from bam file and peak file using deeptools."
Epilog = """Example usage: python compute_frip.py <PEAK_FILE> <BAM_FILE> <NCPU> <DNA_TYPE> <GENOME_NAME> <FRIP_FILE_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('PEAK_FILE', help="Input BED file containing peaks.")
argParser.add_argument('BAM_FILE', help="Input BAM file.")
argParser.add_argument('NCPU', help="Number of CPUs.")
argParser.add_argument('DNA_TYPE', help="ssType1 ; ssType2 ; dsDNA ; unclassified")
argParser.add_argument('GENOME_NAME', help="Genome name (e.g. mm10)")
argParser.add_argument('FRIP_FILE_OUT', help="Output file containing FRIP scores.")
args = argParser.parse_args()


############################################
############################################
## MAIN FUNCTION
############################################
############################################

def get_frip(peaks,bam,cpus,dna_type,genome_name,out):

    ERROR_STR = 'ERROR: Please check input files'
    
    # Open files
    o = open(out, 'w')
	
    # The FRiP score is defined as the fraction of reads that fall into a peak and is often used as a measure of ChIP-seq quality. 
    # For this example, we need a BED file containing the peak regions. Such files are usually computed using a peak caller. 
    cr = crpb.CountReadsPerBin([bam], bedFile=peaks, numberOfProcessors=int(cpus))
    reads_at_peaks = cr.run()
    #The result is a numpy array with a row for each peak region and a column for each BAM file.

    # Now, the total number of reads per peaks per bam file is computed:
    total = reads_at_peaks.sum(axis=0)

    # Next, we need to find the total number of mapped reads in bam file. For this we use the pysam module.
    nReads = pysam.AlignmentFile(bam)

    # Now, b.mapped contains the total number of mapped reads in bam file.
    # Finally, we can compute the FRiP score:
    frip = float(total[0]) / nReads.mapped
	
    # Write FRIP score in ouptut file
    o.write('\t'.join(['FRIP_'+str(dna_type),str(genome_name)+':de-novo',str(frip)]) + '\n')
	
    # Close files
    o.close()
	

############################################
############################################
## RUN FUNCTION
############################################
############################################

get_frip(peaks=args.PEAK_FILE,bam=args.BAM_FILE,cpus=args.NCPU,dna_type=args.DNA_TYPE,genome_name=args.GENOME_NAME,out=args.FRIP_FILE_OUT)
