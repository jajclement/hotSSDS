#!/usr/Rscript 
options(echo=FALSE) # if you want see commands in output file 

args <- commandArgs(trailingOnly = TRUE)

sourceFile <- args[1]
satTable <- args[2]
sampleName <- args[3]

source(sourceFile)

satCurveHS(fIN = satTable, sampleName = sampleName)
