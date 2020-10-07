#!/usr/Rscript
options(echo=FALSE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)

treatment <- args[1]
control   <- args[2]
ncis_path  <- args[3]
output    <- args[4]

print(paste0("Treatment BED :",treatment))
print(paste0("Control BED   :",control))
print(paste0("NCIS Path     :",ncis_path))
print(paste0("Output file   :",output))

##-----------------------------------------------------------------------------------
library("NCIS", lib.loc=ncis_path)
library('rtracklayer')
library("ShortRead")

res <- NCIS(treatment,control,"BED")
write(paste(res$est,
            res$pi0,
            res$binsize.est,
            res$r.seq.depth,
            sep = "\t"),
      output, 
      sep = "\t")
