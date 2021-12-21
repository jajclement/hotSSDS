#!/usr/Rscript

options(echo=FALSE) # if you want see commands in output file

# get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
functionFile <- args[1]
tabFile <- args[2]
sampleName <- args[3]
outputPath <- args[4]

# Source function files
source(functionFile)

# Load libraries
library(ggplot2)
library(scales)
library(gridExtra)

# Read tab file and store content in tab data frame
tab <- read.table(file=tabFile, sep="\t", h=F)

# Get fragment stats
tot <- get_totinfo(tab,sampleName)

# Get different types of fragment length distribution
frag <- get_frag(tab,sampleName)
uH <- get_uH(tab,sampleName)
itr <- get_ITR(tab,sampleName)
off <- get_offset(tab,sampleName)
frip <- get_frip(tab,sampleName)

# Plot stats and distribution
p0 <- plot_barplot_totinfo(tot)
p1 <- plot_barplot_frip(frip)
p2 <- plot_scatter(frag, "Fragment length distribution", "Fragments")
p3 <- plot_scatter(uH, "Micro homology length distribution", "Micro-homologies")
p4 <- plot_scatter(itr, "ITR length distribution", "ITR")
p5 <- plot_scatter(off, "Offset length distribution", "Offset")

# Open pdf file and print final plot
pdf(paste(outputPath,"/all_plots_",sampleName,"_mqc.pdf", sep=""))
grid.arrange(p0,p1,p2,p3,p4, ncol=2)
dev.off()

# Print plots individually in png files
ggsave(paste(outputPath,"/barplot_ssds_stats_",sampleName,"_mqc.png",sep=""),
	 plot = p0, scale = 1, width = 8, height = 4, units = "cm", dpi = 300)

ggsave(paste(outputPath,"/barplot_frip_score_",sampleName,"_mqc.png",sep=""),
          plot = p1, scale = 1, width = 8, height = 4, units = "cm", dpi = 300)

ggsave(paste(outputPath,"/scatter_frag_length__",sampleName,"_mqc.png",sep=""),
          plot = p2, scale = 1, width = 7, height = 5, units = "cm", dpi = 300)

ggsave(paste(outputPath,"/scatter_microh_length_",sampleName,"_mqc.png",sep=""),
          plot = p3, scale = 1, width = 7, height = 5, units = "cm", dpi = 300)

ggsave(paste(outputPath,"/scatter_itr_length_",sampleName,"_mqc.png",sep=""),
          plot = p4, scale = 1, width = 7, height = 5, units = "cm", dpi = 300)

ggsave(paste(outputPath,"/scatter_offset_length_",sampleName,"_mqc.png",sep=""),
          plot = p5, scale = 1, width = 7, height = 5, units = "cm", dpi = 300)

