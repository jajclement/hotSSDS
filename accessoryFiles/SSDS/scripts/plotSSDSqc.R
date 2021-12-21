#!/usr/Rscript

# This script contains R functions for parsing and ploting data from tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# 2021 Pauline Auffret

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline) 
# and returns a tab containing only the data for mapping stats of different fragment types
get_totinfo <- function(tab,name) {
  names(tab) <- c("category","type","number")
  totinfo <- tab[tab$category=="totinfo",]
  tot<-totinfo[which(totinfo$type%in%c("ssDNA_type1_fragments",
                                       "ssDNA_type2_fragments",
                                       "dsDNA_fragments",
                                       "unclassified_fragments")),]
  tot[,4] <- rep(name,dim(tot)[1])
  names(tot) <- c("category","type","number","sample")
  return(tot)
}	

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# and returns a tab containing only the data for fragment length distribution 
get_frag <- function(tab,name) {
    names(tab) <- c("category","type","number")
    tab_frag <- tab[which(tab$category%in%c('Fragments_ssDNA_type1',
                                       'Fragments_ssDNA_type2',
                                       'Fragments_dsDNA',
                                       'Fragments_unclassified')),]
    tab_frag[,4] <- rep(name,dim(tab_frag)[1])
    names(tab_frag) <- c("category","type","number","sample")
    return(tab_frag)
}	

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# and returns a tab containing only the data for micro homology (uH) length distribution 
get_uH <- function(tab,name) {
  names(tab) <- c("category","type","number")
  tab_uH <- tab[which(tab$category%in%c('uH_ssDNA_type1',
                                   'uH_ssDNA_type2',
                                   'uH_dsDNA',
                                   'uH_unclassified')),]
  tab_uH[,4] <- rep(name,dim(tab_uH)[1])
  names(tab_uH) <- c("category","type","number","sample")
  return(tab_uH)
}	

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# and returns a tab containing only the data for offset length distribution
get_offset <- function(tab,name) {
  names(tab) <- c("category","type","number")
  tab_offset <- tab[which(tab$category%in%c('Offset_ssDNA_type1',
                                       'Offset_ssDNA_type2',
                                       'Offset_dsDNA',
                                       'Offset_unclassified')),]
  tab_offset[,4] <- rep(name,dim(tab_offset)[1])
  names(tab_offset) <- c("category","type","number","sample")
  return(tab_offset)
}	

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# and returns a tab containing only the data for ITR length distribution
get_ITR <- function(tab,name) {
  names(tab) <- c("category","type","number")
  tab_itr <- tab[which(tab$category%in%c('ITR_ssDNA_type1',
                                       'ITR_ssDNA_type2',
                                       'ITR_dsDNA',
                                       'ITR_unclassified')),]
  tab_itr[,4] <- rep(name,dim(tab_itr)[1])
  names(tab_itr) <- c("category","type","number","sample")
  return(tab_itr)
}

# This function will parse tab text file (output from ssds_multiqc process in ssdsnextflowpipeline)
# and returns a tab containing only the data for FRIP percentage
get_frip <- function(tab,name) {
  names(tab) <- c("category","type","number")
  tab_frip <- tab[which(tab$category%in%c('FRIP_ssType1',
                                         'FRIP_ssType2',
                                         'FRIP_dsDNA',
                                         'FRIP_unclassified')),]
  tab_frip[,4] <- rep(name,dim(tab_frip)[1])
  names(tab_frip) <- c("category","type","number","sample")
  return(tab_frip)
}		

# This function will plot a barplot for different types of mapping stats
#plot_barplot_totinfo <- function(tab) {
#  p <- ggplot(tab, aes(x = sample, y=number))+
#    geom_col(aes(fill = type), width = 0.7) +
#    xlab("Sample") +
#    ylab("Number of fragments") +
#    ggtitle("SSDS alignment stats") +
#    coord_flip() +
#    scale_y_continuous(labels = scientific)
#  return(p)
#}

# This function will plot a barplot for different types of mapping stats
plot_barplot_totinfo <- function(tab) {
  p <- ggplot(tab, aes(x = sample, y=number))+
    geom_col(aes(fill = type) , width = 0.7) +
    xlab("Samples") +
    ylab("Number of fragments") +
    labs(fill = "Type of fragment") +
    labs(title = "SSDS parsing statistics",
         subtitle = "Number of single stranded type 1 and type 2 fragments,\ndouble stranded fragments and unclassified fragments" ) +
    coord_flip() +
    theme(axis.text.x = element_text(colour = "grey20", size = 5),
          axis.text.y = element_text(colour = "grey20", size = 5),
          text=element_text(size = 8)
    ) +
    theme(
      plot.title = element_text(size = 8),    # Position et taille du titre au centre
      plot.subtitle = element_text(size = 6)
    )
  return(p)
}

# This function will plot a scatterplot for one type of fragment
plot_scatter <- function(tab, name) {
  tab$type<-as.numeric(tab$type)
  tab$number<-as.numeric(tab$number)
  p <- ggplot(tab, aes(x=type, y=number)) + 
    geom_point(aes(col=category),size=1.5) + 
    labs(y="Number of fragments", 
         x="Size of fragments (bp)", 
         title=name) +
  scale_y_continuous(labels = scientific)
  return(p)
}

#♣ This function will plot a barplot for FRIP scores
plot_barplot_frip <- function(tab) {
  p <- ggplot(tab, aes(x = category, y=number, fill=type))+
    geom_bar(stat="identity", position=position_dodge()) +
    xlab("Reference") +
    ylab("FRIP") +
    ggtitle("FRIP score stats") +
    scale_y_continuous(labels = percent) 
  return(p)
}
