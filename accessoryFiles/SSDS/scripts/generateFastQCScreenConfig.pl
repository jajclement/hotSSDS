#!/usr/bin/perl

use strict;

my $genomes2use  = $ARGV[0];
my $nThreads     = $ARGV[1]?$ARGV[1]:1;
#my $genomeFolder = $ENV{NXF_GENOMES};
my $genomeFolder = $ARGV[2];

opendir(DIR,$genomeFolder);

print "# This is a config file for FastQ Screen\n\n";
print "THREADS $nThreads\n\n";

while (my $entry = readdir(DIR)){
  my $fa = $genomeFolder.'/'.$entry.'/BWAIndex/version0.7.10/genome.fa';
  if (-e $fa){
    if ($genomes2use =~ /($entry)/i){
      #print "$entry --> $fa\n" ;
      print "DATABASE $entry $fa\n"
    }
  }
}

#
#
# # This is an example configuration file for FastQ Screen
#
# ############################
# ## Bowtie, Bowtie 2 or BWA #
# ############################
# ## If the Bowtie, Bowtie 2 or BWA binary is not in your PATH, you can set
# ## this value to tell the program where to find your chosen aligner.  Uncomment
# ## the relevant line below and set the appropriate location.  Please note,
# ## this path should INCLUDE the executable filename.
#
# #BOWTIE	/usr/local/bin/bowtie/bowtie
# #BOWTIE2 /usr/local/bowtie2/bowtie2
# #BWA /usr/local/bwa/bwa
#
#
#
# ############################################
# ## Bismark (for bisulfite sequencing only) #
# ############################################
# ## If the Bismark binary is not in your PATH then you can set this value to
# ## tell the program where to find it.  Uncomment the line below and set the
# ## appropriate location. Please note, this path should INCLUDE the executable
# ## filename.
#
# #BISMARK	/usr/local/bin/bismark/bismark
#
#
#
# ############
# ## Threads #
# ############
# ## Genome aligners can be made to run across multiple CPU cores to speed up
# ## searches.  Set this value to the number of cores you want for mapping reads.
#
# THREADS		4
#
#
#
# ##############
# ## DATABASES #
# ##############
# ## This section enables you to configure multiple genomes databases (aligner index
# ## files) to search against in your screen.  For each genome you need to provide a
# ## database name (which can't contain spaces) and the location of the aligner index
# ## files.
# ##
# ## The path to the index files SHOULD INCLUDE THE BASENAME of the index, e.g:
# ## /data/RDCO/genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
# ## Thus, the index files (Homo_sapiens.GRCh37.1.bt2, Homo_sapiens.GRCh37.2.bt2, etc.)
# ## are found in a folder named 'GRCh37'.
# ##
# ## If, for example, the Bowtie, Bowtie2 and BWA indices of a given genome reside in
# ## the SAME FOLDER, a SINLGE path may be provided to ALL the of indices.  The index
# ## used will be the one compatible with the chosen aligner (as specified using the
# ## --aligner flag).
# ##
# ## The entries shown below are only suggested examples, you can add as many DATABASE
# ## sections as required, and you can comment out or remove as many of the existing
# ## entries as desired.  We suggest including genomes and sequences that may be sources
# ## of contamination either because they where run on your sequencer previously, or may
# ## have contaminated your sample during the library preparation step.
# ##
# ## Human - sequences available from
# ## ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/
# #DATABASE	Human	/data/RDCO/genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
# ##
# ## Mouse - sequence available from
# ## ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/
# #DATABASE	Mouse	/data/RDCO/genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37
# ##
# ## Ecoli- sequence available from EMBL accession U00096.2
# #DATABASE	Ecoli	/data/RDCO/genomes/Ecoli/Ecoli
# ##
# ## PhiX - sequence available from Refseq accession NC_001422.1
# #DATABASE	PhiX	/data/RDCO/genomes/PhiX/phi_plus_SNPs
# ##
# ## Adapters - sequence derived from the FastQC contaminats file found at: www.bioinformatics.babraham.ac.uk/projects/fastqc
# #DATABASE	Adapters	/data/RDCO/genomes/Contaminants/Contaminants
# ##
# ## Vector - Sequence taken from the UniVec database
# ## http://www.ncbi.nlm.nih.gov/VecScreen/UniVec.html
# #DATABASE	Vectors		/data/RDCO/genomes/Vectors/Vectors
#
# DATABASE	Human	/data/RDCO/genomes/hg19/BWAIndex/version0.7.10/genome.fa
# DATABASE	Mouse	/data/RDCO/genomes/mm10/BWAIndex/version0.7.10/genome.fa
# DATABASE	Rat	/data/RDCO/genomes/rn5/BWAIndex/version0.7.10/genome.fa
# DATABASE	Yeast	/data/RDCO/genomes/sacCer3/BWAIndex/version0.7.10/genome.fa
# DATABASE	Ecoli	/data/RDCO/genomes/ecoli/BWAIndex/version0.7.10/genome.fa
# DATABASE	Bsub	/data/RDCO/genomes/bsub/BWAIndex/version0.7.10/genome.fa
# DATABASE	UniVec	/data/RDCO/genomes/UniVec_2015_noIllumina/BWAIndex/version0.7.10/genome.fa
# DATABASE	PhiX	/data/RDCO/genomes/phiX/BWAIndex/version0.7.10/genome.fa
# DATABASE	Adapter	/data/RDCO/genomes/illuminaAdapters/BWAIndex/version0.7.10/genome.fa
# DATABASE	Dog	/data/RDCO/genomes/canFam3/BWAIndex/version0.7.10/genome.fa
# DATABASE	Opossum	/data/RDCO/genomes/monDom/BWAIndex/version0.7.10/genome.fa
#
# #generateFastQCScreenConfig.pl
