#!/usr/bin/perl 
use strict;
use Getopt::Long;

#Perl script adapted from callSSDSpeaks pipeline (https://github.com/kevbrick/callSSDSpeaks) by Kevin Brick
#P. Auffret 2020 

my $tf;
my $dir;

GetOptions ('tf=s'	=> \$tf,
	    'dir=s'	=> \$dir);

open TMP, '>', $tf or die $!;

#my $cmd="cd $dir; wc -l *peaks_sc.bed | grep -v total | sort -k1n,1n";
my $cmd="cd $dir; wc -l *N*pc*peaks_sc.bed | grep -v total | sort -k1n,1n"; #add more precise regex 

print TMP join("\t","reads","pc","hs")."\n"; 	
open my $IN, '-|', $cmd;
while (<$IN>){
    chomp;
    next if ($_ =~ /\stotal\s*$/);	
    $_ =~ /^\s*(\d+).+\.N(\d+)_([\d\.]+)pc.+$/;
    my ($HS,$N,$pc) = ($1,$2,$3*100);
    print TMP join("\t",$N,$pc,$HS)."\n";
}

close TMP;
close $IN;
