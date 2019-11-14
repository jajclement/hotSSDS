#!/usr/bin/perl
use strict;

my ($samIn,$samOut) = @ARGV;

my $tmp = '/lscratch/'.$ENV{SLURM_JOBID}.'/'.(rand()*1000000000).'.txt';

open OUT, '>', $tmp;

## Count instances of each read (Supplementary alignments will result in >2 entries per pair
my $cmd = 'cut -f1 '.$samIn.' |uniq -c ';
open my $PIPE, '-|', $cmd;

while(<$PIPE>){
	chomp; 
	$_ =~ s/^\s+(\d+)\s+(\S+)/$1\t$2/; 
	my ($n,$nm) = ($1,$2);
	
	my $v = (($nm =~ /^\@/)?"":"ss:Z:$n");

	for my $i(1..$n){
		print OUT $v."\n"
	}
}

close OUT; close $PIPE;

## ADD The ss:Z flag to SAM and un-check flag 2 where necessary
open OUTsam,  '>', $samOut;

my $cmd2 = 'paste '.$samIn.' '.$tmp;
open my $PIPE2, '-|', $cmd2;

while(<$PIPE2>){
	chomp;
	my @F = split(/\t/,$_);

	## SET THE "MAPPED IN PROPER PAIR" FLAG OFF
	if ($_ !~ /ss:Z:2\s*$/ && $_ !~ /^\@/){
		if (($F[1] & 2) == 2){
			$F[1] = $F[1] ^ 2;
		}
	}
	
	print OUTsam join("\t",@F)."\n";
}

close OUTsam; 
