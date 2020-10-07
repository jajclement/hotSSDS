#!/usr/bin/perl

my $n = `wc -l $ARGV[0]`;chomp $n;

my $pc = $ARGV[1]/$n;

open IN, $ARGV[0];

while (<IN>){
    chomp;
    print $_."\n" if (rand() <= $pc);
}
