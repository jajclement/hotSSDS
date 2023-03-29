#!/usr/bin/perl
use strict;

while (<>){
	$_ =~ s/\+/REPLACEPOS/;
	$_ =~ s/\-/\+/;
	$_ =~ s/REPLACEPOS/\-/;
	print $_;
}

