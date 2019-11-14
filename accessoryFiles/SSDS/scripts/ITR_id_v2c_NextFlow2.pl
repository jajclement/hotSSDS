use strict;
use Bio::DB::Sam;
use Getopt::Long;

my ($in_bam,$g) = @ARGV;

die ("Genomes path NOT defined ... Please set environment variable \$GENOMES ...") unless ($ENV{GENOMES});

my $outname = $in_bam;
my $in_fa   = $ENV{GENOMES}.'/'.$g.'/BWAIndex/version0.7.10/genome.fa';

die ("Invalid genomes path ['.$ENV{GENOMES}.'] ...\n$in_fa does not exist ...") unless (-e ($in_fa));

## SET UP TEMP FILES FOR EACH SPLIT BAM
my $splitBEDss1 	= $outname.'.ssDNA_type1.bed';
my $splitBEDss2		= $outname.'.ssDNA_type2.bed';
my $splitBEDds		= $outname.'.dsDNA.bed';
my $splitBEDdsS		= $outname.'.dsDNA_strict.bed';
my $splitBEDunC		= $outname.'.unclassified.bed';

my $splitBEDSORTss1	= $outname.'.ssDNA_type1.sorted.bed';
my $splitBEDSORTss2	= $outname.'.ssDNA_type2.sorted.bed';
my $splitBEDSORTds	= $outname.'.dsDNA.sorted.bed';
my $splitBEDSORTdsS	= $outname.'.dsDNA_strict.sorted.bed';
my $splitBEDSORTunC	= $outname.'.unclassified.sorted.bed';

my $splitSAMss1		= $outname.'.ssDNA_type1.sam';
my $splitSAMss2		= $outname.'.ssDNA_type2.sam';
my $splitSAMds		= $outname.'.dsDNA.sam';
my $splitSAMdsS		= $outname.'.dsDNA_strict.sam';
my $splitSAMunC		= $outname.'.unclassified.sam';

my $splitBAMss1		= $outname.'.ssDNA_type1.bam';
my $splitBAMss2		= $outname.'.ssDNA_type2.bam';
my $splitBAMds		= $outname.'.dsDNA.bam';
my $splitBAMdsS		= $outname.'.dsDNA_strict.bam';
my $splitBAMunC		= $outname.'.unclassified.bam';

open (OUT_BEDss1, ">", $splitBEDss1) or die  " unable to open output file!\n";
open (OUT_BEDss2, ">", $splitBEDss2) or die  " unable to open output file!\n";
open (OUT_BEDds,  ">", $splitBEDds)  or die  " unable to open output file!\n";
open (OUT_BEDdsS, ">", $splitBEDdsS) or die  " unable to open output file!\n";
open (OUT_BEDunC, ">", $splitBEDunC) or die  " unable to open output file!\n";
open (OUT_SAMss1, ">", $splitSAMss1) or die  " unable to open output file!\n";
open (OUT_SAMss2, ">", $splitSAMss2) or die  " unable to open output file!\n";
open (OUT_SAMds,  ">", $splitSAMds)  or die  " unable to open output file!\n";
open (OUT_SAMdsS, ">", $splitSAMdsS) or die  " unable to open output file!\n";
open (OUT_SAMunC, ">", $splitSAMunC) or die  " unable to open output file!\n";

my $nMM = 1;

my $sam = Bio::DB::Sam->new(
		-bam  => $in_bam, 
		-fasta=>$in_fa,
		-autoindex  => 1,
		-expand_flags => 1,
		) or die " unable to open input bam file '$in_bam'!\n";


## NOT NEEDED: KB 02-04-16
#	if ($output =~ /(sam|all|bedbam)/){
#		my $header = $sam->header();
#		die " no header in input bam file!!!\n" unless ($header);	
#		print OUT_SAM $header->text;
#	}

my $i = 0;
for my $tid (0 .. $sam->n_targets - 1)  {
	my $seq_id = $sam->target_name($tid);
	my $iter_pair = $sam->features(
			'-type'     => 'read_pair',
			'-iterator' => 1,
			'-seq_id'   => $seq_id);

	while (my $pair = $iter_pair->next_seq() ) {
	
		my ($left,$right) = $pair->get_SeqFeatures;
	

		last if $i > 1E9 ;
		next if !defined($left) || !defined($right)|| $left->unmapped || $right->unmapped;
		next if (!$left->proper_pair || $left->seq_id ne $right->seq_id || ($right->end - $left->start >10000));
		
		my $first_mate; my $second_mate;
	
		if ($left->get_tag_values('FIRST_MATE')) {
			$first_mate = $left; $second_mate = $right; 
		} else {
			$first_mate = $right; $second_mate = $left;
		}
	
		my ($ref_2,$matches_2,$query_2) = $second_mate->padded_alignment;
		my ($ref_1,$matches_1,$query_1) = $first_mate->padded_alignment;
	
		my ($ITR_len, $aln_offset, $uH_len, $aln_correction) = find_offset($ref_2,$query_2,$first_mate->target->dna,$second_mate->strand,1);
	
		# adjust coord/adjust offset in seq_coord; # deal with soft-clipping? - incorporate "super-strict" filtering
		# 5' of read 2 must be shifted for fill-in part
	
		# generate output;
		my ($strand,$start,$end,$chrom,$f_qual,$s_qual);
	
		$strand = ($first_mate->strand == 1)? "+": "-";
		$start  = ($first_mate->strand == 1)? $left->start - 1 : $left->start + $aln_correction - 1;
		$end    = ($first_mate->strand == 1)? $right->end - $aln_correction: $right->end;
		$chrom  = $first_mate->seq_id;
		$f_qual = $first_mate->qual;
		$s_qual = $second_mate->qual;
	
		next if ($left->get_tag_values('QC_FAILED') || $right->get_tag_values('QC_FAILED')); ## KB
	
		## Kevin 11/17/11 - output to SAM
		my @s_tam = split(/\t/,$second_mate->tam_line());
		my @f_tam = split(/\t/,$first_mate->tam_line());
	
		$f_tam[7] = $second_mate->start;
		$s_tam[7] = $first_mate->start;
				
		my ($s_line,$nI,$nD,$nS) = generateSAMline(2,\@s_tam,$ITR_len-$uH_len,$second_mate->strand);
		my ($f_line,$nA,$nD,$nS) = generateSAMline(1,\@f_tam,$ITR_len-$uH_len,$first_mate->strand,$nI,$nD,$nS);
	
		if ($f_line && $s_line){
			my $lineOut; 

			if ($ITR_len > 5 && $aln_offset > 2 && $uH_len >= 0){
				print OUT_BEDss1  join("\t", ($chrom,$start,$end,$f_qual."_".$s_qual,$ITR_len."_".$uH_len."_".$aln_offset,$strand)), "\n" ;
				print OUT_SAMss1 $f_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				print OUT_SAMss1 $s_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				$lineOut++;
			}
			
			if ($ITR_len > 5 && $aln_offset < 3 && $uH_len >= 0){
				print OUT_BEDss2  join("\t", ($chrom,$start,$end,$f_qual."_".$s_qual,$ITR_len."_".$uH_len."_".$aln_offset,$strand)), "\n" ;
				print OUT_SAMss2 $f_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				print OUT_SAMss2 $s_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";   
				$lineOut++;
			}

			if ($ITR_len < 3){
				print OUT_BEDds   join("\t", ($chrom,$start,$end,$f_qual."_".$s_qual,$ITR_len."_".$uH_len."_".$aln_offset,$strand)), "\n" ;
				print OUT_SAMds  $f_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				print OUT_SAMds  $s_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				$lineOut++;
			}

			if ($ITR_len < 1 && $uH_len < 1){
				print OUT_BEDdsS   join("\t", ($chrom,$start,$end,$f_qual."_".$s_qual,$ITR_len."_".$uH_len."_".$aln_offset,$strand)), "\n" ;
				print OUT_SAMdsS  $f_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				print OUT_SAMdsS  $s_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				$lineOut++;
			}

			unless ($lineOut){
				print OUT_BEDunC  join("\t", ($chrom,$start,$end,$f_qual."_".$s_qual,$ITR_len."_".$uH_len."_".$aln_offset,$strand)), "\n" ;
				print OUT_SAMunC  $f_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				print OUT_SAMunC  $s_line."\tit:i:$ITR_len\tuh:i:$uH_len\tos:i:$aln_offset\tmm:i:$nMM\n";
				$lineOut++;
			}

		}			
	
		$i++;
	}
}

close OUT_BEDss1;
close OUT_BEDss2;
close OUT_BEDds;
close OUT_BEDdsS;
close OUT_BEDunC;
close OUT_SAM;
close OUT_SAMss1;
close OUT_SAMss2;
close OUT_SAMds;
close OUT_SAMdsS;
close OUT_SAMunC;


################################################################################
sub find_offset {
# $ref_2,$query_2
my ($seqA, $seqB, $f_read, $strand,$max_mm)=@_;

	if ($strand == -1) {
		$seqA = scalar reverse $seqA;
		$seqB = scalar reverse $seqB;
		$f_read =~ tr/ATGC/TACG/;
	}
	
	my @arrA = split(//,$seqA);
	my @arrB = split(//,$seqB);
	my @arrT = split(//,$f_read);
	my @arrC;
	my $j = 0;
	my $i = 0;
	while ($i<=$#arrT && $j <= $#arrB ) {
		if ($arrB[$i] eq "-") {
			$arrC[$i] = $arrB[$i];
		} else {
			$arrC[$i] = $arrT[$j]; $j++;
		}
		$i++;
	}
	
	$j = 0; my $mm = 0; my $n_gap = 0;
	
	while ($mm <= $max_mm && $j <= $#arrC) {
		$mm++ if $arrB[$j] ne $arrC[$j];
		$n_gap++ if $arrB[$j] eq "-";
		$j++;
	}
	
	my $aln_offset = $j-1;
	
	while ($arrB[$aln_offset] ne $arrC[$aln_offset]) {$aln_offset--;}
	
	my $seq_offset = $aln_offset - $n_gap + 1;
	
	$seq_offset = 0 if $seq_offset < 0;
	
	my @cmpAB; my @match_pat; my %score; my %rev_score;  my $max_score = 0;
	
	for my $i(0..$#arrA){ 
		$cmpAB[$i] = ($arrA[$i] eq $arrB[$i] || $arrA[$i] eq "N" || $arrB[$i] eq "N") ? 1 : -1; 
		$cmpAB[$i] = -2 if ($arrA[$i] ne $arrB[$i] && $arrB[$i] eq $arrC[$i]);
	}
	
	my $min_i = 0;      while ($min_i < $#arrA && $cmpAB[$min_i+1] == -1) {$min_i++;} # min start of offset 0-based
	my $max_i = $#arrA; while ($max_i > 0      && $cmpAB[$max_i-1] ==  1) {$max_i--;}
	
	for my $i(0..$#arrA){ $match_pat[$i] = $i >= $min_i ? 1 : -1;} # 
	foreach my $offset ($min_i..$max_i) {
		foreach my $i ($min_i..$max_i) {$score{$offset}+=$cmpAB[$i]*$match_pat[$i];}
		push @{$rev_score{$score{$offset}}}, $offset;
		$max_score = $score{$offset} if $score{$offset}>$max_score;
		$match_pat[$offset] = -1;
	}
	
	my $max_offset = 0;
	
	foreach my $offset (@{$rev_score{$max_score}}) {$max_offset = $max_offset<$offset? $offset : $max_offset;}
	my $aln_correction = $max_offset; 
	
	for my $i(0..$max_offset){$aln_correction-- if $arrA[$i] eq "-";}
	
	# ITR length, aln offset, uH
	my $uH_len = $aln_offset-$max_offset+1; $uH_len = 0 if $uH_len < 0;
	return $seq_offset, $max_offset, $uH_len ,$aln_correction;
}

################################################################################
sub parseSAMflag{
	my $flag =shift;
	
	my $B = dec2bin($flag);
	
	while (length($B) < 11){$B = '0'.$B}
	
	$B = reverse($B);

	my %r1;
	
	($r1{'Read_Paired'}, $r1{'Pair_OK'}, $r1{'R1Unmapped'}, $r1{'R2Unmapped'},
	 $r1{'R1RC'}, $r1{'R2RC'}, $r1{'R1'}, $r1{'R2'}, $r1{'SecondaryAln'}, $r1{'QCfail'}, $r1{'Duplicate'}) = split(//,$B);
	
	return %r1;
	
}

################################################################################
sub SAMstrand{
    my $ss = shift;
    my %f  = parseSAMflag($ss);
    my $strand = 1;
    $strand = -1 if ($f{'R1RC'} && $f{'R1'});
    $strand = -1 if ($f{'R1RC'} && $f{'R2'});
    
    return $strand;
}

################################################################################
sub generateSAMline{
    my ($nRead,$rln,$itrLen,$iStrand,$numI,$numD,$numS) = @_;
        
    my ($Ibefore,$Iafter) = (0,0);
    my ($Dbefore,$Dafter) = (0,0);
    my ($Sbefore,$Safter) = (0,0);
    
    if ($nRead == 2){
	
	$Ibefore += $1 while($$rln[5] =~ /(\d+)I/g);
	$Dbefore += $1 while($$rln[5] =~ /(\d+)D/g);
        $Sbefore += $1 while($$rln[5] =~ /(\d+)S/g);
	
	if ($iStrand<0){
            $$rln[5] = trimCigar($$rln[5],$itrLen,'right');
            $$rln[9] =~ s/(\S+)\S{$itrLen}/$1/;
            $$rln[10] =~ s/(\S+)\S{$itrLen}/$1/;
        }else{
            $$rln[5] = trimCigar($$rln[5],$itrLen,'left');
            $$rln[9] =~ s/\S{$itrLen}(\S+)/$1/;
            $$rln[10] =~ s/\S{$itrLen}(\S+)/$1/;
        }
	
	$Iafter += $1 while($$rln[5] =~ /(\d+)I/g);
	$Dafter += $1 while($$rln[5] =~ /(\d+)D/g);
	$Safter += $1 while($$rln[5] =~ /(\d+)S/g);
    }

    my $nD = $numD?$numD:($Dbefore - $Dafter);
    my $nI = $numI?$numI:($Ibefore - $Iafter);
    my $nS = $numS?$numS:($Sbefore - $Safter);
    
    return '' if ($$rln[5] !~ /[MIDN]/);

    if ($iStrand < 0 && $nRead == 1){
        $$rln[7] += $itrLen-$nI+$nD-$nS;
        $$rln[8] += $itrLen-$nI+$nD-$nS;
    }

    if ($iStrand > 0 && $nRead == 2){
        $$rln[3] += $itrLen-$nI+$nD-$nS;
        $$rln[8] -= $itrLen-$nI+$nD-$nS;
    }
    
    my $retLn;
    
    for my $v(@{$rln}){
        my $vOLD = $v;
        
        if ($nRead == 2){
            next if ($v =~ /NM:i:/);
            
            my ($NM,$MDtrim);
            
            if ($v =~ s/(MD:\S+:)(\S+)/$1/){
                my $txt = $1;
                $MDtrim = trimXM($2,$itrLen-$nI+$nD-$nS,'left')  if ($iStrand > 0);
                $MDtrim = trimXM($2,$itrLen-$nI+$nD-$nS,'right') if ($iStrand < 0);
                $v = $txt.$MDtrim;
                $NM++ while ($MDtrim =~ /[GATC]/g); 
                $vOLD =~ s/MD/YK/;
                $v = join("\t","NM:i:".($NM?$NM:"0"),$v,$vOLD);
            }
        }
        $retLn .= $v."\t";
    }
    
    chop($retLn);
    
    return ($retLn,$nI,$nD,$nS);    
}

################################################################################
sub cigar2str{
    my ($cigSTR) = shift;
    my $s;
    while ($cigSTR =~ s/(\d+)(\S)//){
        $s .= "$2" x $1;
    }
    return $s;
}

################################################################################
sub str2cigar{
    my ($STR) = shift;
    my ($s,$cig);
    my $prev = '';
    
    while ($STR =~ s/^(\S)//){
        if ($1 eq $prev){
            $cig .= $1;
        }else{
            $s .= length($cig).$prev if ($prev);
            $cig = $1;
        }
        $prev = $1;
    }
    
    $s .= length($cig).$prev if ($prev);
    return $s;
}

################################################################################
sub trimCigar{
    my ($cigarString,$nTrim,$side) = @_;
    
    my $cigStr = cigar2str($cigarString);
    my @cigArr = split(//,$cigStr);
    my $retStr;
    
    @cigArr = reverse(@cigArr) if ($side eq 'right');
    
    for my $cS(@cigArr){
	$nTrim-- unless ($cS eq 'D');
	$retStr .= $cS if ($nTrim < 0);
    }
    
    $retStr = reverse($retStr) if ($side eq 'right');
    
    return str2cigar($retStr);
}


################################################################################
sub XM2str{
    my ($xmSTR) = shift;
    my $s;
    while ($xmSTR =~ s/(\d+|[GATC]|\^[GATC]+)//){
        my $n = $1;
        next if ($n eq '0');
	if ($n =~ /\d+/){$s .= "M" x $n; next}
	$s .= $n;
        
    }
    return $s;
}

################################################################################
sub str2XM{
    my ($STR) = shift;
    my ($s,$xm);
    my $prev = '';
    
    my $cnt = 0;
    while ($STR =~ s/^(M|[GATC]|\^[GATC]+)//){
        my $x = $1;
	
	if ($x eq 'M'){
            $cnt++; next;
        }else{
            $s .= ($cnt>0?$cnt:'0').$x;
            $cnt = 0;
        }
    }        
    
    $s .= $cnt if ($cnt > 0);
    
    return $s;
}

################################################################################
sub trimXM{
    my ($xmString,$nTrim,$side) = @_;
    
    my $xmStr = XM2str($xmString);
    my @XMArr = split(//,$xmStr);
    my $retStr;
    
    @XMArr = reverse(@XMArr) if ($side eq 'left');
    
    my $XMdel = 0;
    
    $XMdel += (length($1)+1) while ($xmStr =~ /\^([GATC]+)/g);
    my $XMLen = length($xmStr) - $XMdel;
    my $XMOKsz = $XMLen - $nTrim;
    
    my ($noTrim,$nGood);
    
    for my $cS(@XMArr){
	$noTrim = 1 if ($cS eq '^');
	$noTrim = 0 if ($cS eq 'M');
	$nGood++ unless ($noTrim);
	$retStr .= $cS ;
	last if ($nGood == $XMOKsz);
    }
    
    $retStr = reverse($retStr) if ($side eq 'left');
    
    return str2XM($retStr);
}

################################################################################
sub bin2dec {
	return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

################################################################################
sub dec2bin {
	my $str = unpack("B32", pack("N", shift));
	$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros
	return $str;
}
