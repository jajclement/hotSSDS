use strict;
use Bio::DB::Sam;
use Getopt::Long;

# Adapted from ITR_id_v2c_NextFlow2.pl (K. Brick, original SSDS pipeline, 2012).
# That script assumed R2 was aligned with bwa aln + a custom "bwa-ra" binary
# that forces an alignment through the non-homologous 5' ITR fill-in part as
# mismatches. This version instead reads BAM files produced by plain bwa mem,
# where the ITR is represented as a soft-clip (long ITR) and/or leading
# mismatches in the CIGAR M block (short ITR, when the clip penalty exceeds
# the mismatch penalty). Only find_offset() is replaced (see find_offset_mem
# below) -- classification thresholds and SAM/BED output format are unchanged
# so downstream processes (parseITRs in main.nf) do not need to change.

my ($in_bam,$in_fa) = @ARGV;
die ("Invalid genomes path ...\n$in_fa does not exist ...") unless (-e ($in_fa));

my $outname = $in_bam;

## SET UP TEMP FILES FOR EACH SPLIT BAM
my $splitBEDss1 	= $outname.'.ssDNA_type1.bed';
my $splitBEDss2		= $outname.'.ssDNA_type2.bed';
my $splitBEDds		= $outname.'.dsDNA.bed';
my $splitBEDdsS		= $outname.'.dsDNA_strict.bed';
my $splitBEDunC		= $outname.'.unclassified.bed';

my $splitSAMss1		= $outname.'.ssDNA_type1.sam';
my $splitSAMss2		= $outname.'.ssDNA_type2.sam';
my $splitSAMds		= $outname.'.dsDNA.sam';
my $splitSAMdsS		= $outname.'.dsDNA_strict.sam';
my $splitSAMunC		= $outname.'.unclassified.sam';

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
		next if ($left->get_tag_values('QC_FAILED') || $right->get_tag_values('QC_FAILED'));

		my $first_mate; my $second_mate;

		if ($left->get_tag_values('FIRST_MATE')) {
			$first_mate = $left; $second_mate = $right;
		} else {
			$first_mate = $right; $second_mate = $left;
		}

		my @s_tam = split(/\t/,$second_mate->tam_line());
		my @f_tam = split(/\t/,$first_mate->tam_line());

		my ($ITR_len, $aln_offset, $uH_len, $aln_correction) =
		    find_offset_mem(\@s_tam,\@f_tam,$second_mate->strand,1);

		# Guard: for a short fragment, $aln_correction (trimmed off both start and end,
		# see below) can exceed half the fragment span and invert start/end into a
		# malformed BED interval (start > end). Clamp it so the reported span never
		# collapses to zero or negative.
		my $frag_span = $right->end - $left->start;
		my $max_correction = $frag_span > 1 ? int(($frag_span - 1) / 2) : 0;
		$aln_correction = $max_correction if $aln_correction > $max_correction;

		# generate output;
		my ($strand,$start,$end,$chrom,$f_qual,$s_qual);

		$strand = ($first_mate->strand == 1)? "+": "-";
		$start  = ($first_mate->strand == 1)? $left->start - 1 : $left->start + $aln_correction - 1;
		$end    = ($first_mate->strand == 1)? $right->end - $aln_correction: $right->end;
		$chrom  = $first_mate->seq_id;
		$f_qual = $first_mate->qual;
		$s_qual = $second_mate->qual;

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
close OUT_SAMss1;
close OUT_SAMss2;
close OUT_SAMds;
close OUT_SAMdsS;
close OUT_SAMunC;


################################################################################
# parse_cigar_tokens : "10S90M2D40M" -> ([10,'S'],[90,'M'],[2,'D'],[40,'M'])
# in the order the CIGAR is stored (i.e. reference-forward order, NOT
# necessarily the read's own biological 5'->3' order -- see find_offset_mem).
################################################################################
sub parse_cigar_tokens {
	my ($cig) = @_;
	my @tok;
	while ($cig =~ /(\d+)([MIDNSHP=X])/g) { push @tok, [$1,$2]; }
	return @tok;
}

################################################################################
# find_offset_mem : replaces the legacy find_offset() for bwa-mem-aligned R2.
#
# bwa mem soft-clips the 5' ITR fill-in of R2 instead of forcing it through
# bwa aln's edit-distance model (that was the job of the modified bwa-ra
# binary). For a short ITR, bwa mem may instead leave it as leading
# mismatches inside the CIGAR M block (soft-clip penalty > mismatch penalty),
# so the candidate ITR region is read off *both* the CIGAR soft-clip and the
# leading mismatches in the MD tag immediately following it.
#
# The fill-in part of that candidate region is then measured by comparing it
# to the complement of R1's own 5' end (Khil et al. 2012: "The fill-in part
# of the ITR must be complementary to the 5'-end of the first read but not
# homologous to the target genomic DNA"), tolerating $max_mm mismatches while
# extending from the boundary. Agreement beyond that point (microhomology) is
# recorded separately.
#
# Known limitation (flag for validation against real data, see plan Part D):
# a CIGAR insertion landing exactly at the ITR/genome boundary is not
# accounted for when walking the MD tag (MD does not encode insertions).
#
# Returns ($ITR_len, $aln_offset, $uH_len, $aln_correction) with the same
# roles as the legacy find_offset(), so downstream classification and
# generateSAMline() are unchanged.
################################################################################
sub find_offset_mem {
	my ($s_tam, $f_tam, $strand, $max_mm) = @_;

	my $cigar2 = $$s_tam[5];
	my $seq2   = $$s_tam[9];

	return (0,0,0,0) if ($cigar2 eq '*' || $seq2 eq '*');

	my $md2 = '';
	for my $idx (11..$#$s_tam) { $md2 = $1 if ($$s_tam[$idx] =~ /^MD:Z:(\S+)/); }

	# Orient R2's CIGAR tokens and sequence to its own biological 5'->3' order
	my @tok = parse_cigar_tokens($cigar2);
	@tok = reverse(@tok) if ($strand == -1);

	my $oseq = $seq2;
	if ($strand == -1) {
		$oseq = reverse($oseq);
		$oseq =~ tr/ACGTacgt/TGCAtgca/;
	}

	# 1) Leading soft-clip length (5'-of-R2-oriented)
	my $clip_len = (@tok && $tok[0][1] eq 'S') ? $tok[0][0] : 0;

	# 2) Leading mismatches from MD, immediately after the clip
	my @mdtok;
	my $md_work = $md2;
	while (length($md_work)) {
		if    ($md_work =~ s/^(\d+)//)         { push @mdtok, ['M',$1] if $1 > 0; }
		elsif ($md_work =~ s/^\^([A-Za-z]+)//) { push @mdtok, ['D',$1]; }
		elsif ($md_work =~ s/^([A-Za-z])//)    { push @mdtok, ['X',$1]; }
		else  { last; }
	}
	@mdtok = reverse(@mdtok) if ($strand == -1);

	my $lead_mm = 0;
	for my $t (@mdtok) {
		last unless ($t->[0] eq 'X');
		$lead_mm++;
	}

	my $ITR_len = $clip_len + $lead_mm;

	# 3) Compare R2's candidate ITR bases to R1's own 5' end, in R1's own 5'->3' order.
	# Empirically verified against real bwa-mem-aligned read pairs (see
	# CHANGELOG_bwa_mem_migration.md): a genuine fill-in shows up as a direct base-for-base
	# match here, not a complement -- an earlier version of this code applied one extra/missing
	# complement operation, which capped aln_offset near 0 for essentially every pair.
	my $r1seq = $$f_tam[9];
	my $r1oriented;
	if ($strand == -1) {
		# R2 on '-' strand => R1 on '+' strand: SEQ is already R1's own 5'->3'; no transform needed.
		$r1oriented = $r1seq;
	} else {
		# R2 on '+' strand => R1 on '-' strand: SEQ is R1's reverse-complement;
		# reverse-complement it back to R1's own 5'->3' order.
		($r1oriented = reverse($r1seq)) =~ tr/ACGTacgt/TGCAtgca/;
	}

	my @r2 = split(//, substr($oseq, 0, $ITR_len));
	my @r1 = split(//, $r1oriented);

	my $j = 0; my $mm = 0;
	while ($j <= $#r2 && $j <= $#r1 && $mm <= $max_mm) {
		$mm++ if ($r2[$j] ne $r1[$j]);
		$j++;
	}
	my $aln_offset = $j - 1;
	$aln_offset-- while ($aln_offset >= 0 && $r2[$aln_offset] ne $r1[$aln_offset]);
	$aln_offset++;
	$aln_offset = 0 if ($aln_offset < 0);

	# 4) Microhomology: extra agreement beyond aln_offset (ambiguous zone)
	my $uH_len = 0;
	while (($aln_offset + $uH_len) <= $#r2 && ($aln_offset + $uH_len) <= $#r1
	       && $r2[$aln_offset + $uH_len] eq $r1[$aln_offset + $uH_len]) {
		$uH_len++;
	}

	my $aln_correction = $ITR_len;

	return ($ITR_len, $aln_offset, $uH_len, $aln_correction);
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
    if ($itrLen >=0) {
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
