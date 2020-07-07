#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;
use Getopt::Long;
use Math::Round; 

GetOptions ('bam=s' 	    => \(my $bam),
	        'v+'            => \(my $verbose),
            'g=s'           => \(my $g),
            'gIdx=s'        => \(my $gIdx),
            'gw=s'          => \(my $gwFile),
            's=i'           => \(my $step = 100),
            'w=i'           => \(my $win = 1000),
            'o=s'           => \(my $out = 'BAM2BW.out'),
            'z+'            => \(my $showZeros), ## print Zero value intervals
            'keepBG+'       => \(my $keepBG), ## print Zero value intervals
            'p+'            => \(my $pointIntervals), ## use midpoints for BG and BW interval definitions
            'revStrand+'    => \(my $reverseStrands), ## REVERSE Strand order
            'noSortedChk+'  => \(my $noSortedCheck));

if (`samtools view $bam |head -n 100 |wc -l` < 100){
	my $bw = 'EMPTY_'.$bam.'.bigwig';
	sysAndPrint("touch $bw");
	exit
}

my $tmpFolder = $ENV{SCRATCH}.'/';
my $tmpID     = int(rand()*100000000000000000);
my $tmpFileBase = $tmpFolder.$tmpID;

my $in = $tmpFileBase.'_1';
bamToTempBed($bam,$in,$tmpFileBase);

$gIdx = $ENV{NXF_GENOMES}.'/'.$g.'/genome.fa.fai';

## Get or generate genomic windows
my $newGW    = $tmpFileBase.'_2';
my $newGWtmp = $tmpFileBase.'_2tmp';
my $gw = (($gwFile && (-e $gwFile))?$gwFile:genGenomicWindows($g,$step,$win,$gIdx,$newGWtmp,$newGW));

## Ensure GW file and reads file are ordered the same way
die("Input file ($in) must be sorted with -k1,1 -k2n,2n -k3n,3n ... \n") unless ($noSortedCheck || compareOrder($gw,$in));

## Split Input bed to F/R reads (fragments)
my $fwdFile = $tmpFileBase.'_3F';
my $revFile = $tmpFileBase.'_3R';
my ($readCountTot,$readCountF,$readCountR) = splitFR($in,$fwdFile,$revFile,$reverseStrands);

## Generate file names and open pipes
my ($bedgraphFR,$bedgraphF,$bedgraphR,$bedgraphTot,$bwFR,$bwF,$bwR,$bwTot) = genFileNames($out);

open FRBG,  '>', $bedgraphFR;
open FBG,   '>', $bedgraphF;
open RBG,   '>', $bedgraphR;
open TOTBG, '>', $bedgraphTot;

## Count reads in intervals
my $ibFwd  = $tmpFileBase.'_4F';
my $ibRev  = $tmpFileBase.'_4R';
my $ibBoth = $tmpFileBase.'_4Both.tab';

#open my $fwdCmd, '-|', "intersectBed -sorted -a $gw -b $fwdFile -c";
#open my $revCmd, '-|', "intersectBed -sorted -a $gw -b $revFile -c";
#while (defined(my $lineF=<$fwdCmd>) && defined(my $lineR=<$revCmd>)){
#    chomp $lineF; chomp $lineR;
#    
#    my ($csF,$fromF,$toF,$valF) = split(/\t/,$lineF);
#    my ($csR,$fromR,$toR,$valR) = split(/\t/,$lineR);
#    
#    my $fName = join(":",$csF,$fromF,$toF);
#    my $rName = join(":",$csR,$fromR,$toR);

sysAndPrint("intersectBed -sorted -a $gw    -b $fwdFile -c >$ibFwd");
sysAndPrint("intersectBed -sorted -a $ibFwd -b $revFile -c >$ibRev");

if ($showZeros){
    sysAndPrint("mv $ibRev $ibBoth");
}else{
    sysAndPrint('grep -vP "^chr\S+\t\d+\t\d+\t0\t0$" '.$ibRev.' >'.$ibBoth);
}

open IB, $ibBoth;

while (<IB>){
    chomp $_;

    ## ENSURE NO DOUBLE SPACES SCREW UP COLUMN ORDER
    $_ =~ s/\s+(\S)/\t$1/g;

    my ($cs,$from,$to,$valF,$valR) = split(/\t/,$_);
       
    my $fName = join(":",$cs,$from,$to);
    my $rName = join(":",$cs,$from,$to);
    
    ## Ensure that the lines agree - they always should
    die unless($fName eq $rName);
    
    ## Ensure intervals don't overlap in BedGraph file
    my $locus = round($from + $win/2);
    
    my $locusF = $pointIntervals?($locus-1):($locus-round($step/2));
    my $locusT = $pointIntervals?$locus:($locus+round($step/2)-1);
    
    my $fpkmFRBG  = inFPKM($valF-$valR,$win,$readCountTot);
    my $fpkmFBG   = inFPKM($valF,$win,$readCountF);
    my $fpkmRBG   = inFPKM($valR,$win,$readCountR);
    my $fpkmTOTBG = inFPKM($valF+$valR,$win,$readCountTot);
    
    ## Print BedGraph outputs
    print FRBG  join("\t",$cs,$locusF,$locusT,$fpkmFRBG)."\n"  unless (not($showZeros) && not($fpkmFRBG));
    print FBG   join("\t",$cs,$locusF,$locusT,$fpkmFBG)."\n"   unless (not($showZeros) && not($fpkmFBG));
    print RBG   join("\t",$cs,$locusF,$locusT,$fpkmRBG)."\n"   unless (not($showZeros) && not($fpkmRBG));
    print TOTBG join("\t",$cs,$locusF,$locusT,$fpkmTOTBG)."\n" unless (not($showZeros) && not($fpkmTOTBG));
}

close FRBG; close FBG; close RBG; close TOTBG;

## Make BigWigs
sysAndPrint("bedGraphToBigWig $bedgraphFR  $gIdx $bwFR");
sysAndPrint("bedGraphToBigWig $bedgraphF   $gIdx $bwF");
sysAndPrint("bedGraphToBigWig $bedgraphR   $gIdx $bwR");
sysAndPrint("bedGraphToBigWig $bedgraphTot $gIdx $bwTot");

my $allOK;
$allOK = 1 if (-e $bwFR && -e $bwF && -e $bwR && -e $bwTot);

## Remove Bedgraph Files
unless ($keepBG){
    if ($allOK){
        sysAndPrint("rm $bedgraphFR");
        sysAndPrint("rm $bedgraphF");
        sysAndPrint("rm $bedgraphR");
        sysAndPrint("rm $bedgraphTot");
        sysAndPrint("rm $tmpFileBase*");
    }else{
        print STDERR "WARNING ********************************************\n";
        print STDERR "Something isn't right ... output bigwigs are missing\n";
        print STDERR "BedGraphs and temp files are retained ... \n";
    }
}

################################################################################
sub bamToTempBed{
    my ($pcBam,$pcBed,$plTmpName) = @_;
    
	my $outTMP = $plTmpName."_btt1";
	open OUT, '>', $outTMP;
	open my $PIPE, '-|', "bedtools bamtobed -i $pcBam";
	while (<$PIPE>){
		chomp; 
		my @F = split(/\t/,$_);
		if ($F[3] =~ /2$/){
			if ($F[5] eq '+'){
				$F[5] = '-';
			}else{
				$F[5] = '+';
			}
		}
		print OUT join("\t",@F)."\n";
	}
	close OUT; 
	sysAndPrint("sort -k1,1 -k2n,2n -k3n,3n $outTMP >$pcBed");
        
    return 1;
}

################################################################################
sub genGenomicWindows{
    my ($gwG,$gwStep,$gwWin,$gwIdx,$gwTmpFile,$gwFile) = @_;
    
    ## get CS sizes from fa.fai
    my %csSize;
    open IDX, $gwIdx;
    
    while (<IDX>){
        my @F = split("\t");
        next if ($F[0] =~ /(Rand|rand|Un|un|hap|chrM|chrY|cont|part)/);
        $csSize{$F[0]} = $F[1];
    }

    open GW, '>', $gwTmpFile; 
    ## Generate Windows
    for my $c (sort keys(%csSize)){
        my $lastF;
        
        for (my $f = 0; $f <= ($csSize{$c}-$gwWin); $f+=$gwStep){
            print GW join("\t",$c,$f,$f+$gwWin-1)."\n";
            $lastF = $f;
        }
    }
    close GW;
    
    sysAndPrint('sort -k1,1 -k2n,2n -k3n,3n '.$gwTmpFile.' >'.$gwFile);
    
    return $gwFile;
}

################################################################################
sub compareOrder{
    my ($coWins,$coIn) = @_;
    
    my $coInData  =`sort -k1,1 -k2n,2n -k3n,3n -c $coIn 2>&1`;
    die ("Input reads not sorted ($in) !!\n") if ($coInData);
    
    my $coWinData =`sort -k1,1 -k2n,2n -k3n,3n -c $coWins 2>&1`;
    die ("Window data not sorted !!\n") if ($coWinData);
        
    return 1;
}

################################################################################
sub splitFR{
    my ($sIn,$sFwd,$sRev,$invertStrand) = @_;
    
    my $readCnt  = 0;
    my $readCntF = 0;
    my $readCntR = 0;
    
    open OUTFWD, '>', $sFwd;
    open OUTREV, '>', $sRev;
    
    open IN, $sIn;
    while (<IN>){
	chomp;
	my @Fbed = split(/\t/,$_);
	
	$Fbed[3] = '.';
	$Fbed[4] = '.';

	## KB 180928 - for OKSeq sim
	if ($invertStrand){
		$Fbed[5] =~ s/\+/MINUS/;
		$Fbed[5] =~ s/\-/\+/;
		$Fbed[5] =~ s/MINUS/\-/;
	}

        $readCnt++;
        if ($Fbed[5] =~ /\+/){
            print OUTFWD join("\t",@Fbed)."\n" ;
            $readCntF++;
        }else{
            print OUTREV join("\t",@Fbed)."\n";
            $readCntR++;
        }
    }
    
    close OUTFWD; close OUTREV; close IN;
    return ($readCnt,$readCntF,$readCntR);
}

################################################################################
sub genFileNames{
    my ($outName) = shift;
    
    my $outStem;
    
    if ($outName =~ /^(.+)\.\S+$/){
        $outStem = $1; 
    }else{
        $outStem = $outName;
        die if (-e $outName);
    }
    
    my $of_bedgraphFR   = "$outStem\.FR\.bedgraph";
    my $of_bedgraphF    = "$outStem\.F\.bedgraph";
    my $of_bedgraphR    = "$outStem\.R\.bedgraph";
    my $of_bedgraphTot  = "$outStem\.Tot\.bedgraph";
    my $of_bwFR         = "$outStem\.FR\.bigwig";
    my $of_bwF          = "$outStem\.F\.bigwig";
    my $of_bwR          = "$outStem\.R\.bigwig";
    my $of_bwTot        = "$outStem\.Tot\.bigwig";
    
    return ($of_bedgraphFR,$of_bedgraphF,$of_bedgraphR,$of_bedgraphTot,$of_bwFR,$of_bwF,$of_bwR,$of_bwTot);
}

################################################################################
sub inFPKM {
    my ($sc,$sz,$tot) = @_;
    my $fpkm = ((($sc/$sz)*1000) / $tot) * 1000000;
    return $fpkm?sprintf('%4.3f',$fpkm):0;
}


################################################################################
sub sysAndPrint{
	my $c = shift;
	print STDERR $c."\n" if ($verbose);
	system($c);
}
