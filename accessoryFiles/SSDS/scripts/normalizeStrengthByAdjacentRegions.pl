#!/usr/bin/perl
use strict;
use List::Util qw/max min/;
use Getopt::Long;
use Statistics::Descriptive;
use Math::Round;
use ssPipeline_biowulf;

GetOptions ('bed=s' 	=> \(my $inbed),
            'in=s'  	=> \(my $fragFile),
	    	'input=s'   => \(my $inputFile),
	    	'NCIS=s'   	=> \(my $NCIS),
	    	'v+'    	=> \(my $verbose),
	    	'rc+'  	=> \(my $recenter),
	    	'nt+'   	=> \(my $noTmp),
	    	'ac+'   	=> \(my $altCenter),
	    	'noCheck+'	=> \(my $noCheck),
			'revFrag+' 	=> \(my $reverseFragments),
			'SNS+'   	=> \(my $forSNS),
	    	'h+'    	=> \(my $help),
	    	'help+' 	=> \(my $helpMe),
            'out=s' 	=> \(my $out));

## FOR SNS, reverse fragments
$reverseFragments++ if ($forSNS);

## Show Help
if ($help || $helpMe){
	showHELP();
	exit;
}

die("No Good ($inbed)\n") unless ($inbed);
die("No Good ($inbed)\n") unless (-e $inbed);

unless ($out){$out = $inbed; $out =~ s/^(.+)\.\S+$/$1.normStrength.out/};

my $outTab 	= $out; $outTab =~ s/^(.+)\.\S+$/$1.tab/;
my $outBG  	= $out; $outBG =~ s/^(.+)\.\S+$/$1.bedgraph/;
my $outBG  	= $out; $outBG =~ s/^(.+)\.\S+$/$1.bedgraph/;
my $outTab  	= $out; $outTab =~ s/^(.+)\.\S+$/$1.tab/;

my (%hotspots,$NCISfactor);
my ($tempFlanks,$hFlanks)          = genTempFile();
my ($tempHS,$temphHS)              = genTempFile();
my ($fragstmpFile,$tempFragsFileH) = genTempFile();

$tempFlanks   = 'tempFlanks.bed'    if ($noTmp);
$tempHS       = 'tempHS.bed'        if ($noTmp);
$fragstmpFile = 'tempFragsFile.bed' if ($noTmp);

## FOR SNS
if ($reverseFragments){
	sysAndPrint('perl '.$ENV{NXF_PIPEDIR}.'/reverseStrandsForOriCalling.pl <'.$fragFile.'>'.$fragstmpFile);
}else{
	sysAndPrint("cp $fragFile $fragstmpFile");
}

$fragFile = $fragstmpFile;

my $lc = lineCount($fragFile);
$NCISfactor = getNCISfactor($NCIS) if ($NCIS && $inputFile);

## check that peaks and frags are sorted by -k1,1 -k2n,2n
## if not, the $inbed and $fragFile are directed to temp files with sorted data
my ($tfX1,$tfhX1) = genTempFile();
my ($tfX2,$tfhX2) = genTempFile();
($inbed,$fragFile) = checkAndSort($inbed,$fragFile,$tfX1,$tfX2) unless ($noCheck);

if ($recenter){
	recenterPeaks($tempHS,$inbed,$fragFile);
}else{
	system("cut -f1-3 $inbed >$tempHS");
}

makeFlanksFile($tempHS,$tempFlanks);
genCoverage($tempFlanks,$fragFile,\%hotspots);

printHS(\%hotspots,$lc);

################################################################################
sub checkAndSort{
    my ($coPeaks,$coTags,$tfPk,$tfTag) = @_;
    
    my $coPeakData  =`sort -k1,1 -k2n,2n -c $coPeaks 2>&1`;
    if ($coPeakData){
	sysAndPrint("sort -k1,1 -k2n,2n $coPeaks >$tfPk",1);
	$coPeaks = $tfPk;
    }
    
    my $coTagData =`sort -k1,1 -k2n,2n -c $coTags 2>&1`;
    if ($coTagData){
	sysAndPrint("sort -k1,1 -k2n,2n $coTags >$tfTag",1);
	$coTags = $tfTag;
    }

    return($coPeaks,$coTags);
}
################################################################################
sub recenterPeaks{
    my ($outFile,$rcPeaks,$rcTags,$rcType) = @_;
    
    my ($tempOut,$tempOutFh) = genTempFile();
    $tempOut = 'tempOut.bed' if ($noTmp);
    open OUT, '>', $tempOut;

    my ($tempHS3Col,$temphHS3Col) = genTempFile();
    $tempHS3Col = 't3Col.bed' if ($noTmp);

    system("cut -f1-3 $rcPeaks >$tempHS3Col");
    
    my (@hsArr,%allHSCheck); 
    open TMPHS, $tempHS3Col;
    while (<TMPHS>){
	chomp;
	$_ =~ s/(\s|\t)+(\S)/:$2/g;
	push @hsArr, join(":",$_); 
	$allHSCheck{join(":",$_)}->{val}++; 
    }
    close TMPHS; 

    #my $cmd = "intersectBed -a $rcTags -b $tempHS3Col -wo |sort -k7,7 -k8rn,8rn -k2n,2n ";
    my $cmd = "intersectBed -sorted -a $rcTags -b $tempHS3Col -wo";
    open my $PIPE, '-|', $cmd;
    
    my (%lastHS,%hs,@pos,@neg);
    my ($cs,$from,$to,%hsDone);
    
    while (<$PIPE>){
	
        my @F = split(/\t/,$_);
        ($cs,$from,$to) = ($F[0],$F[1],$F[2]);
	
	## control for some odd reads.
	next if ($F[2]-$F[1] > 1000);
	
        $hs{name} = join(":",$F[6],$F[7],$F[8]);
	$hs{cs}   = $F[6];
	$hs{from} = $F[7];
	$hs{to}   = $F[8];
	
        if ($lastHS{name}){
            if ($hs{name} ne $lastHS{name}){
	        printRecenteredHS(\@pos,\@neg,$lastHS{cs},$lastHS{from},$lastHS{to}) unless ($hsDone{$hs{name}}++);
		$allHSCheck{$lastHS{name}}->{'OK'}++;
                @pos = (); @neg = ();                
            }
            
	    my $arr = (($F[5] eq '+')?\@pos:\@neg);

            for my $nt(max($F[1],$hs{from})..min($F[2],$hs{to})){
                push @$arr, $nt;
            }    
        }
        %lastHS = %hs;
    }

    printRecenteredHS(\@pos,\@neg,$lastHS{cs},$lastHS{from},$lastHS{to},$altCenter);
    $allHSCheck{$lastHS{name}}->{'OK'}++;

    for my $kH (keys(%allHSCheck)){
	my ($kCS,$kFrom,$kTo) = split(/:/,$kH);
	print OUT join("\t",$kCS,$kFrom,$kTo,"ZERO")."\n" unless ($allHSCheck{$kH}->{'OK'});
    }

    close OUT; 
    system ("sort -k1,1 -k2n,2n $tempOut |cut -f1-4 >$outFile");
}  

################################################################################
sub printRecenteredHS{
    my ($pPos,$pNeg,$pCS,$pFrom,$pTo,$altCenter) = @_;

    #my $posMed = $statPos->median;
    #my $negMed = $statNeg->median;

    my $posMed = median($pPos);
    my $negMed = median($pNeg);
    
	if ($pCS eq "chr12" && $pTo >= 110844691){
		my $rr = 1;
	}
    my $midPoint;
    if ($posMed && $negMed && not($altCenter)){
	$midPoint = round(($posMed + $negMed)/2);
    }else{
   	my $statPos = Statistics::Descriptive::Full->new();
   	my $statNeg = Statistics::Descriptive::Full->new();
    	$statPos->add_data(@$pPos) if (@$pPos);
    	$statNeg->add_data(@$pNeg) if (@$pNeg);

	if ($posMed && $negMed && $altCenter){$midPoint = round(($statPos->percentile(80) + $statNeg->percentile(20))/2)};
    	if ($posMed && not($negMed)){$midPoint = round($statPos->percentile(80))};
    	if ($negMed && not($posMed)){$midPoint = round($statNeg->percentile(20))};
    }
     
    my $width    = max(abs($midPoint-$pFrom),abs($midPoint-$pTo));
       
    #my $warning  = ($statPos->percentile(80) > $statNeg->percentile(20))?"PgtN":"OK";
    my $warning = "na";
    print OUT join("\t",$pCS,$midPoint-$width,$midPoint+$width,$warning)."\n";
}
################################################################################
sub makeFlanksFile{
	
	my ($bed,$tfFinal) = @_;
	
	my ($tf1,$tfH1) = genTempFile();
	my ($tf2,$tfH2) = genTempFile();
	
	($tf1,$tf2) = ('tf1.bed','tf2.bed') if ($noTmp);
	
	open TMPFLANKS, ">", $tf1;
	open TMPHS, ">", $tf2;
	
	open IN, $bed;
	
	while (<IN>){
		chomp;
		my @F = split(/\t/,$_);
		my $midPt     	= round(($F[1]+$F[2])/2);
		my $halfWidth  	= round(($F[2]-$F[1])/2);
		my $lastXpc  	= round($halfWidth*.4);
		
		my $fName  	= join(":",$F[0],$F[1],$F[2],"F");
		my $hsName 	= join(":",$F[0],$F[1],$F[2],"HS");
		my $rhsName 	= join(":",$F[0],$F[1],$F[2],"RHS");
		my $lhsName 	= join(":",$F[0],$F[1],$F[2],"LHS");
		my $rhsPCName 	= join(":",$F[0],$F[1],$F[2],"RHSPC");
		my $lhsPCName 	= join(":",$F[0],$F[1],$F[2],"LHSPC");
		
		print TMPFLANKS   printLocus($F[0],$F[1],$midPt,$lhsName);
		print TMPFLANKS   printLocus($F[0],$midPt+1,$F[2],$rhsName);
		print TMPFLANKS   printLocus($F[0],$F[1],$F[1]+$lastXpc,$lhsPCName);
		print TMPFLANKS   printLocus($F[0],$F[2]-$lastXpc,$F[2],$rhsPCName);
		print TMPFLANKS   printLocus($F[0],$F[1]-10001,$F[1]-1,$fName);
		print TMPFLANKS   printLocus($F[0],$F[2]+1,$F[2]+10001,$fName);
		print TMPHS 	  printLocus($F[0],$F[1],$F[2],$hsName);
	}
	
	close IN; close TMPFLANKS; close TMPHS;
	
	##################### STEP 2
	my ($tf3,$tfH3) = genTempFile();
	$tf3 = 'tf3.bed' if ($noTmp);
	
	open TMPNOHS, ">", $tf3;
	
	my $cmd = 'intersectBed -a '.$tf1.' -b '.$bed.' -wao';
	open my $PIPE, '-|', $cmd;
	
	while (<$PIPE>){
		my ($cs,$from,$to,$name,$hsCS,$hsFrom,$hsTo,@others) = split(/\t/,$_);
		
		## Keep ALL Left / Right Hotspot Halves
		if ($name =~ /[LR]HS/){
			print TMPNOHS join("\t",$cs,$from,$to,$name)."\n" ;
			next;
		}
		
		if ($hsFrom eq '-1'){
			print TMPNOHS join("\t",$cs,$from,$to,$name)."\n" ;
		}else{
			next if ($hsFrom <= $from && $hsTo >= $to);
			
			if ($hsFrom <= $from && $hsTo < $to){
				print TMPNOHS join("\t",$cs,$hsTo+1,$to,$name)."\n" ;
			}
			
			if ($hsFrom > $from && $hsTo >= $to){
				print TMPNOHS join("\t",$cs,$from,$hsFrom-1,$name)."\n" ;
			}
			
			if ($hsFrom > $from && $hsTo < $to){
				print TMPNOHS join("\t",$cs,$from,$hsFrom-1,$name)."\n" ;
				print TMPNOHS join("\t",$cs,$hsTo+1,$to,$name)."\n" ;
			}
		}
	}
	
	close TMPNOHS;
	
	sysAndPrint("sort -k1,1 -k2n,2n $tf2 $tf3 >$tfFinal",1);
}

################################################################################
sub genCoverage{
	my ($flankFile,$coverageFile,$HS) = @_;
	
	my ($coverageFileF,$tfH1) = genTempFile();
	my ($coverageFileR,$tfH2) = genTempFile();
	
	genFR($coverageFile,$coverageFileF,$coverageFileR);
	
	my %cmd;
	#$cmd{fwd} = 'coverageBed -a '.$coverageFileF.' -b '.$flankFile;
	#$cmd{rev} = 'coverageBed -a '.$coverageFileR.' -b '.$flankFile;
	
	$cmd{fwd} = 'intersectBed -sorted -b '.$coverageFileF.' -a '.$flankFile.' -c';
	$cmd{rev} = 'intersectBed -sorted -b '.$coverageFileR.' -a '.$flankFile.' -c';
	
	$cmd{input} = 'intersectBed -sorted -b '.$inputFile.' -a '.$flankFile.' -c';
	
	for my $strand ('fwd','rev','input'){

		next if ($strand eq 'input' && not ($inputFile));		

		open my $PIPE, '-|', $cmd{$strand};
		
		while (<$PIPE>){	
			chomp;
			#my ($cs,$from,$to,$name,$cover,$nBPCover,$lenB,$pcBPCover) = split(/\t/,$_);
			my ($cs,$from,$to,$name,$cover) = split(/\t/,$_);
			
			$name =~ s/:(F|HS|LHS|RHS|LHSPC|RHSPC)$//; my $type = $1;
			
			my $midPt     	= round(($from+$to)/2);
			my $width 	= $to-$from;

			if ($type eq 'HS'){
				$$HS{$name}->{cover}->{$strand} = $cover;
				$$HS{$name}->{cover}->{total}   += $cover;
				$$HS{$name}->{size}		= $width;
				$$HS{$name}->{cs}		= $cs;
				$$HS{$name}->{from}		= $from;
				$$HS{$name}->{to}		= $to;
				$$HS{$name}->{totinput}		= $cover if ($strand eq 'input');
				$$HS{$name}->{NCISinput}	= $cover*$NCISfactor if ($strand eq 'input');
			}
			
			if ($type eq 'F'){
				$$HS{$name}->{flankcover}->{input}   += $cover*$NCISfactor if ($strand eq 'input');
				$$HS{$name}->{flankcover}->{$strand} = $cover;
				$$HS{$name}->{flankcover}->{total}   += $cover;
				$$HS{$name}->{flanksize}  	     += $width;
			}
			
			next if ($strand eq 'input');
			
			if ($type eq 'LHS'){
				my $type = ($strand eq 'fwd')?'signal':'noise';
				$$HS{$name}->{LHS}->{$type}        = $cover;
				$$HS{$name}->{LHSsize}  	   = $width;
			}

			if ($type eq 'RHS'){
				my $type = ($strand eq 'rev')?'signal':'noise';
				$$HS{$name}->{RHS}->{$type}        = $cover;
				$$HS{$name}->{RHSsize}  	   = $width;
			}
			
			if ($type eq 'LHSPC'){
				my $type = ($strand eq 'fwd')?'signal':'noise';
				$$HS{$name}->{LHSPC}->{$type}      = $cover;
				$$HS{$name}->{LHSPCsize}  	   = $width;
			}

			if ($type eq 'RHSPC'){
				my $type = ($strand eq 'rev')?'signal':'noise';
				$$HS{$name}->{RHSPC}->{$type}      = $cover;
				$$HS{$name}->{RHSPCsize}  	   = $width;
			}
			
		}
	}
}

################################################################################
sub printHS{
	my ($HS,$fragCnt) = @_;
	
	my ($prelimOut,$poH) = genTempFile();
	open TMPFINAL, '>', $prelimOut;
	my ($totStrength,$NCISfactor);
	
	#$NCISfactor = getNCISfactor($NCIS) if ($inputFile && $NCIS);
	
	for my $name (sort (keys %{$HS})){
		next unless ($$HS{$name}->{size} >0);

		my ($adjacentNCISCover,$inputCover,$TOTinputCover);
		my ($adjacentNCISCover_perNT,$inputCover_perNT,$TOTinputCover_perNT);
		my ($adjacentNCISCover_HS,$inputCover_HS,$TOTinputCover_HS);
			
		if ($inputFile && $NCIS){
			$adjacentNCISCover  	 = ($$HS{$name}->{flankcover}->{input});
			$inputCover  	 	 = ($$HS{$name}->{NCISinput});
			$TOTinputCover  	 = ($$HS{$name}->{input});
			
			$adjacentNCISCover_perNT = $$HS{$name}->{flanksize}?$adjacentNCISCover/$$HS{$name}->{flanksize}:$adjacentNCISCover;
			$inputCover_perNT 	 = $inputCover/$$HS{$name}->{size};
			$TOTinputCover_perNT  	 = $TOTinputCover/$$HS{$name}->{size};
			
			$adjacentNCISCover_HS  	 = $adjacentNCISCover_perNT*$$HS{$name}->{size};
			$inputCover_HS   	 = $inputCover_perNT*$$HS{$name}->{size};
			$TOTinputCover_HS  	 = $TOTinputCover_perNT*$$HS{$name}->{size};
		}
		
		my $cover           	 = $$HS{$name}->{cover}->{total};
		my $adjacentCover  	 = ($$HS{$name}->{flankcover}->{total});
		my $antiSenseBG    	 = ($$HS{$name}->{RHS}->{noise} + $$HS{$name}->{LHS}->{noise});
		my $antiSensePCBG    	 = divZero(($$HS{$name}->{RHSPC}->{noise} + $$HS{$name}->{LHSPC}->{noise}),($$HS{$name}->{RHSPCsize} + $$HS{$name}->{LHSPCsize}))*$$HS{$name}->{size};
		my $leftSignal    	 = $$HS{$name}->{LHS}->{signal};
		my $rightSignal    	 = $$HS{$name}->{RHS}->{signal};
		my $leftNoise    	 = $$HS{$name}->{LHS}->{noise};
		my $rightNoise    	 = $$HS{$name}->{RHS}->{noise};

		my $leftFwd	    	 = $$HS{$name}->{LHS}->{signal};
		my $rightRev    	 = $$HS{$name}->{RHS}->{signal};
		my $leftRev	    	 = $$HS{$name}->{LHS}->{noise};
		my $rightFwd    	 = $$HS{$name}->{RHS}->{noise};		

		my $cover_perNT          = divZero($cover,$$HS{$name}->{size});
		my $adjacentCover_perNT  = $$HS{$name}->{flanksize}?divZero($adjacentCover,$$HS{$name}->{flanksize}):$adjacentCover;
		my $antiSenseBG_perNT    = divZero($antiSenseBG,($$HS{$name}->{RHSsize} + $$HS{$name}->{LHSsize}));
		my $antiSensePCBG_perNT  = divZero($antiSensePCBG,($$HS{$name}->{RHSPCsize} + $$HS{$name}->{LHSPCsize}));
		my $leftSignal_perNT  	 = divZero($$HS{$name}->{LHS}->{signal},$$HS{$name}->{LHSsize});
		my $rightSignal_perNT 	 = divZero($$HS{$name}->{RHS}->{signal},$$HS{$name}->{RHSsize});
		
		my $cover_HS          	 = $cover_perNT*$$HS{$name}->{size};
		my $adjacentCover_HS  	 = $adjacentCover_perNT*$$HS{$name}->{size};
		my $antiSenseBG_HS 	 = $antiSenseBG_perNT*$$HS{$name}->{size};
		my $antiSensePCBG_HS 	 = $antiSensePCBG_perNT*$$HS{$name}->{size};
		
		my $strengthRawRawNoBG 	 = $cover-$antiSenseBG;	
		my $strengthRaw 	 = $cover_HS;
		my $strengthNormAdj 	 = round($cover_HS-$adjacentCover_HS);
		my $strengthNormAS 	 = round($cover_HS-$antiSenseBG_HS);
		my $strengthPCNormAS 	 = round($cover_HS-$antiSensePCBG_HS);
		my $strengthNormInput 	 = $inputFile?round($cover_HS-$inputCover_HS):"0";
		
		my ($cs,$from,$to) = ($$HS{$name}->{cs},$$HS{$name}->{from},$$HS{$name}->{to});
		
		#print TMPFINAL join("\t",$cs,$from,$to,$cover,$strengthRawRawNoBG,$strengthRaw,$strengthNormAdj,$strengthNormAS,$strengthPCNormAS,$strengthNormInput,round($adjacentCover_HS),round($antiSenseBG_HS),round($antiSensePCBG_HS),round($inputCover_HS),round($adjacentNCISCover_HS),$leftSignal,$rightSignal,$leftNoise,$rightNoise,$fragCnt)."\n";
		#print STDERR join("\t",$cs,$from,$to,$strengthRaw,$strengthNormAdj,$strengthNormAS,$strengthPCNormAS,round($adjacentCover_HS),round($antiSenseBG_HS),round($antiSensePCBG_HS),$leftSignal,$rightSignal,$leftNoise,$rightNoise,$fragCnt)."\n" if ($verbose);
		print TMPFINAL join("\t",$cs,$from,$to,($cover-round($antiSensePCBG)),$cover,round($antiSensePCBG),$leftFwd,$rightRev,$leftRev,$rightFwd."\n");
		$totStrength += ($cover-round($antiSensePCBG));
	}
	
	close TMPFINAL;
	
	my ($outTmpA,$oTA) = genTempFile();
	my ($outTmpB,$oTB) = genTempFile();
	my ($outTmpC,$oTC) = genTempFile();
	my ($outTmpD,$oTD) = genTempFile();

	## Generate hotspot output Files
	# make initial sorted output file
	sysAndPrint("sort -k1,1 -k2n,2n $prelimOut >$outTmpA",1);

	# Generate bedgraph file with ranks
	my $rank;
	open my $RPIPE, '-|', "sort -k4rn,4rn $outTmpA";
	open RANKTMP, '>', $outTmpB;
	
	while (<$RPIPE>){
		chomp;
		my @X = split(/\t/,$_);
		print RANKTMP join("\t",$X[0],$X[1],$X[2],++$rank)."\n";
	}
	
	close RANKTMP; 
	
	sysAndPrint("sort -k1,1 -k2n,2n $outTmpB >$outTmpC",1);
	sysAndPrint("paste $outTmpA $outTmpC >$outTmpD",1);

	sysAndPrint("cut -f1-4 $outTmpA >$outBG");
#	sysAndPrint("cut -f1-4,8-100 $outTmpD >$outTab");

	open TABIN,  $outTmpD; 
	open TABOUT, '>', $outTab;
	print TABOUT join("\t","cs","from","to","strength","strengthPC","strengthRank","signalReads","bgReads","leftFwd","rightRev","leftRev","rightFwd"."\n");
	while (<TABIN>){
		chomp; 
		my @F = split(/\t/,$_);
		my $pc = sprintf("%6.5f",$F[3]/$totStrength*100);
		print TABOUT  join("\t",@F[0..3],$pc,$F[13],@F[4..9])."\n";
	}
	close TABIN;


}

################################################################################
sub divZero{
	my ($dzNum,$dzDivisor) = @_;

	if ($dzDivisor && $dzDivisor >0){
		return ($dzNum/$dzDivisor);
	}
	return $dzNum;
}
	

################################################################################
sub getNCISfactor{
	my $inNCIS = shift;
	
	my $NCISdata = `cat $inNCIS`;
	my @fNCIS = split(/\s+/,$NCISdata);
	
	return $fNCIS[0];
}
################################################################################
sub printLocus{
	my ($plCS,$plFrom,$plTo,$plName) = @_;
	
	$plFrom = 1 if($plFrom < 1);
	return "" if ($plCS !~ /^chr([0-9]|X|Y)+$/);
	return "" if ($plFrom > $plTo);
	return (join("\t",$plCS,$plFrom,$plTo,$plName)."\n");
}
################################################################################
sub genFR{
	
	my ($fAll,$fF,$fR) = @_;
	open TMPF, ">", $fF;
	open TMPR, ">", $fR;
	
	open IN, $fAll;
	
	while (<IN>){
		print TMPF if ($_ =~ /\s\+\s*$/);
		print TMPR if ($_ =~ /\s\-\s*$/);
	}
	
	close TMPF; close TMPR;
}

################################################################################
sub lineCount{
	my $f = shift;
	my $fSz = `wc -l $f`;
	chomp $fSz;
	$fSz =~ s/^(\d+).+$/$1/;
	return $fSz;
}	

################################################################################
sub median { $_[0]->[ @{$_[0]} / 2 ] }

################################################################################
sub showHELP{

print <<HELP

normalizeStrengthByAdjacentRegions.pl
K. Brick: July 16th 2014
------------------------------------------------------------------------------
ARGS:
'bed=s' => \(my \$inbed),
'in=s'  => \(my \$fragFile),
'v+'    => \(my \$verbose),
'rc+'   => \(my \$recenter),
'nt+'   => \(my \$noTmp),
'h+'    => \(my \$help),
'help+' => \(my \$helpMe),
'out=s' => \(my \$out));

WHAT IT DOES:
This script takes a set of ssDNA defined peaks, a ssDNA fragment bed file and
does the following:
1. Recenters the peaks by the median of the F/R dists
2. Calculate the in-peak background in a number of ways.
3. Output recentered peaks with strength. Strength is normalized by the number
   of wrong-direction fragments (i.e. REV to left, FWD to right) and expressed
   as RPKM, using the total normalized in-hotspot tag count as a denominator.
   
HELP
}
