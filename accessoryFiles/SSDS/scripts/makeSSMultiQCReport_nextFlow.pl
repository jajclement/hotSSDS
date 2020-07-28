use strict; 
use Time::HiRes qw/ time sleep /;
use Getopt::Long;

GetOptions ('g=s' => \(my $genomeName));
GetOptions ('h=s' => \(my $hotspots));

my ($initialBAM,@bedFiles) = @ARGV;

my $ssdsFile; 
for my $currentBED(@bedFiles){
	$ssdsFile = $currentBED if ($currentBED =~ /ssDNA_type1.bed/);
}

my $ssdsStem = $ssdsFile;
$ssdsStem =~ s/^(.+)\.\S+\.(ssDNA_type1.bed|ssDNA_type1.bam|ssDNA_type1.final.bam)/$1/;

my $outFile = $ssdsStem.'.v170814.SSDSreport.tab';

my $t1Bed   = $ssdsStem.'*ssDNA_type1.bed';
my $t2Bed   = $ssdsStem.'*ssDNA_type2.bed';
my $dsBed   = $ssdsStem.'*dsDNA.bed';
my $unBed   = $ssdsStem.'*unclassified.bed';
my $initBAM = $ssdsStem.'.*.bam';

die if ($outFile  eq $ssdsFile);
die if ($ssdsStem eq $ssdsFile);

## Sort all bed files correctly
my $randStem  = 'ssds_stat_'.int(rand()*100000000);
my $t1Tmp     = $randStem.'_t1.sorted.bed';
my $t2Tmp     = $randStem.'_t2.sorted.bed';
my $dsTmp     = $randStem.'_ds.sorted.bed';
my $unTmp     = $randStem.'_un.sorted.bed';
my $tmpOut    = $randStem.'_output.tab';

system('sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 '.$t1Bed.' >'.$t1Tmp) if ($t1Bed);
system('sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 '.$t2Bed.' >'.$t2Tmp) if ($t2Bed);
system('sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 '.$dsBed.' >'.$dsTmp) if ($dsBed);
system('sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 '.$unBed.' >'.$unTmp) if ($unBed);

## open OUTFILE
open OUT, '>', $tmpOut;

## Get species
my $species = $genomeName?$genomeName:getSpecies($ssdsFile);
print $species."\n";

## generate overall Stats
my $t1Sz = printBEDcounts($t1Tmp,$t1Bed,'ssDNA_type1_fragments');
my $t2Sz = printBEDcounts($t2Tmp,$t2Bed,'ssDNA_type2_fragments');
my $dsSz = printBEDcounts($dsTmp,$dsBed,'dsDNA_fragments');
my $unSz = printBEDcounts($unTmp,$unBed,'unclassified_fragments');
printBAMcounts($initialBAM,'total_fragments');

getStats($t1Tmp,$t1Bed,'ssDNA_type1');
getStats($t2Tmp,$t2Bed,'ssDNA_type2');
getStats($dsTmp,$dsBed,'dsDNA');
getStats($unTmp,$unBed,'unclassified');

## Get FRIPs
getFRIPs($t1Tmp,$t1Bed,$t1Sz,$species,'ssType1');
getFRIPs($t2Tmp,$t2Bed,$t2Sz,$species,'ssType2');
getFRIPs($dsTmp,$dsBed,$dsSz,$species,'dsDNA');
getFRIPs($unTmp,$unBed,$unSz,$species,'unclassified');

close OUT; 

system("mv $tmpOut $outFile");
system('rm '.$randStem.'*');

############################################################
sub getFile{
	my ($getFolder,$getMatch) = @_;
	
	my @f = `ls $getFolder\/$getMatch`;
	
	chomp @f; 
	
	if ($#f == 0 && -e ($f[0])){
		return $f[0];
	}else{
		return '';
	}
}

############################################################
sub printBEDcounts{
	my ($bed_file,$bedCheck,$title) = @_;
	
	## if no bed exists, print filler so multiQC doesn't crash
	if (not ($bedCheck)){
		print OUT join("\t","totinfo",$title,"0\n");
		print OUT join("\t","totinfo","unique_".$title,"0\n");
		return 0;
	}
	
	my $cnt_line        = `wc -l $bed_file`;
	my $unique_cnt_line = `uniq $bed_file |wc -l`;
	
	$cnt_line        =~ s/^(\d+).*/$1/; chomp $cnt_line;
	$unique_cnt_line =~ s/^(\d+).*/$1/; chomp $unique_cnt_line;
	
	print OUT join("\t","totinfo",$title,$cnt_line."\n");
	print OUT join("\t","totinfo","unique_".$title,$unique_cnt_line."\n");
	
	return $cnt_line;
}

############################################################
sub printBAMcounts{
	my ($bam_file,$title) = @_;
	my $cnt_line = countBAM($bam_file);
	
	## modified this to only count unaligned adapter !! Right or wrong ??
	my $seqsFile = 'ssds_sequences_for_adapterCheck.txt';
	my $adapter  = `samtools view $bam_file |cut -f10 >$seqsFile`;
	my $adapter  = getAdapterReads($seqsFile);

	chomp $adapter; 
	
	$adapter    = int($adapter/2); ## because its PE
	
	print OUT join("\t","totinfo",$title,$cnt_line."\n");
	print OUT join("\t","totinfo",'adapter',$adapter."\n");
}

############################################################
sub countBAM{
	my $inBAM = shift;
	
	my ($bamCount);
	
	if (-e $inBAM.'.bai'){
		open my $PIPE, '-|', 'samtools idxstats '.$inBAM;
		while (<$PIPE>){
			chomp;
			my @F = split(/\t/,$_);
			$bamCount += ($F[2]+$F[3]);
		}
	}else{
		$bamCount = `samtools view -c $inBAM`;
		$bamCount    =~ s/^(\d+).*/$1/;
	}
		
	$bamCount   /= 2; ## because its PE
	
	return $bamCount;
}

############################################################
sub getStats{
	my ($bed_file,$bedCheck,$title) = @_;
	
	if (not ($bedCheck)){
		print OUT join("\t","ITR_".$title,"0","0\n");
		print OUT join("\t","uH_".$title,"0","0\n");
		print OUT join("\t","Offset_".$title,"0","0\n");
		print OUT join("\t","Fragments_".$title,"0","0\n");
		return ;
	}
	
	my (%itr,%uH,%fillin,%fragLen);
	my $cmd = 'uniq '.$bed_file.' |head -n 1000000 ';
	open my $PIPE, '-|', $cmd; 
	
	while (<$PIPE>){
		chomp;
		my @F = split(/\t/,$_);
		my @X = split(/_/,$F[4]); 
		$itr{$X[0]}++;
		$uH{$X[1]}++;
		$fillin{$X[2]}++;
		$fragLen{$F[2]-$F[1]}++ if (($F[2]-$F[1]) < 501);
	}

	toHist(\%itr    ,"ITR_".$title);
	toHist(\%uH     ,"uH_".$title);
	toHist(\%fillin ,"Offset_".$title);
	toHist(\%fragLen,"Fragments_".$title);
}

############################################################
sub toHist{
	
	my ($refHash,$type) = @_;

	my $min = 9e999;
	my $max = -9e999;
	
	for my $k(sort {$a <=> $b} keys(%{$refHash})){
		print OUT join("\t",$type,$k,$$refHash{$k}."\n");
	}
}

############################################################
sub getSpecies{
	my $inName = shift; 
	
	$inName =~ s/^.+\.(\S+?)\.ssDNA_type1.+/$1/;
	#$inName =~ s/(Bsub168|E_Coli_K12|E_coli_K12|MonDom5|bSub168|bsub|bwa_hg19|canfam3|eColiK12|ecoli|hg19|hg38|mm10|mm10TC1|rn5|saccer3)//;
	
	return lc($inName);
}
	
############################################################
sub getFRIPs{
	
	my ($inBed,$bedCheck,$totalFrags,$species,$lbl) = @_;

	my %hs; 
	
	#open my $PIPE, '-|', 'find '.$ENV{'HOTSPOTS'}.' -name "*bed"';
	open my $PIPE, '-|', 'find '.$hotspots.' -name "*bed"';	

	while (<$PIPE>){
		chomp; 
		$_ =~ /^.+\/(\S+)\/(\S+)\.bed/;
		my ($genome,$name) = ($1,$2);
		$hs{$genome}->{$name} = $_;
	}
	close $PIPE; 
	
	for my $sp (sort keys(%hs)){
		for my $hs (sort keys(%{$hs{$sp}})){	
			if (not ($bedCheck)){
				## Default filler if no file exists (for Unc)
				print OUT join("\t","FRIP\_$lbl",$sp.':'.$hs,"0")."\n";
			}else{
				if ($sp eq $species){
					#my $hsFile = $hsFolder{$sp}.'/'.$hs{$sp}->{$hs};
					my $hsFile = $hs{$sp}->{$hs};
					my $inHS = `intersectBed -a $inBed -b $hsFile -sorted |wc -l `;
					$inHS =~ s/^(\d+)\s+.+/$1/;
					my $frip = $totalFrags?$inHS/$totalFrags*100:0;
					print OUT join("\t","FRIP\_$lbl",$sp.':'.$hs,$frip)."\n";
				}else{
					print OUT join("\t","FRIP\_$lbl",$sp.':'.$hs,"0")."\n";
				}
			}
		}
	}
}
	
############################################################
sub getAdapterReads{

	my $inSeqFile = shift;

	#revComp of universal TruSeq LT adapter 
	#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
	my $univ = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT';

	#TruSeq multiplex adapter
	#GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	my $truseq = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC';

	my $nMM = 2;

	my @seqs;
	addSeqs(\@seqs,$univ,5);
	addSeqs(\@seqs,$truseq,5);

	my $lineNum;
	my $adapterCount = 0;

	open INSEQS, $inSeqFile;
	## Fastest search first
	while (<INSEQS>){
		$lineNum++;
		chomp; 
		my $mySeq = $_;
		if ($mySeq =~ /(GATCGGAAGAGCACACGTC|CACGTCTGAACTCCAGTCAC|GATCGGAAGAGCGT|AAGAGTGTAGATCTCGGTG)/){
			$adapterCount++;
			next;
		}

		## OK Is it the universal adapter ?
		my $uSeq = substr($mySeq,1,30);
		if ($univ =~ /($uSeq)/){
			$adapterCount++;
			next;
		}	

		## OK Is it the truseq adapter ?
		my $tSeq = substr($mySeq,1,30);
		if ($univ =~ /($tSeq)/){
			$adapterCount++;
			next;
		}

		## OK .... the slow option ... allow $nM mismatches to each
		my $mySeq2 = $mySeq;
		if ($mySeq2 =~ s/(AGATC|GGAAG|AGCGT|CGTGT|AGGGA|AAGAG|TGTAG|ATCTC|GATCG|GAAGA|GCACA|CGTCT|GAACT|CCAGT)/NNNNN/){
			if ($mySeq2 =~ s/(AGATC|GGAAG|AGCGT|CGTGT|AGGGA|AAGAG|TGTAG|ATCTC|GATCG|GAAGA|GCACA|CGTCT|GAACT|CCAGT)/NNNNN/){
				if ($mySeq2 =~ s/(AGATC|GGAAG|AGCGT|CGTGT|AGGGA|AAGAG|TGTAG|ATCTC|GATCG|GAAGA|GCACA|CGTCT|GAACT|CCAGT)/NNNNN/){
					for my $i(@seqs){
						$adapterCount++;
					}
				}
			}
		}

	}
	return $adapterCount;
}

########################################################################
sub addSeqs{
	my ($arr,$s,$n) = @_;
	while ($s =~ s/^(.{$n})(.{15})/$2/){
		push @{$arr}, $1.$2;
	}
}
