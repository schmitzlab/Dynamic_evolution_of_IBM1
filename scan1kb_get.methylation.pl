#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use Math::CDF qw(:all);
#use Statistics::Multtest qw(:all);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
my $day=`date +"%D"`;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$gff,$fOut,$tsv);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"fa:s"=>\$fa,
				"gff:s"=>\$gff,
				"tsv:s"=>\$tsv,
				) or &USAGE;
&USAGE unless ($fa and $gff and $fOut and $tsv);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   05/31/23
Description:	This script scans for mCHG-enriched introns longer than 1kb and generates fake isoform CDS sequences to create a new reference that includes these isoform transcripts. It enables the analysis of isoform ratio changes in genes with and without methylation changes to assess how methylation reduction impacts transcription. 
  Options:
  -fa <file>  input file,fa,forced  
  -gff <file>  input file,gff,forced 
  -tsv <file>  input dir,tsv,forced 
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
my %context=(
"CGA"=>"CG","CGC"=>"CG","CGG"=>"CG","CGT"=>"CG","CG"=>"CG","CGN"=>"CG","CHG"=>"CHG","CHH"=>"CHH",
"CAG"=>"CHG","CCG"=>"CHG","CTG"=>"CHG","CAA"=>"CHH","CAC"=>"CHH","CAT"=>"CHH",
"CCA"=>"CHH","CCC"=>"CHH","CCT"=>"CHH","CTA"=>"CHH","CTC"=>"CHH","CTT"=>"CHH",
);
my @code=("CG","CHG","CHH");

######################################1.reading gff and protein_primaryTranscriptOnly.fa files to get the cds location of primary transcripts
$/="\>";
open (IN, $fa) or die $!;
my %primaryTrans;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\n/,$_,2;
	my $name=(split/\s+/,$lines[0])[0];
	$primaryTrans{$name}++;#
}
close IN;
#######################################2. find the location of >1kb introns and create fake isoform transcripts based on that.
my %cdsloc;###rearrange the cds order based on strand state
#my %cdsend;
#my %cdsstart;
my %bcdsloc;###recording cds in original order of gff files
my %ordergene;
my %gene;
$/="\n";
my %strand;
open (IN, $gff) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
	my @lines=split/\s+/,$_;
	#print $_;die;
	if (($lines[2]=~/CDS/||$lines[2]=~/UTR/)&&$lines[0]=~/\d+/) {
		my $gene;
        
		if ($lines[8]=~/ID\=(.*?)\.Araport11.447/) {
			$gene=$1;
            # if($gene=~/AT1G11270/){
            #     print $gene;die;
            # }
			if ($primaryTrans{$gene}) {
				if (!$gene{$gene}) {
					push @{$ordergene{$lines[0]}},$gene;
					$gene{$gene}++;
					$strand{$gene}=$lines[6];
				}

				if(!$bcdsloc{$lines[0]}{$gene}){
					$bcdsloc{$lines[0]}{$gene}="$lines[3] $lines[4]";
				}
				else {
					$bcdsloc{$lines[0]}{$gene}.="\t$lines[3] $lines[4]";
				}
			}
		}
	}
}
close IN;

################# need to include the intron between UTR and cds.but also need to solve the fragment orders.
foreach my $chr (sort {$a<=> $b} keys %bcdsloc) {
	foreach my $gene (@{$ordergene{$chr}}) {
		if (!$bcdsloc{$chr}{$gene}) {
			print "$chr\t$gene\n";
		}
		my @a = split/\t/,$bcdsloc{$chr}{$gene};
		my $temp;
		for(my $i=0;$i<@a-1;$i++){				
			for (my $j=$i+1;$j<@a;$j++){
				#print "$i\t$j\n";
				my ($s1,$e1)=split/\s+/,$a[$i];
				my ($s2,$e2)=split/\s+/,$a[$j];
				if($s1 > $s2){
					$temp=$a[$i];
					$a[$i]=$a[$j];
					$a[$j]=$temp;
					# print "$i\t$j\t$s1\t$s2\n";
					# print Dumper @a;
					# print "############\n";
				}
			}
		}
		@{$cdsloc{$chr}{$gene}}=@a;
	}
}

#################################find the location of intron region and intron length
######record the location of long introns for each genes of arabidopsis
open (OUT, ">$fOut/intron.test.text") or die $!;
my %lintron;
my %filter;
my @orintron;#######record the gene order of intron 
foreach my $chr (sort {$a<=> $b}keys %cdsloc) {
	foreach my $gene (@{$ordergene{$chr}}) {
		for(my $i=1;$i<@{$cdsloc{$chr}{$gene}};$i++){
			my($s1,$e1)=split/\s+/,$cdsloc{$chr}{$gene}[$i-1];##location of start and end point for previous cds
			my($s2,$e2)=split/\s+/,$cdsloc{$chr}{$gene}[$i]; ##location of start and end point for the current cds
			my $end=$s2;
			my $start=$e1;
			my $len=$end-$start;
			if ($len > 1000){
				if (!$filter{$gene}) {
					push @orintron,$gene;
					$filter{$gene}++;
				}
				push @{$lintron{$chr}{$gene}}, "$start $end";
				print OUT "$chr\t$gene\t$start\t$end\n";
			}
		}
	}
}
close OUT;

###################################read methylation and calculated methylation level on >1kb intron region.
%bcdsloc=();
my $name=basename($tsv);
my %record=();
my %pos=();
open (IN, $tsv) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/||$.==1);
	my @lines=split/\s+/,$_;
	$lines[0]=~s/Chr//;
	if ($lines[5]>=3) {
		$record{$lines[0]}{$lines[1]}=$_;
		push @{$pos{$lines[0]}},$lines[1];
	}
}
close IN;

################################################3. evaluation the methylation level of each introns 
my %Mstat=();my %Tstat=();
my %MAstat=();my %TAstat=();
my %Mreads=();my %Treads=();
my %MAreads=();my %TAreads=();
foreach my $chr (sort {$a<=> $b}keys %lintron) {
	my $shift=0;
	foreach my $gene (@orintron) {
		MM:foreach my $intron (@{$lintron{$chr}{$gene}}) {
			my ($start,$end)=split/\s+/,$intron;
			if ($pos{$chr}&&$shift<@{$pos{$chr}}) {
				for (my $i=$shift;$i<@{$pos{$chr}};$i++) {
					if ($pos{$chr}[$i]>=$start&&$pos{$chr}[$i]<=$end) {
						my @info=split/\s+/,$record{$chr}{$pos{$chr}[$i]};
						if ($context{$info[3]}) {
							$Tstat{$gene}{$intron}{$context{$info[3]}}++;
							$TAstat{$context{$info[3]}}++;
							$Mreads{$gene}{$intron}{$context{$info[3]}}+=$info[4];
							$Treads{$gene}{$intron}{$context{$info[3]}}+=$info[5];
							$MAreads{$context{$info[3]}}+=$info[4];
							$TAreads{$context{$info[3]}}+=$info[5];
							if ($info[-1]==1) {
								$Mstat{$gene}{$intron}{$context{$info[3]}}++;
								$MAstat{$context{$info[3]}}++;
							}
						}
						$shift=$i+1;
					}
					elsif ($pos{$chr}[$i]>$end) {
						$shift=$i+1;
						next MM;
					}
					elsif ($pos{$chr}[$i]<$start) {
						$shift=$i+1;
					}
				}
			}
		}
	}
}
open (OUT, ">$fOut/$name.1kbintron.raw.out") or die $!;
print OUT "TOTAL count";
foreach my $cont (@code) {
	if (!$MAstat{$cont}) {
		$MAstat{$cont}="0";
	}
	if (!$TAstat{$cont}) {
		$TAstat{$cont}="0";
	}
	print OUT "\t$cont\t$MAstat{$cont}\t$TAstat{$cont}\t$MAreads{$cont}\t$TAreads{$cont}";
}
print OUT "\n";
foreach my $chr (sort {$a<=> $b}keys %lintron) {
	foreach my $gene (@orintron) {
		if ($Tstat{$gene}) {	
			foreach my $intron (@{$lintron{$chr}{$gene}}) {
				print OUT "$chr\t$gene\t$strand{$gene}\t$intron";
				foreach my $cont (@code) {
					if (!$Mstat{$gene}{$intron}{$cont}) {
						$Mstat{$gene}{$intron}{$cont}="0";
					}
					if (!$Tstat{$gene}{$intron}{$cont}) {
						$Tstat{$gene}{$intron}{$cont}="0";
					}
					if (!$Mreads{$gene}{$intron}{$cont}) {
						$Mreads{$gene}{$intron}{$cont}="0";
					}
					if (!$Treads{$gene}{$intron}{$cont}) {
						$Treads{$gene}{$intron}{$cont}="0";
					}
					print OUT "\t$cont\t$Mstat{$gene}{$intron}{$cont}\t$Tstat{$gene}{$intron}{$cont}\t$Mreads{$gene}{$intron}{$cont}\t$Treads{$gene}{$intron}{$cont}";
				}
				print OUT "\n";
			}
		}
	}
}

#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#���б��е����ֵ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#���б��е���Сֵ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


