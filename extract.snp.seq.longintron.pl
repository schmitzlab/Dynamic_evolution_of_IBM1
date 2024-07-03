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
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fOut,$fIn,$fIn2,$fIn3);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"g:s"=>\$fIn2,
				"m:s"=>\$fIn3,
				) or &USAGE;
&USAGE unless ($fOut and $fIn and $fIn2 and $fIn3);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to get the snp
  Options:
  -i <file>  input file,snp files, forced 
  -o <dir>  output file, forced  
  -g <file>  output file,all long intron file, forced  
  -m <file>  output file,methylated long intron file, forced
  -h         Help

USAGE
	print $usage;
	exit;
}

#######################################read gff files
my @introns;#record all long intron location
my %length;###########length of sequence investigate
my %nsites;###############snp number on interested regions
my $tlen;##total gene length;
my $length;

# read the location of introns
open (IN, $fIn2) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	#my @lines=split/\s+/,$_;
	push @introns,$_;
}
close IN;

# read the location of methylated introns
my %mintrons; # record methylated introns
open (IN, $fIn3) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\s+/,$_;
	$mintrons{$lines[1]}{$lines[3]}++;
}
close IN;

#print $length;
#$length{'long_intron'}=$tlen-$length;
#$length{'methyl_intron'}=$target;
############################################extract SNP
my $shift=0;
my %index;
my %snp;##########snp in the introns
#my %status;##########record methylation/unmethylation status of each introns
#open (OUT, ">$fOut") or die $!;
open (IN, $fIn) or die $!;##snp list
MM:while (<IN>) {
	chomp;
	next if (/^$/);
    if($.==1){
        my @head=split/\s+/,$_;
		for (my $i=2;$i<@head;$i++){
			$index{$i}=$head[$i];
		}
		#print OUT "$_\n";
        next;
    }
	my @lines=split/\s+/,$_;
	#############################
	for (my $i=$shift;$i<@introns;$i++){
		my ($chr,$gene,$start,$end)=split/\s+/,$introns[$i];

		if($lines[0]==$chr&&($lines[1]>=$start&&$lines[1]<=$end)){
			#print $introns[$i];die;
			for (my $j=2;$j<@lines;$j++) {
				push @{$snp{$introns[$i]}{$j}},$lines[$j];
			}
			$shift=$i;
			next MM;
		}
	}
}
close IN;

#output the snp fasta files
mkdir $fOut if (! -d $fOut);
foreach my $intron (keys %snp){
	my ($chr,$gene,$start,$end)=split/\s+/,$intron;

	if($mintrons{$gene}{$start}){
		open (OUT, ">$fOut/methylated.$gene.$start.snp.list") or die $!;
		foreach my $sam (sort keys %{$snp{$intron}}){
			my $seq=join("",@{$snp{$intron}{$sam}});
			print OUT ">$sam\n";
			print OUT "$seq\n";
		}
		close OUT;
	}
	else{
		open (OUT, ">$fOut/unmethylated.$gene.$start.snp.list") or die $!;
		foreach my $sam (sort keys %{$snp{$intron}}){
			my $seq=join("",@{$snp{$intron}{$sam}});
			print OUT ">$sam\n";
			print OUT "$seq\n";
		}
		close OUT;
	}
	$nsites{$intron}=$#{$snp{$intron}{'88'}}+1;
}

#die;

# output the statistics
open (OUT, ">$fOut/intron.length.snpsites.txt") or die $!;
print OUT "geneid\tchr\tstart\tend\tstatus\tseq_len\tnsites\n";

foreach my $in (@introns){

	my ($chr,$gene,$start,$end)=split/\s+/,$in;
	my $length=$end-$start+1;
	my $status;
	if($mintrons{$gene}{$start}){

		$status="methy";
	}
	else {$status="unmethy"}
	if($nsites{$in}){
		print OUT "$gene\t$chr\t$start\t$end\t$status\t$length\t$nsites{$in}\n";
	}
}
close OUT;
#print $#{$msnp{'88'}};die;
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


