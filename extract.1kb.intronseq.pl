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
my ($fa,$in,$min,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"fa:s"=>\$fa,
				"i:s"=>\$in,
				"m:s"=>\$min,
				) or &USAGE;
&USAGE unless ($fa and $in and $fOut and $min);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   05/31/23
Description:	this script is used to extract the intron sequence based on intron location list
  Options:
  -fa <file>  input file,genome assembly fasta file,forced 
  -i <file>  input file,long intron location,forced 
  -m <file>  input dir,methylated long intron location,forced 
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

###################################read intron location
my %introns;
open (IN, $in) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
    my @lines=split/\s+/,$_;
	push @{$introns{$lines[0]}},$_;
}
close IN;
##################################read methylated intron location
my %mintrons;
open (IN, $min) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
	my @lines=split/\s+/,$_;
    $mintrons{$lines[1]}{$lines[3]}="$lines[3]\t$lines[4]";
}
close IN;
################################ read genome assembly and extract intron sequence
mkdir $fOut if (! -d $fOut);
open (OUT, ">$fOut/long.unmethylated.intron.fa") or die $!;
open (OUT2, ">$fOut/long.methylated.intron.fa") or die $!;

$/="\>";
open (IN, $fa) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\n/,$_,2;
    my @info=split/\s+/,$lines[0];
    if($introns{$info[0]}){
        # extract intron sequence for each record
        foreach my $loc (@{$introns{$info[0]}}){
            my @record = split/\s+/,$loc;
            my $seq=substr($lines[1],$record[2],$record[3]-$record[2]+1);
			$seq=~s/\n//g;
			
			if($mintrons{$record[1]}{$record[2]}){
				print OUT2 ">$loc\n$seq\n";
			}
			else{
				print OUT ">$loc\n$seq\n";
			}
        }

    }
}
close IN;
close OUT;
close OUT2;
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


