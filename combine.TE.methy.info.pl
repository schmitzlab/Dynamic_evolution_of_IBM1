#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
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
my ($fIn,$fOut,$fIn2,$fIn3);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"i2:s"=>\$fIn2,
				"i3:s"=>\$fIn3,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $fIn2 and $fIn3);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2022/09/21
Description:	this program is used to convert dmr_original list to standard gabit trait form
Usage:
  Options:
  -i <file>  input file,RMRB.embl,forced  
  -i2 <file>  input file,34.intron.seq.fa,forced 
  -i3 <file>  input file,34.intron.seq.fa.out,forced 
  -o <file>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
$/="\//\n";
open (IN, $fIn) or die $!;
my %type;
while (<IN>) {
	chomp;
	my @lines=split/\n/,$_;
    my $id;
    my $type="unknown";
    my $subtype="unknown";
    #print $_;die;
    foreach my $l (@lines){
        #print $l,"\n";
        #print "###################\n";
        if($l=~/ID\s+(.*?)\;/){
            $id=$1;
        }
        elsif($l=~/\sType\:\s+(.*)/){
            $type=$1;
        }
        elsif($l=~/\sSubType:\s+(.*)/){
			if($1){
				$subtype=$1;
			}
            
        }
    }
   $type{$id}="$type\t$subtype";
   #print "$id\t$type\t$subtype\n";
}
close IN;
#print Dumper %type;die;
###########################################extract id convector from intron.fa
$/="\>";
my %convert;
open (IN, $fIn2) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\n/,$_,2;
	my @head=split/\s+/,$lines[0],2;
	$convert{$head[0]}=$head[1];
}
close IN;
#print Dumper %convert;die;
##########################################extract TE location from repeatmasker output file and give the type annotation
open (OUT, ">$fOut") or die $!;
$/="\n";
my %table;
print OUT "Geneid_species\tinstart\tinend\ttestart\tteend\tmotif\ttype\tsubtype\n";
open (IN, $fIn3) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	if($_=~/\s+(\w+\_\w+\_\d+)\s+(\d+)\s+(\d+)\s+.*?\s+.*?\s+(.*?)\s+(.*?)\s+/){
		my @a=($1,$2,$3,$4,$5);
		if($convert{$a[0]}){
			my @info=split/\s+/,$convert{$a[0]};
			my $comb="$info[0]\_$info[1]";
			print OUT "$comb\t$info[2]\t$info[3]\t$a[1]\t$a[2]\t$a[3]";
			if($a[4]=~/Unspecified/){
				if($type{$a[3]}){
					print OUT "\t$type{$a[3]}\n";
				}
				else {
					print OUT "\tnotfound\tnotfound\n";
				}
			}
			else {
				print OUT "\tSimple_repeat\tSimple_repeat\n";
			}
		}
	}
}
close IN;
close OUT;
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


