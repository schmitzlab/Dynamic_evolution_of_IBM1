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
my ($fOut,$fIn,$gff);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"g:s"=>\$gff,
				) or &USAGE;
&USAGE unless ($fOut and $fIn and $gff);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to find the methylation level of each gene
  Options:
  -i <file>  input file,snp files, forced 
  -o <dir>  output file, forced  
  -g <file>  output file,gff file, forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
my @type=("cds","noncoding","methyl_intron");
# ###IBM1 intron methylation start end
my $ms=2430000;
my $me=2430640;
my $target=$me-$ms;
#######################################read gff files
my @cds;
my %length;###########length of sequence investigate
my %nsites;###############snp number on interested regions
my $tlen;##total gene length;
my $length;
open (IN, $gff) or die $!;##snp list
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\s+/,$_;
    if($lines[8] && $lines[8]=~/AT3G07610\.1/){
      if($lines[2]=~/mRNA/){
	  	$tlen=$lines[4]-$lines[3];
	  }
	  if($lines[2]=~/CDS/){
			push @cds,"$lines[0]\t$lines[3]\t$lines[4]";	   	  
	   		$length+=$lines[4]-$lines[3];
	   }
    }
}
close IN;
#print $length;
$length{'noncoding'}=$tlen-$length;
$length{'cds'}=$length;
$length{'methyl_intron'}=$target;
############################################extract SNP
my %index;
my %snp;##########snp in the methylation region of intron
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
	my $count=0;
    foreach my $loc (@cds){
		my ($chr,$s,$e)=split/\s+/,$loc;
		if($lines[1]>=$s&&$lines[1]<=$e){
			$count++;
			for (my $i=2;$i<@lines;$i++) {
				push @{$snp{'cds'}{$i}},$lines[$i];
			}
			next MM;
		}
	}
	if($count<1){
		for (my $i=2;$i<@lines;$i++) {
		  	push @{$snp{'noncoding'}{$i}},$lines[$i];
		}
	}
	if($lines[1]>=$ms && $lines[1]<=$me){
	#print $_,"\n";
		for (my $i=2;$i<@lines;$i++) {
			push @{$snp{'methyl_intron'}{$i}},$lines[$i];	
		}
	}
}
close IN;

#print Dumper @type;die;
mkdir $fOut if (! -d $fOut);
foreach my $t (@type){
	open (OUT, ">$fOut/$t.snp.list") or die $!;
	foreach my $sam (sort keys %{$snp{$t}}){
		my $seq=join("",@{$snp{$t}{$sam}});
		print OUT ">$sam\n";
		print OUT "$seq\n";
	}
	close OUT;
	$nsites{$t}=$#{$snp{$t}{'88'}}+1;
}
open (OUT, ">$fOut/target.seq.length.txt") or die $!;
my $nonlen=$tlen-$length;
print OUT "type\tseq_len\tnsites\n";
foreach my $t (@type){
	print OUT "$t\t$length{$t}\t$nsites{$t}\n";
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


