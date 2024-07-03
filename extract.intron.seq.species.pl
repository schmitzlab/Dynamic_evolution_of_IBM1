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
my ($fIn,$fIn2,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
                "i2:s"=>\$fIn2,
				) or &USAGE;
&USAGE unless ($fIn and $fIn2 and $fOut);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2022/09/15
Description:	find the intron region and extract the intronic sequence for 34 species.
Usage:
  Options:
  -i <list>  input file,for location of intron intron.methy.info.csv  
  -i2 <dir>  input dir,for getting gff file and genome reference file,forced   
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
###the overall goals of this scripts:
#####################1. read gene and species name
#####################2. obtain gene absolute location from gff file
#####################3. extract the whole gene sequences from genome file
#####################4. obtain intron sequence based on relative location in intron.methy.info.csv and extract intron sequence
############################1. read gene and species name
$/="\n";
open (IN, $fIn) or die $!;
my %inloc;###relative location of each intron sequence  on genes
while (<IN>) {
	chomp;
    next if (/^$/||$.==1);
	my @lines=split/\,/,$_;
	my @label=split/\_/,$lines[0];
	my $sp="$label[-2]\_$label[-1]";
	$lines[0]=~s/\_$sp//;
	push @{$inloc{$sp}{$lines[0]}},"$lines[1] $lines[2]";
	#print $lines[0];die;
}
close IN;
#print Dumper %inloc;
# ###read the gff files
my @gffs=glob"$fIn2/*/*gff*";###get the gff files
my @fa=glob"$fIn2/*/*add*.fa";
my %seq;
foreach my $f (@gffs){
	
	my $sp=(split/\//,$f)[-2];
	#print $sp,"\n";
	my %geneloc;
	my %chr;
	$/="\n";
	my %filter;
	open (IN, $f) or die $!;
	while (<IN>) {
		chomp;
		next if(/^#/);
		my @lines=split/\t/,$_;
		#print Dumper %{$inloc{$sp}};die;
		foreach my $id (keys %{$inloc{$sp}}){
			if($lines[8]=~/$id/ && $lines[2]=~/mRNA/){
				if(!$filter{$id}){
					$filter{$id}++;
					$geneloc{$id}="$lines[0] $lines[3] $lines[4] $lines[6]";
					$chr{$lines[0]}++;
				}
				#print $id,"\n";
				#print $lines[0],"\n";
			}
		}
	}
	close IN;
	#####record sequence of target chr/scaffold
	my $genome;
	foreach my $fa (@fa){
		if($fa=~/$sp/){
			$genome=$fa;
			last;
		}
	}
	#print $genome,"\n";
	$/="\>";
	open (IN, $genome) or die $!;
	while (<IN>) {
		chomp;
		next if(/^$/);
		my @lines=split/\n/,$_,2;
		my $scafid=(split/\s+/,$lines[0])[0];
		#print "$scafid\n";
		if($chr{$scafid}){
			#print $scafid;die;
			$lines[1]=~s/\n//g;
			foreach my $id (keys %geneloc){
				my ($chr,$s,$e,$strand)=split/\s+/,$geneloc{$id};
				if($chr eq $scafid){
					my $seq=substr($lines[1],$s-1,$e-$s+1);
					if($strand eq "-"){
						$seq=&revcom($seq);
					}
					$seq{$sp}{$id}=$seq;#######record the whole gene sequence for each gene of each species
					#print "$sp\t$id\n$seq\n";
				}
			}
		}
	}
	close IN;
	$/="\n";
}
################################################extract intron seq
open (OUT, ">$fOut") or die $!;
my $count=1;
foreach my $sp (sort keys %inloc){
	foreach my $id (keys %{$inloc{$sp}}){
		###########output intron sequence
		foreach my $loc (@{$inloc{$sp}{$id}}){########intron location
			#print "$sp\t$id\n$seq{$sp}{$id}";die;
			my ($s,$e)=split/\s+/,$loc;
			my $seq;
			if(!$seq{$sp}{$id}){
				print "$sp\t$id";
			}
			else{
				$seq=substr($seq{$sp}{$id},$s-1,$e-$s+1);
			}
			
			print OUT ">$sp\_$count $id $sp $s $e\n";
			print OUT "$seq\n";	
			$count++;		
		}
	}
	
}
close OUT;
# #######################################################################################
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


