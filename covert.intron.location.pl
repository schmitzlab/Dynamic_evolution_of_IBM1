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
my ($fIn,$fIn2,$fIn3,$fOut,$list);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"l:s"=>\$list,
				"i:s"=>\$fIn,
                "i2:s"=>\$fIn2,
                "i3:s"=>\$fIn3,
				) or &USAGE;
&USAGE unless ($fIn and $fIn2 and $fIn3  and $fOut and $list);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2019/11/22
Description:	this program is used to find the location of intron, exon, and connection between protein domain and exon location
                v2.0, used for gff with target gene in negative string, and in an order from 5' to 3' UTR,solve the problem that the 
                script too rely on exist of UTR. this new version change the way of how id are read to fit with new id format. make some
                input file reading changes to fit the format of id in the tree file in this round of analysis
Usage:
  Options:
  -l <file>  input file,species info,forced 
  -i <file>  input file,genetree list,forced  
  -i2 <dir>  input directory,gff files,forced  
  -i3 <file>  input file,jmjC pep domain analysis tsv file,forced  
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if(!-d $fOut);
###########read species strand information
open (IN, $list) or die $!;
my %direction;
while (<IN>) {
	chomp;
    next if (/^$/||$.==1);
	my @lines=split/\,/,$_;
	$direction{$lines[0]}=$lines[-1];
}
close IN;
#print Dumper %direction;die;
######################extract gff files
my @gffs=glob"$fIn2/*/*gff*";
foreach my $sp (keys %direction){
	my $gff;
	#print $sp;die;
	my $nod=$direction{$sp};####record the order of CDS coordinator in negative strand: 53 (5'->3') or 35(3'->5')
	foreach my $g (@gffs){
		if($g=~/$sp/){
			$gff=$g;#############find the gff file for each species
			last;
		}
	}
	######################read gene id 
	open (IN, $fIn) or die $!;
	my %gene;
	while (<IN>) {
		chomp;
		next if (/^$/);
		while($_=~/(.*?)\:\d\.\d+/g){
			my $term=$1;
			if($term=~/$sp/){
				$term=~s/\,//g;
				$term=~s/\(//g;
				my $id=(split/\_$sp/,$term)[0];
				$gene{$id}++;
				#print $id;die;
			}
		}	
	}
	close IN;
	# ############################2 get the location and length of jmjC domain on each specie's IBM1.
	my %jd;##jmjC domain location
	open (IN, $fIn3) or die $!;
	while (<IN>) {
	chomp;
	my @lines=split/\t/,$_;
	foreach my $id (keys %gene){
		(my $name=$id)=~s/\.\d+$//;
		#print $name;
		if($lines[0]=~/$name/){
			if($lines[3]=~/ProSiteProfiles/ && $lines[5]=~/JmjC/){
				my $start=$lines[6]*3-2;
				my $end=$lines[7]*3-2;
				$jd{$id}="$start\t$end";####convert to cds location'
				#print "$id\t$id\t$start\t$end\n";
				}
			}
		}
	}
	close IN;
	# ############################3 extract cds intron location convert to positive strand direction

	open (IN, $gff) or die $!;
	my %reladis;#mRNA length;
	my %tcds;##total CDS length;the length that remove intron region
	my %begin;##record the start end point
	my %tail;##record the coordinate end point
	my %filter;
	while (<IN>) {
		chomp;
		next if($_=~/^\#/||$_=~/^$/);
		my @lines=split/\t/,$_;
		if($lines[2]=~/mRNA/){
			if($lines[6] eq "+"){
				foreach my $id (keys %gene){
					if($lines[8]=~/$id/){
						#print $_;die;
						if(!$filter{$id}){
							$filter{$id}++;
							$reladis{$id}=$lines[4]-$lines[3];
							$begin{$id}=$lines[3];
							$tail{$id}=$lines[4];
						}
					}
				}
			}
			else{
				foreach my $id (keys %gene){
					if($lines[8]=~/$id/){
						#print "$id\n$_\n";
						if(!$filter{$id}){
							$filter{$id}++;
							$reladis{$id}=$lines[4]-$lines[3];
							$begin{$id}=$lines[3];
							$tail{$id}=$lines[4];
						}
					}
				}
			}
		}
		if($lines[2]=~/CDS/){
			foreach my $id (keys %gene){
				if($lines[8]=~/$id/){
					$tcds{$id}+=$lines[4]-$lines[3]+1;
				}
			}
		}
	}
	close IN;
	
	open (OUT, ">$fOut/$sp.gene.structure.bed") or die $!;
	open (OUT2, ">$fOut/$sp.domain.location.bed") or die $!;
	open (IN, $gff) or die $!;####read gff file###Record location information, convert negative to positive direction
	my ($s,$e);
	my ($es,$ed);####the start and end point of each exon.
	my %count;
	my ($exons,$exone,$dis);

	while (<IN>) {
		chomp;
		next if($_=~/^\#/||$_=~/^$/);
		my @lines=split/\t/,$_;
		if($lines[2]=~/UTR/||$lines[2]=~/CDS/){
			foreach my $id (keys %gene){
				#print "$id\n$lines[8]";die;
				if($lines[8]=~/$id/){
					if(!$count{$id}){
						$count{$id}=1;
					};
					if($lines[6] eq "+"){
						$s=$lines[3]-$begin{$id};
						$e=$lines[4]-$begin{$id};
						if($lines[2]=~/CDS/){
							my $frag=$e-$s;
							if($count{$id}==1){
								$es=0;
								$ed=$es+$frag;
							}
							else{
								$es=$ed+1;
								$ed=$es+$frag;
							}
							my ($start,$end)=split/\t/,$jd{$id};
							
							if($start>=$es&&$start<=$ed){
								$dis=$start-$es;
								$exons=$s+$dis;###the domain's start point based on the exon intron coordinate.
								$exone=$e;
								#print "$start\t$end";die;
								print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
							}
							elsif($es>=$start&&$ed<=$end){
								$exons=$s;
								$exone=$e;
								print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
							}
							elsif($end>=$es&&$end<=$ed){
								$dis=$end-$es;
								$exons=$s;
								$exone=$s+$dis;
								print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
							}
							$count{$id}++;
						}
						if($lines[2]=~/UTR/){
							print OUT "$id\_$sp\t$s\t$e\tUTR\n";
						}
						else{
							print OUT "$id\_$sp\t$s\t$e\tCDS\n";
						}
					}
					else {
						if($nod eq "35"){
							#print "......";
							$e=abs($lines[3]-$begin{$id}-$reladis{$id});
							$s=abs($lines[4]-$begin{$id}-$reladis{$id});
							if($lines[2]=~/CDS/){
								my $frag=$e-$s;
								if($count{$id}==1){
									$ed=$tcds{$id};
									$es=$ed-$frag;
								}
								else{
									$ed=$es-1;
									$es=$ed-$frag;
								}
								
								my ($start,$end)=split/\t/,$jd{$id};
								#print "$start\t$end\t$tcds{$id}\n";die;
								if($start>=$es&&$start<=$ed){
									$dis=$start-$es;
									$exons=$s+$dis;###the domain's start point based on the exon intron coordinate.
									$exone=$e;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								elsif($es>=$start&&$ed<=$end){
									$exons=$s;
									$exone=$e;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								elsif($end>=$es&&$end<=$ed){
									#print "$es\t$ed\n";
									#print "$s\t$e\n";die;
									$dis=$end-$es;
									$exons=$s;
									$exone=$s+$dis;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								$count{$id}++;
							}
							if($lines[2]=~/UTR/){
								print OUT "$id\_$sp\t$s\t$e\tUTR\n";
							}
							else{
								print OUT "$id\_$sp\t$s\t$e\tCDS\n";
							}
						}
						else{
							$s=abs($lines[4]-$tail{$id});
							$e=abs($lines[3]-$tail{$id});
							if($lines[2]=~/CDS/){
								my $frag=$e-$s;
								if($count{$id}==1){
									$es=0;
									$ed=$es+$frag;
								}
								else{
									$es=$ed+1;
									$ed=$es+$frag;
								}
								my ($start,$end)=split/\t/,$jd{$id};
								#print "$start\t$end\t$es\t$ed\n";
								if($start>=$es&&$start<=$ed){
									$dis=$start-$es;
									$exons=$s+$dis;###the domain's start point based on the exon intron coordinate.
									$exone=$e;
									#print "$start\t$end\t$es\t$ed\n";die;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								elsif($es>=$start&&$ed<=$end){
									$exons=$s;
									$exone=$e;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								elsif($end>=$es&&$end<=$ed){
									$dis=$end-$es;
									$exons=$s;
									$exone=$s+$dis;
									print OUT2 "$id\_$sp\t$exons\t$exone\tjmjC-domain\n";
								}
								$count{$id}++;
							}
							if($lines[2]=~/UTR/){
								print OUT "$id\_$sp\t$s\t$e\tUTR\n";
							}
							else{
								#print $_,"\n";
								print OUT "$id\_$sp\t$s\t$e\tCDS\n";
							}
						}
					}
				}
			}
		}
	}
	close OUT;
	close OUT2;
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


