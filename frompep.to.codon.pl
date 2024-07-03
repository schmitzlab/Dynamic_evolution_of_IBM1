#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
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
Program Date:   2019/10/23
Description:	this program is used to convert peptide to codon
Usage:
  Options:
  -i  <file>  input file,pep,fa,forced  
  -i2 <file>  input file,cds,fa,forced
  -i3 <file>  input file,pep align,fa,forced  
  -o  <dir>  output dir,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
##################read pep
$/="\>";
my %pep;
open (IN, $fIn) or die $!;
while (<IN>) {
	chomp;
    next if (/^$/);
	my @lines=split/\n/,$_,2;
    $lines[1]=~s/\n//g;
    my @amino=split//,$lines[1];
    for (my $i=0;$i<@amino;$i++){
        $pep{$lines[0]}[$i]=$amino[$i];
    }
}
close IN;
##################read CDS
my %codon;
open (IN, $fIn2) or die $!;
while (<IN>) {
	chomp;
    next if (/^$/);
	my @lines=split/\n/,$_,2;
    $lines[1]=~s/\n//g;
    my @base=split//,$lines[1];
    for (my $i=0;$i<@base;$i+=3){
        my $bs=$base[$i].$base[$i+1].$base[$i+2];
        my $c=$i/3;
        $codon{$lines[0]}[$c]=$bs;
    }
}
close IN;
mkdir $fOut if(! -d $fOut);
my $fname=basename($fIn2);
$fname=~s/check.fa/aln/;
open (OUT, ">$fOut/$fname") or die $!;
open (OUT2, ">$fOut/$fname.name.list") or die $!;
my $count=1;
open (IN, $fIn3) or die $!;######read pep alignment
while (<IN>) {
	chomp;
    next if (/^$/);
	my @lines=split/\n/,$_,2;
    $lines[1]=~s/\n//g;
    my @a=split//,$lines[1];
    my $index=0;
    print OUT ">seq$count\n";
    print OUT2 "seq$count\t$lines[0]\n";
    $count++;
    for (my $i=0;$i<@a;$i++){
        if($a[$i] eq "-"){
           print OUT "---"; 
        }
        else{
            print OUT "$codon{$lines[0]}[$index]";
            $index++;
        }
    }
    print OUT "\n";
}
close IN;
close OUT;
close OUT2;
#######################################################################################

