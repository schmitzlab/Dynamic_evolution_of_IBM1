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
my ($fIn,$fOut,$fIn2);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"i2:s"=>\$fIn2,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $fIn2);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2019/10/15
Description:	this program is used to check if amino acid are match with the reading frame  of cds...,filter the bad seq and give the final statistics
Usage:
  Options:
  -i  <file>  input file,pep,selected,fa,forced  
  -i2 <file>  input file,cds,fa,forced  
  -o  <file>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
##################read pep
$/="\>";
open (IN, $fIn) or die $!;
my %pep;
while (<IN>) {
	chomp;
    next if (/^$/);
	my @lines=split/\n/,$_,2;
    my @info=split/\-/,$lines[0];
    $info[0]=~s/.p$//;
    $lines[0]=join("-",@info);
    $lines[1]=~s/\n//g;
    $lines[1]=~s/\*$//;
    $lines[1]=uc $lines[1];
    $pep{$lines[0]}=$lines[1];
}
close IN;
##################read CDS
open (IN, $fIn2) or die $!;
my %cds;
while (<IN>) {
	chomp;
    next if (/^$/);
	my @lines=split/\n/,$_,2;
    $lines[1]=~s/\n//g;
    $lines[1]=uc $lines[1];
    if($pep{$lines[0]}){
        $cds{$lines[0]}=$lines[1];
    }
}
close IN;
########################
 my(%gcode) = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

####check if the position between amic acid and base are right between codon and peptide
my %selectaa;### the record the pep sequence that pass the check
my %selectcc;### the record the condon sequence that pass the check
my %s1;### record sequence without CDS seq which cause match error between pep and codon seq filter
###s2;### record pep sequence using X to mark stop codon  cause match error between pep and codon seq, address
my %s3;### record sequence whose pep and codon are not match with no other reason, only a few species, guess from some error, just filter
###s4 if it has N in condon sequence,some sequence using X to represent unknow amino, filter those sites when do the screen.
my %frtout;###record the number of sequence filtered out from this script
###########
foreach my $name (keys %pep){
    my @amino=split//,$pep{$name};
    my @base=split//,$cds{$name};
    my $count=0;
    my $aaseq;###record the amino acid that pass the check
    my $codonseq;###record the codon that pass the check
    #4 species has weird stuations, just exclude them. record this situation, do the filter
    for (my $i=0;$i<@amino;$i++){
        my $c=$i*3;###index for condon
        my $as;my $bs;
        if(!$base[$c+1]||!$base[$c+2]){
            $s1{$name}++;
            $frtout{$name}++;
            last;
        }
        else {
            $bs=$base[$c].$base[$c+1].$base[$c+2];
            $as=$amino[$i];
            if($gcode{$bs}){
                if($as eq $gcode{$bs}){
                    $aaseq.=$as;
                    $codonseq.=$bs;
                    $count++;
                }
                elsif($i==$#amino&&$as eq "X"){
                    $count++;
                    
                }
                else{
                    $s3{$name}++;
                    $frtout{$name}++;
                    print "$name\n";
                    print $aaseq,"\n";
                    print $codonseq,"\n";
                    print $as,"\n";
                    print $bs,"\n";
                    last;
                }
            }
        }       
    }
    $selectaa{$name}=$aaseq;
    $selectcc{$name}=$codonseq;
}

mkdir $fOut if(! -d $fOut);
my $fname=basename($fIn);
$fname=~s/\.pep\.fa//;
open (OUT, ">$fOut/$fname.pep.check.fa") or die $!;
open (OUT2, ">$fOut/$fname.codon.check.fa") or die $!;
open (OUT3, ">$fOut/$fname.filterout.stat") or die $!;
foreach my $id (keys %selectaa){
    if(!$frtout{$id}){
        print OUT ">$id\n";
        print OUT "$selectaa{$id}\n";
        print OUT2 ">$id\n";
        print OUT2 "$selectcc{$id}\n";
    }
}
close OUT;
close OUT2;

my $out=keys %frtout;
my $all=keys %selectaa;
my $keep=$all-$out;
print OUT3 "All:$all\n";
print OUT3 "Keep:$keep\n";
print OUT3 "OUT:$out\n";
print OUT3 "###################without CDS seq\n";
foreach my $id (keys %s1){
    print OUT3 "$id\n";
}
print OUT3 "###################not match\n";
foreach my $id (keys %s3){
    print OUT3 "$id\n";
}
close OUT3;
#######################################################################################

