#!/usr/bin/perl
#
# Modified by miaoxin 16-05-17
#
# for small RNA solexa sequencing

## small RNA classes:
# miRNA
# miRNA.novel, novel miRNAs, need further validation
# rRNAetc, include rRNA tRNA snRNA snoRNA scRNA srpRNA ...
# repeat, match known repeat elements
# mRNA, maybe siRNA candidate or mRNA degradation species
# piRNA, piwi-interacting RNA
# intron ?
# unann, unannotated

use strict;
use warnings;
use Data::Dumper;
#use FindBin qw($Bin $Script);

die "usage: $0 <tag_txt> <config> <outfile>\n" if(@ARGV<3);
my $tag_file = shift;
my $pipe_info = shift;
my $outfile = shift;

my %pipe_info;
open my $in_fh,"$pipe_info" or die "$!:$pipe_info\n";
while(<$in_fh>){
    chomp;
    my @aa = split/\t/,$_;
    $pipe_info{$aa[0]} = $aa[1];
}
close $in_fh;
    
my @keys = qw/exist_mirna exist_mirna_edit ncgb    rfam    known_mirna novel_mirna repeat exon_sense genome/;
my @cols = qw/5           3                5       5       5           5           3      6    2/;
my @words= qw/exist_mirna exist_mirna_edit rRNAetc rRNAetc known_mirna novel_mirna repeat exon_sense genome_others/;

my %annot;
for my $i(0..$#keys){
    my $key = "result_$keys[$i]_result";
    if(exists $pipe_info{$key}){
        print "Read $key [$pipe_info{$key}] [$cols[$i]]\n";
    }else{
        print "next $key\n";
        next;
    }
#    next;
    my $file = $pipe_info{$key};
    open my $in_fh,"$file" or die "$!:$file\n";
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        next if(exists $annot{$aa[0]});
        my $col = $cols[$i] - 1;
        my $class;
        if($keys[$i] eq 'genome'){
            $aa[$col] = $words[$i];
            $class = $words[$i];
        }else{
            $class = Class($aa[$col],$words[$i]);
        }
        $annot{$aa[0]} = "$aa[$col]\t$class\t$words[$i]";
    }
    close $in_fh;
}

open TAG,$tag_file or die $!;
open OUT,">$outfile" or die "$!:$outfile";
print OUT "id\tannot\tclass\tsource\n";
my $id;
<TAG>;
while(<TAG>){
    chomp;
    my @aa = split/\t/,$_;
    $id = $aa[0];
    $annot{$id} //= "unann\tna\tna";
    print OUT "$id\t$annot{$id}\n";
}
close TAG;
close OUT;

sub Class{
    my $ann = shift;
    my $cls = shift;

    if ($ann=~/piRNA/i) { 
        $cls="piRNA"; 
    }elsif($ann=~/rRNA/){
        $cls="rRNA";
    }elsif($ann=~/tRNA/){
        $cls="tRNA";
    }elsif($ann=~/snRNA/){
        $cls="snRNA";
    }elsif($ann=~/snoRNA/){
        $cls="snoRNA";
    }elsif($ann=~/srpRNA/){
        $cls="srpRNA";
    }elsif($ann=~/scRNA/){
        $cls="scRNA";
    }

    return $cls;
}
