#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: get_sample_count_for_exp.pl
#
#        USAGE: ./get_sample_count_for_exp.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: smhong (), smhong@genedenovo.com
# ORGANIZATION: R&D
#      VERSION: 1.0
#      CREATED: 06/09/2020 03:18:34 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use FindBin qw/$Bin $Script/;
#use File::Basename qw/basename dirname/;
#use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 



    e.g: perl $0 <> <>

USAGE

die "$usage" if(@ARGV<1); 

my $annot_stat = shift;
my $count_list = shift;

my %count_list = map{$_=>1}split/,/,$count_list;

open my $in_fh,"$annot_stat" or die "$!:$annot_stat\n";
my @cols;
while(<$in_fh>){
    chomp;
    my @aa = split/\t/,$_;
    # sample    total   rRNA    scRNA   snRNA   snoRNA  tRNA    exon_sense  exon_antisense  intron_sense    intron_antisenserepeat  exist_mirna exist_mirna_edit    known_mirna novel_mirna transcriptome   unann
    # all   47503187    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    0(0.00%)    32887124(69.23%)    2454651(5.17%)  1422344(2.99%)  0(0.00%)0(0.00%)    0(0.00%)
    if($. == 1){
        for my $i(0..$#aa){
            $aa[$i] =~ s/_abundance//;
            push @cols,$i if(exists $count_list{$aa[$i]});
        }
        next;
    }
    print join("\t",$aa[0],Sum(@aa[@cols]))."\n";
}
close $in_fh;

sub Sum{
    my @values = @_;

    my $sum = 0;
    for my $v(@values){
        $v =~ s/\s*\(.*//;
        $sum += $v;
    }

    return $sum;
}

