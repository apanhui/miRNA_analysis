#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: remove_clean_fa_some_tag.pl
#
#        USAGE: ./remove_clean_fa_some_tag.pl  
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
#      CREATED: 06/06/2020 01:37:48 PM
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

my $clean_fa = shift;
my $pipe_dir = shift;
my $type = shift; # known novel 
my $filter_num = shift; ## < sample_num

my %dir;
open my $in_fh,"$pipe_dir" or die "$!:$pipe_dir\n";
while(<$in_fh>){
    chomp;
    my @aa = split/\t/,$_;
    $dir{$aa[0]} = $aa[1];
}
close $in_fh;

my @keys = qw/result_rfam_result result_ncgb_result result_exist_mirna_result result_exist_mirna_edit_result/; ## known_mirna 
if($type eq 'novel_mirna'){
    push @keys,qw/result_exon_sense_result result_repeat_result result_known_mirna_result/;
}

my %remove_list;
for my $key(@keys){

    if(!exists $dir{$key}){
        print STDERR "next $key\n";
        next;
    }else{
        print STDERR "Read $key [$dir{$key}]\n";
    }
    open my $in_fh,"$dir{$key}" or die "$!:$dir{$key}\n";
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        $remove_list{$aa[0]} = 1;
    }
    close $in_fh;
}

open FA,$clean_fa or die $!;
my $flag;
while(<FA>){
    chomp;
    if(/>(\S+)\s+(\d+)/){
        $flag = exists $remove_list{$1} ? 0 : 1;
        $flag = 0 if($flag == 1 && $2 < $filter_num);
    }
    next if($flag == 0);
    print "$_\n";
}
close FA;
