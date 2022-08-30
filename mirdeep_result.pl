#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: mirdeep_result.pl
#
#        USAGE: ./mirdeep_result.pl  
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
#      CREATED: 07/13/2020 03:55:22 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use FindBin qw/$Bin $Script/;
#use File::Basename qw/basename dirname/;
#use Getopt::Long;
#use lib $Bin;
use Data::Dumper;
use utf8;
use Carp qw(carp confess croak);

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 



    e.g: perl $0 <output.mrd> <result_csv> <outdir>

USAGE

die "$usage" if(@ARGV<1); 

my $output_mrd = shift;
my $result_csv = shift;
my $pdf_dir = shift;
my $outdir = shift;
my $run = shift;

$run //= 1;

my $name = "novel_mirna";
my %idmap = MrdResult($output_mrd,"$outdir/$name.list","$outdir/$name.tag_annot.txt","$outdir/$name.aln.txt");
#MirnaFa($mirna_fa,"$outdir/$name.fa",\%idmap);
TableResult($result_csv,"$outdir/$name.table.xls","$outdir/$name.fa",\%idmap);
Structure($pdf_dir,"$outdir/structure",\%idmap);

sub MrdResult{
    my $file = shift;
    my $outfile = shift;
    my $outfile2 = shift; ## $outdir/novel_mirna.tag_annot.txt
    my $out_aln_txt = shift;

    open my $in_fh,"$file" or die "$!:$file\n";
    open my $out_fh,">$outfile" or die "$!:$outfile\n";
    open my $out_aln_fh,">$out_aln_txt" or die "$!:$out_aln_txt\n";
    my $id;
    my ($mature_start,$mature_end,$seq);
    my ($total,$tag_total);
    my ($left_cut,$right_cut);
    my @tags;
    my $arm;
    my $mirna_id_tmp = "m0000";
    my $mirna_id;
    my %tag_info;
    my %idmap;
    my $score = 0;
    $/ = ">";
    <$in_fh>;
    while(<$in_fh>){
        chomp;
#        print "$_\n";
# exp              ffffffffffffffffffffMMMMMMMMMMMMMMMMMMllllllllSSSSSSSSSSSSSSSSSSSSSSffffffffffffffffffffffffffffff
# pri_seq          guggcaggaacauagaauggaggggaaucugacugucuaaugacugacaaugauugguuccaaccuucaugaaggugguguguggcaagagcacagug
# pri_struct       ((.(((.(..(((...(((((((((((((....((((........))))......))))))..)))))))....))).).))).))....(((((...
# ttt_00592819_x5  ...................gaggggaaucugacugucu............................................................
        my %info = GetInfo($_);
        my $score = $info{'score total'};
        next if($score < 0);
        if($info{exp} =~ /(f*)(M+)(l*)(S+)(f*)/){
            $mature_start = length($1)+1;
            $mature_end = $mature_start + length($2)-1;
            ($left_cut,$right_cut) = (2,5);
            $arm = "5p";
#            print join("\t",$mature_start,$mature_end,$1,$2,$3,$4,$5)."\n";
        }elsif($info{exp} =~ /(f*)(S+)(l*)(M+)(f*)/){
            $mature_start = length($1)+length($2)+length($3) + 1;
            $mature_end = $mature_start + length($4)-1;
            ($left_cut,$right_cut) = (2,5);
            $arm = "3p";
#            print join("\t",$mature_start,$mature_end,$1,$2,$3,$4,$5)."\n";
        }
        $mirna_id_tmp++;
        my $id = $info{id};
        my $mirna_id = "novel-$mirna_id_tmp-$arm";
        $idmap{$id} = $mirna_id;
        $idmap{structure}{$id} = $info{pri_struct};
        my $total = $info{'mature read count'};
        my $seq = FormatSeq(substr($info{pri_seq},$mature_start-1,$mature_end-$mature_start+1));
        my $tag_total = 0;
        my @tags;
        for my $tag(@{$info{tag}}){
            $tag =~ /^ttt_(\d+)_x(\d+)\s+(\.*)([augc]+)/;
            my $tag_start = length($3) + 1;
            my $tag_end = $tag_start + length($4) - 1;
            my $tag_id = "t$1";
            my $tag_seq = $4;
#            print join("\t",$tag_start,$tag_end,$1,$2,$3,$4,$mirna_id)."\n";
            if($tag_start >= $mature_start - $left_cut && $tag_end <= $mature_end + $right_cut){
                $tag_total += $2;
                push @tags,$tag_id;
            }else{
#                print STDERR join("\t",$id,$mature_start,$mature_end,$tag_start,$tag_end)."\n";
            }
            $tag_info{$tag_id}{num} = $2;
            $tag_info{$tag_id}{len} = $tag_end - $tag_start + 1;
            $tag_info{$tag_id}{seq} = uc($tag_seq);
            $tag_info{$tag_id}{seq} =~ s/U/T/g;
            push @{$tag_info{$tag_id}{mirnas}},$mirna_id;
        }
#        die "$id\t$total\t$tag_total" if($total != $tag_total);
        print $out_fh join("\t","$mirna_id",$seq,join(",",sort @tags))."\n";
        print $out_aln_fh ">$mirna_id\n$info{line}\n";
    }
    $/ = "\n";
    close $in_fh;
    close $out_fh;
    close $out_aln_fh;
#    print STDERR "$id\t$mirna_id\n";

    open my $out2_fh,">$outfile2" or die "$!:$outfile2\n";
    for my $tag_id(sort keys %tag_info){
        my $mirnas = join(",",sort @{$tag_info{$tag_id}{mirnas}});
        print $out2_fh join("\t",$tag_id,(map{$tag_info{$tag_id}{$_}}qw/len num seq/),"novel_mirna",$mirnas)."\n";
    }
    close $out2_fh;

    return %idmap;
}

sub GetInfo{
    my $line = shift;

    my @infos = split/\n/,$line;
    my %hash;
    $hash{id} = shift @infos;
    for my $info(@infos){
        next if($info eq '');
        next if($info =~ /^\s*$/);
#        next if($info =~ /^obs/);
        if($info =~ /^ttt_\d+/){
            push @{$hash{tag}},$info;
            $info =~ s/ttt_/t/;
            $info =~ s/\s/    /;
        }else{
            my ($key,$value);
            if($info =~ /\t\s/){
                ($key,$value) = $info =~ /(.*?)\t\s+(\S+)/;
            }else{
                ($key,$value) = split/\s+/,$info;
            }
            die "$info" if(!defined $key); 
            $hash{$key} = $value;
        }
        $hash{line} .= "$info\n";
    }
#    print Dumper(%hash);die;

    return %hash;
}

sub MirnaFa{
    my $infa = shift;
    my $out_fa = shift;
    my $idmap_ref = shift;

    open my $in_fh,"$infa" or die "$!:$infa\n";
    open my $out_fh,">$out_fa" or die "$!:$out_fa\n";
    while(<$in_fh>){
        chomp;
        if(/>(\S+)/){
            print $out_fh ">$idmap_ref->{$1}\n";
        }else{
            print $out_fh "$_\n";
        }
    }
    close $in_fh;
    close $out_fh;

    return;
}

sub TableResult{
    my $csv = shift;
    my $outfile = shift;
    my $outfa = shift;
    my $idmap_ref = shift;

    open my $in_fh,"$csv" or die "$!:$csv\n";
    open my $out_fh,">$outfile" or die "$!:$outfile\n";
#    open my $out_fa_fh,">$outfa" or die "$!:$outfa\n";
    my $no_need = 1;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
# provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive       rfam alert      total read count        mature read count       loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed        UCSC browser    NCBI blastn     consensus mature sequence       consensus star sequence consensus precursor sequence    precursor coordinate
# 7_2130  145.4           -       283     277     0       6       yes     -       -       -       -       uggccuaaccgaacccuguuc   aacaggggcggggcgacg      aacaggggcggggcgacgucaccuuaccuuuucugguggccuaaccgaacccuguuc       7:28378500..28378557:-
        if(/provisional/){
            $no_need = 0;
            print $out_fh join("\t",qw/mature_id hairpin_id genomic_id hairpin_start hairpin_end hairpin_strand hairpin_length(nt) miRDeep2_score hairpin_GC mature_seq mature_length(nt) hairpin_seq hairpin_struct/)."\n";
            next;
        }
        next if($no_need == 1);
        die "No mirna_id info of [$aa[0]]" if(!exists $idmap_ref->{$aa[0]});
        my $mirna_id = $idmap_ref->{$aa[0]};
        my ($hairpin_id,$arm) = $mirna_id =~ /(.*)-([35]p)/;
        my @pos = $_ =~ /(\S+):(\d+)\.\.(\d+):(\S+)/; ## chr start end strand
        my $len = $pos[2] - $pos[1] + 1;
        my $mature_seq = FormatSeq($aa[13]);
        my $mirna_len = length($mature_seq);
        my $hairpin_seq = FormatSeq($aa[15]);
        my $hairpin_GC = CalGC($hairpin_seq);
        print $out_fh join("\t",$mirna_id,$hairpin_id,@pos,$len,$aa[1],$hairpin_GC,$mature_seq,$mirna_len,$hairpin_seq,$idmap_ref->{structure}{$aa[0]})."\n";
#        print $out_fa_fh ">$mirna_id\n$mature_seq\n";
    }
    close $in_fh;
    close $out_fh;
#    close $out_fa_fh;

    return;
}

sub Structure{
    my $pdf_dir = shift;
    my $outdir = shift;
    my $idmap_ref = shift;

    mkdir $outdir if(!-e $outdir);
    my @pdfs = glob("$pdf_dir/*pdf");

    for my $pdf(@pdfs){
#        print "$pdf\n";
        my $name = $pdf;
        $name =~ s/.*\///;
        $name =~ s/.pdf$//;
        next if(!exists $idmap_ref->{$name});
        Excu("ln -sf $pdf $outdir/$idmap_ref->{$name}.pdf");
    }

    return;
}

sub CalGC{
    my $seq = shift;
    
    my $gc_num = $seq =~ tr/GC/GC/;

    my $len = length($seq);

    my $gc = sprintf("%.2f",$gc_num/$len*100);

    return $gc;
}

sub FormatSeq{
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/U/T/;

    return $seq;
}

sub FilterLine{
    my $line = shift;

    my @aa = split/\s+/,$line;
    my %cc = map{$_=>1} qw/score total loop star pri_struct obs/;

    my $out = exists $cc{$aa[0]} ? 1 : 0;

    return $out;
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    if($run == 1){
        my $ret = system($cmd);
        Error("Error $ret : [$cmd]") if($ret != 0);
    }
    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}

