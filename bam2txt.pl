#!/usr/bin/perl
# for small RNA solexa sequencing

use strict;

die "usage: $0 <clean.fas> <*.soap>\n" unless @ARGV==2;

my $tag_file=shift;
my $soap_file=shift;

my %tagId_tag;
my %tagId_number;
&read_tag_file($tag_file,\%tagId_tag,\%tagId_number);

my $tot_n=0;
my $tot_c=0;
foreach(keys %tagId_number){
	$tot_n++;
	$tot_c+=$tagId_number{$_};
}

my %tagId_locs;
my $mat_n=0;
my $mat_c=0;
my $tag_name="NO_1";
open IN, $soap_file or die $!;
# t00000001 GTCAACATCAGTCTGATAAGCTA IIIIIIIIIIIIIIIIIIIIIII 1   a   23  -   11  86584119     23M    23
# bam 
# t00000041 16  8   123548954   255 19M *   0   0   AAAGCCTACAGCACCCGGT IIIIIIIIIIIIIIIIIII XA:i:0  MD:Z:19 NM:i:0
while (<IN>) {
	chomp;
	my @d=split;
    next if($d[1] == 4); ## unmap 
	if($tag_name ne $d[0]){
		$mat_n++;
		$mat_c+=$tagId_number{$d[0]};
	};
	$tag_name=$d[0];
    my $chr = $d[2];
    my $start = $d[3];
    my $seq = $tagId_tag{$d[0]};
    my $end = $start + length($seq) - 1;
    my $strand = $d[1] == 16 ? "-" : "+";

    # t00000048 11  74319892    74319911    +   ATACCGGGTGCTGTAGGCTT    181026
	print join ("\t",$d[0],$chr,$start,$end,$strand,$seq,$tagId_number{$d[0]}),"\n";
}
close IN;

print STDERR "\tunique\ttotal\n";
print STDERR "total\:\t$tot_n\t$tot_c\n";
print STDERR "match\:\t$mat_n\t$mat_c\n";

sub read_tag_file{
	my $infile=shift;
	my $tagId_tag=shift;
	my $tagId_number=shift;
	
	open IN, $infile || die $!;
	my $unique=0; my $total=0;
	while (<IN>) {
		chomp;
		if(/^>(\S+)\s+(\d+)/){
			my $id=$1;
			my $number=$2;
			my $tag = <IN>;
			chomp $tag;
			$tagId_tag->{$id}=$tag;
			$tagId_number->{$id}=$number;
		}
	}
	close IN;
}
