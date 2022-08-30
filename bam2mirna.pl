#!usr/bin/perl -w
use strict;

my $mirbase_file=shift;
my $sam_file=shift;
my $clean_file=shift;
my $title=shift;
my $output = shift;

open(IN,$mirbase_file)||die"cannot open:$!";
my %hash;   my %mature_seq;
while(<IN>) {
	chomp;   my @a=split(/\s+/,$_);
	my $mature_start=index($a[4],$a[3])+1;
	my $mature_end=$mature_start+length($a[3])-1;
	if($mature_start>=0) {
		$hash{$a[1]}{$a[0]}="$mature_start\t$mature_end";
	}
	$mature_seq{$a[0]}=$a[3];
}
close IN;

open(IN,$clean_file)||die"cannot open:$!:$clean_file";
my %count;   my %fa;    my $key;
while(<IN>) {
	chomp;
	if(/>/) {
		s/>//;    my @a=split(/\s+/,$_);   $key=$a[0];
		$count{$a[0]}=$a[1];
	} else {
		$fa{$key}.=$_;
	}
}

my $range=2;
my %annot_tag;
my %annot_tag_mature;
open(IN,$sam_file)||die"cannot open:$!:$sam_file";
while(<IN>) {
	chomp;
	next if(/^@/);
	my @a=split(/\t/,$_);
	next if($a[2] eq "*");
	if(exists $hash{$a[2]}) {
		my $align_start=$a[3];    my $align_end=$align_start+length($a[9])-1;
		foreach my $out(keys %{$hash{$a[2]}}) {
			my @b=split(/\t/,$hash{$a[2]}{$out});
			if($align_start>=$b[0]-$range && $align_end<=$b[1]+$range) {
#				print "$_\t$mature_seq{$out}\t$out\t$count{$a[0]}\n";
			    push @{$annot_tag{$a[0]}},$a[2];
			    push @{$annot_tag_mature{id}{$a[0]}},$out;
			    push @{$annot_tag_mature{seq}{$a[0]}},$mature_seq{$out};
			    push @{$annot_tag_mature{pos}{$a[0]}},$a[3];
			}
#			push @{$annot_tag{$a[0]}},$a[2];
		}
	}
}
close IN;

open(OUT,">$output")||die"cannot open:$!";
foreach my $out(sort keys %annot_tag) {
#	my @uniq_arr=sort keys %{ {map {$_ => 1} @{$annot_tag{$out}}} };
    my @uniq_arr = @{$annot_tag{$out}};
    my $mature = join(",",@{$annot_tag_mature{id}{$out}});
    my $seq = join(",",@{$annot_tag_mature{seq}{$out}});
    my $pos = join(",",@{$annot_tag_mature{pos}{$out}});
	my $len=length $fa{$out};
	print OUT "$out\t$len\t$count{$out}\t$fa{$out}\t$title\t",(join",",@uniq_arr),"\t$mature\t$seq\t$pos\n";
}
close OUT;
