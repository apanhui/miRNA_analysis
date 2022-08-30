#!usr/bin/perl -w
use strict;
my $file=shift;
open(IN,$file)||die"cannot open:$!";
my %hash;    my %fa;    my %count;   my %tag;
while(<IN>)
{
	chomp;    my @a=split(/\t/,$_);
#	if($a[2]=~/^\w\w\w-(MIR\d+)/)
#	{
#		$family=$1;
#	}
#	elsif($a[2]=~/^\w\w\w-(\w\w\w-\d+)/)
#	{
#		$family=$1;
#	}
#	elsif($a[2]=~/^(\w+)-(.*)/)
#	{
#		$family=$2;
#	}
#	push @{$tag{$a[0]}},$family;
    # t00000001 23  14670439    TAGCTTATCAGACTGATGTTGAC exist_mirna mmu-mir-21a mmu-miR-21a-5p  TAGCTTATCAGACTGATGTTGA  18
    my @ids = split/,/,$a[-3];
    my @seqs = split/,/,$a[-2];
    for my $i(0..$#ids){
        push @{$hash{$ids[$i]}},$a[0];
        $fa{$ids[$i]} = $seqs[$i];
    }
	$count{$a[0]}=$a[2];
}
close IN;

foreach my $out(sort keys %hash)
{
	print "$out\t";
	my @uniq_arr=keys %{ {map {$_ => 1} @{$hash{$out}}} };
	my @arr=sort {$count{$b}<=>$count{$a}} @uniq_arr ;
	print "$fa{$out}\t";
	foreach my $in(@arr)
	{
		print "$in,";
	}
	print "\n";
}

exit(0);
open(OUT,">miRNA-variant.out")||die"cannot open:$!";
foreach my $out(sort keys %hash)
{
	my @id=sort keys %{{map {$_=>1} @{$hash{$out}}}};
	my $total;
	foreach my $in(@id)
	{
		$total+=$count{$in}
	}
	if(@id>1)
	{
		print OUT "$out\t$total\t$fa{$id[0]}\t$id[0]\t$count{$id[0]}\t$fa{$id[1]}\t$id[1]\t$count{$id[1]}\n";
	}
	else
	{
		print OUT "$out\t$total\t$fa{$id[0]}\t$id[0]\t$count{$id[0]}\t-\t-\t-\n";
	}
}
close OUT;

open(OUT,">miRNA-variant.dat")||die"cannot open:$!";
my $total_count=0;
foreach my $out(sort {$count{$b}<=>$count{$a}} keys %tag)
{
	my @uniq=sort keys %{{map {$_=>1} @{$tag{$out}}}};
	print OUT "$out\t$count{$out}\t$fa{$out}\t",(join",",@uniq),"\n";
	$total_count+=$count{$out}
}
close OUT;

my $uniq_count=scalar keys %fa;
open(OUT,">miRNA-variant.stat")||die"cannot open:$!";
print OUT "mismatch: 2\n";
print OUT "number of uniq sequences: $uniq_count\n";
print OUT "number of total sequences: $total_count\n";
close OUT;
