#!usr/bin/perl -w
use strict;
my $clean_txt=shift;
my $all_mirna_file=shift;
my $SamplesOrder=shift;
my @SamplesOrderLst=split(/,/,$SamplesOrder);
my $out_txt = shift;
my $out_fa = shift;

if(-e $SamplesOrder){
    @SamplesOrderLst = ();
    open my $in_fh,"$SamplesOrder" or die "$!:$SamplesOrder\n";
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        push @SamplesOrderLst,$aa[0];
    }
    close $in_fh;
}

my %e;
open ALL,$all_mirna_file or die $!;
while(<ALL>){
    chomp;
    my @aa = split;
    my @bb = split/,/,$aa[2];
    map{$e{$_} = 1} @bb;
}
close ALL;


open(IN,$clean_txt)||die "cannot open:$!";
my %sample;    my %tag;
while(<IN>)
{
	chomp;   my @a=split(/\s+/,$_);
	if($.==1)
	{
		foreach my $out(4 .. $#a)
		{
			$sample{$out}=$a[$out];
		}
	}
	else
	{
        next if(!exists $e{$a[0]});
		foreach my $out(4 .. $#a)
		{
			$tag{$a[0]}{$sample{$out}}=$a[$out];
		}
	}
}
close IN;

open(IN,$all_mirna_file)||die"cannot oepn:$!";
my %count;    my %total;     my %annot;
my %seq;
while(<IN>)
{
	chomp;    my @a=split(/\s+/,$_);
	my @b=split(/,/,$a[2]);
	foreach my $out(@b)
	{
		foreach my $in(keys %{$tag{$out}})
		{
			$count{$a[0]}{$in}+=$tag{$out}{$in};
			$total{$a[0]}+=$tag{$out}{$in}
		}
	}
	my $len=length $a[1];
	$annot{$a[0]}="$len\t$a[1]";
    $seq{$a[0]} = $a[1];
}
close IN;

open my $out_fh,">$out_txt" or die "$!:$out_txt\n";
open my $out_fa_fh,">$out_fa" or die "$!:$out_fa\n";

print $out_fh "id\tlength\tseq\ttotal\t",(join"\t",@SamplesOrderLst),"\n";
foreach my $out(sort {$total{$b}<=>$total{$a} or $a cmp $b} keys %count)
{
	print $out_fh "$out\t$annot{$out}\t$total{$out}";
    print $out_fa_fh ">$out   $total{$out}\n$seq{$out}\n";
	foreach my $in(@SamplesOrderLst){
        $count{$out}{$in} //= 0;
		print $out_fh "\t$count{$out}{$in}";
	}
	print $out_fh "\n";
}
close $out_fh;
close $out_fa_fh;
