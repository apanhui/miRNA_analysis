#!usr/bin/perl -w
use strict;
my $txt_list=shift;
my $sample_list=shift;
my $out_txt = shift;
my $out_fa = shift;
my $low_cutoff = shift;
my $type = shift;

$low_cutoff //= 0;
my $len_cut = 18;
my @files=split(/,/,$txt_list);
my @sample=split(/,/,$sample_list);

if(-e $sample_list){
    open ID,$sample_list or die $!;
    @sample = <ID>;
    chomp(@sample);
    close ID;
    @files = map{"$txt_list/$_/clean.txt"}@sample;
}
#$low_cutoff = $low_cutoff eq 'yes' ? scalar(@sample)/2 : 1;
if($low_cutoff eq 'yes'){
    $low_cutoff = $type eq 'plant' ? 1 : scalar(@sample)/2;
}
print "low_cutoff $low_cutoff\n";

my %hash;    my %total;
foreach my $out( 0 .. $#files) {
	open(IN,$files[$out])||die"cannot open:$!";
	while(<IN>) {
		chomp;    my @a=split(/\s+/,$_);
		$hash{$a[2]}{$sample[$out]}+=$a[1];
		$total{$a[2]}+=$a[1];
	}
	close IN;
}
open OUT,">$out_txt" or die $!;
open FA,">$out_fa" or die $!;
open FA2,">$out_fa.filter" or die $!;
print OUT "id\tlength\tseq\ttotal_count\t",((join"\t",@sample)),"\n";
#open LOW,">$out_txt.low_cutoff" or die $!;
#print LOW "id\tlength\tseq\ttotal_count\t",((join"\t",@sample)),"\n";
my $id="t00000001";
foreach my $out(sort {$total{$b} <=> $total{$a}} sort keys %total) {
	my $len=length $out;
	my $out_line =  "$id\t$len\t$out\t$total{$out}";
    foreach my $in (@sample) {
        if(exists $hash{$out}{$in}) {
            $out_line .= "\t$hash{$out}{$in}";
        } else {
            $out_line .= "\t0";
        }
    }

#    print OUT "$out_line\n";
    if($len >= $len_cut){
        if($total{$out} > $low_cutoff){
            print OUT "$out_line\n";
            print FA ">$id\t$total{$out}\n$out\n";
        }else{
            print OUT "$out_line\n";
        }
    }else{
        print FA2 ">$id\t$total{$out}\n$out\n";
#        print LOW "$out_line\n";
    }
    $id++;
}
close IN;
close OUT;
close FA;
close FA2;
#close LOW;
