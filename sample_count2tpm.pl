#!usr/bin/perl -w
use strict;
my $count_file=shift;
my $sample_count_txt=shift;

open(IN,$sample_count_txt)||die"cannot open:$!";
my %sample_total_count;
while(<IN>)
{
	chomp;    my @a=split(/\t/,$_);    $sample_total_count{$a[0]}=$a[1];
}
close IN;

open(IN,$count_file)||die"cannot open:$!";
my %sample;
while(<IN>)
{
	chomp;    my @a=split(/\t/,$_);
	if($.==1)
	{
        my $count_name;
        my $tpm_name;
		foreach my $out(4 .. $#a)
		{
			$sample{$out}=$a[$out];
            $count_name .= "\t$a[$out]\_count";
			$tpm_name .= "\t$a[$out]\_TPM";
		}
		print join("\t",@a[0..3])."$count_name$tpm_name\n";
	}
	else{
        print "$_";
		foreach my $out(4 .. $#a)
		{
			if($a[$out]>0)
			{
				my $tmp=sprintf("%.4f",$a[$out]/$sample_total_count{$sample{$out}}*1000000);
				print "\t$tmp";
			}
			else
			{
				print "\t0.01";
			}
		}
		print "\n";
	}
}
close IN;
