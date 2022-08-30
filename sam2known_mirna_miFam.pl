#!usr/bin/perl -w
use strict;
use warnings;

die "\n    Usage: perl $0 <mirbase_file> <sam result>   >result\n\n" if(@ARGV < 2);

my $mirbase_file=shift;
#my $exist_file = shift;
my $file=shift;
my $out_mirna_list = shift;
#my $out_mirna_fa = shift;

#my %exist_info = ReadExist($exist_file);

open(IN,$mirbase_file)||die"cannot open:$!";
my %hash1;   my %hairpin_seq;
while(<IN>) {
	chomp;   my @a=split(/\s+/,$_);
	my $mature_start=index($a[4],$a[3])+1;
	my $mature_end=$mature_start+length($a[3])-1;
	if($mature_start>=0) {
		$hash1{$a[1]}{$a[0]}="$mature_start\t$mature_end";
	}
	$hairpin_seq{$a[1]}=$a[4];
}
close IN;

open(IN,$file)||die"cannot open:$!";
my %fa;    my %count;   my %tag1;   my %tag2;
while(<IN>) {
	chomp;  
    # t00000012 21  413993  TCCCAAATGTAGACAAAGCAT   known_mirna ath-MIR158a,aly-MIR158a ath-miR158a-3p,aly-miR158a-3p   TCCCAAATGTAGACAAAGCA,TCCCAAATGTAGACAAAGCA   63,61
    my @a=split(/\t/,$_);    
    my @mature = split/,/,$a[-4];
    my @mirnas = split/,/,$a[-3];
#    next if(Exist(\%exist_info,@mirnas));
#    print "$_\n";
    for my $i(0..$#mature){
        my $id = $mature[$i];
        my $family;
    	if($id=~/^\w+-(MIR\d+)/) {
    		$family=$1;
    		$family =~ s/MIR/miR/;  ##  plant name rule  --hsm
    	} elsif($id=~/^\w+-(mir-\d+)/) {
    		$family=$1;
    		$family =~ s/mir/miR/; ## animal name rule --hsm
        }elsif($id =~ /^\w+-(mir-\w.*)/){ # fungi hvt-mir-H14
            $family = $1;
            $family =~ s/mir/miR/; ## fungi name rule --hsm
        }elsif($id =~ /^\w+-(miR-\w.*)/){ # fungi kshv-miR-K12-2
            $family = $1;         
#            print "$id\t$family\n";
    	} elsif($id=~/^\w+-(bantam)/) {
    		$family=$1;
    	} elsif($id=~/^\S*(mir-iab-\d)/) {
    		$family=$1;
    	} elsif($id=~/^\w+-(\w\w\w-\d+)/) {
    		$family=$1;
    	} else {
    		print STDERR "error: sam2known_mirna.pl    $id not match\n";
    		exit(0);
    	}
        my $mirna = $mirnas[$i];
    	if($mirna=~/-3p/) {
    		push @{$tag1{$a[0]}{$family}},"3p";
    	} elsif($mirna=~/-5p/) {
    		push @{$tag1{$a[0]}{$family}},"5p";
    	} else {
    		push @{$tag1{$a[0]}{$family}},"np";
    	}
    }
    $fa{$a[0]}=$a[3];
    $count{$a[0]}=$a[2];
}
close IN;

my %hash;
foreach my $out(sort keys %tag1) {
	foreach my $in(sort keys %{$tag1{$out}}) {
		my @arms1=sort keys %{ {map {$_ => 1} @{$tag1{$out}{$in}}} };
		if(@arms1 ==1 && $arms1[0] eq "np") {
			push @{$hash{$in}{"np"}},$out;
		} elsif(@arms1 == 1 && $arms1[0] ne "np") {
			push @{$hash{$in}{$arms1[0]}},$out;
		} elsif(@arms1 == 2 && $arms1[1] eq "np") {
			push @{$hash{$in}{$arms1[0]}},$out;
		} elsif(@arms1 == 2 && $arms1[1] ne "np") {
			push @{$hash{$in}{"np"}},$out;
		} elsif(@arms1 ==3) {
			push @{$hash{$in}{"np"}},$out;
		}
	}
}

my %id=("5p"=>"x","3p"=>"y","np"=>"z");

open OUT,">$out_mirna_list" or die $!;
#open FA,">$out_mirna_fa" or die $!;
foreach my $out(sort keys %hash) {
	foreach my $in(sort keys %{$hash{$out}}) {
		print OUT "$out-$id{$in}\t";
		my @uniq_arr=keys %{ {map {$_ => 1} @{$hash{$out}{$in}}} };
		my @arr=sort @uniq_arr;
		print OUT "$fa{$arr[0]}\t";
#        print FA ">$out-$id{$in}\n$fa{$arr[0]}\n";
		foreach my $in(sort @arr) {
			print OUT "$in,";
		}
		print OUT "\t$out\n";
	}
}
close OUT;
#close FA;

sub ReadExist{
    my $exist_file = shift;

    my %exist;
    open EX,$exist_file or die $!;
    while(<EX>){
        chomp;
        my @aa = split/\t/,$_;
        $exist{$aa[1]} = 1;
    }
    close EX;

    return %exist;
}

sub Exist{
    my $info_ref = shift;
    my @ids = @_;

    my $found = 0;
    for my $id(@ids){
        next if(!exists $info_ref->{$id});
        $found = 1;
    }
    print "$found @ids\n";

    return $found;
}
