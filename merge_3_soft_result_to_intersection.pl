#!usr/bin/perl -w
use strict;

die "Usage:perl $0 <mireap> <miranda> <targetscan> <outfile>\n" unless @ARGV == 5;

my $mireap = shift;
my $miranda = shift;
my $targetscan = shift;
my $outfile = shift;
my $outdir = shift;

my %hash;

if($outfile =~ /.gz$/){
    open OUT,"|gzip >$outfile" or die $!;
}else{
    open OUT,">$outfile" or die $!;
}

## targetscan
if($targetscan ne 'none'){
    open my $out_fh,">$outdir/targetscan.xls" or die "$!:$outdir/targetscan.xls\n";
    open A, $targetscan or die "Can't open $targetscan!";
    <A>;
    while(<A>) {
    # a_Gene_ID       miRNA_family_ID species_ID      MSA_start       MSA_end UTR_start       UTR_end Group_num       Site_type       miRNA in this species   Group_type      Species_in_this_group   Species_in_this_group_with_this_site_type       ORF_overlap
    # ENSMUST00000000033_3UTR let-7-x 10090   2362    2368    2362    2368    2       7mer-1a x       7mer-1a 10090           0
        s/_3UTR//;
        my @a = split(/\t/, $_);
        my @b = split(/\s+/,$a[0]);
        next if($a[0] eq 'a_Gene_ID');
        $hash{"$a[1]\t$b[0]"}{'t'}++;
        print $out_fh "$a[1]\t$b[0]\n";
        if($miranda eq 'none' && $mireap eq 'none'){
            print OUT "$_";
        }
    }
    close A;
    close $out_fh;
}

## miranda
if($miranda ne 'none'){
    open my $out_fh,">$outdir/miranda.xls" or die "$!:$outdir/miranda.xls\n";
    open A, $miranda  or  die "Can't open $miranda!";
    $/ = "\n>";
    while(<A>) {
        s/>$//; s/^>//; s/_3UTR//;
        my @a = split(/\s+/, $_);
    # >miR-21-x       ENSMUST00000000058_3UTR 151.00  -17.43  2 21    1650 1680       26      61.54%  69.23%
    #    miRNA:    3' acagTTGTAGTCAG-------ACTATTCGAt 5'
    #                     ||:||||||:       | ||||||| 
    #    Target:   5' ataaAATATCAGTTTAATAAATTATAAGCTt 3'
        $hash{"$a[0]\t$a[1]"}{'m'}++;
        print $out_fh "$a[0]\t$a[1]\n";
        if($mireap eq 'none'){
            if($targetscan eq 'none'){
                print OUT ">$_";
            }else{
                print OUT ">$_" if(exists $hash{"$a[0]\t$a[1]"}{'t'});
                delete $hash{"$a[0]\t$a[1]"};
            }
        }
    }
    $/ = "\n";
    close A;
    close $out_fh;
}

if($mireap ne 'none'){
    open my $out_fh,">$outdir/RNAhybrid.xls" or die "$!:$outdir/RNAhybrid.xls\n";
    open A, $mireap or die "Can't open $mireap!";
    $/ = "\n>";
#print "miRNA_name\tmiRNA_len\ttarget_name\ttarget_3UTR_len\tmatch_3UTR_pos\tMFE\tp-value\tprediction_value\n";
    while(<A>) {
        s/>$//; s/^>//; s/_3UTR//;
        /^(\S+)\s+\S+\s+(\S+)/;
        print $out_fh "$1\t$2\n";
        if($miranda ne 'none' || $targetscan ne 'none'){
            next unless exists $hash{"$1\t$2"};
            if($miranda ne 'none' && $targetscan ne 'none'){
                print OUT ">$_" if exists $hash{"$1\t$2"}{'m'} && exists $hash{"$1\t$2"}{'t'};
            }else{
                print OUT ">$_" if exists $hash{"$1\t$2"}{'m'} || exists $hash{"$1\t$2"}{'t'};
            }
            delete $hash{"$1\t$2"};
        }else{ ## mireap eq none && targetscan eq 'none'
            print OUT ">$_";
        }
    }
    $/ = "\n";
    close A;
    close $out_fh;
}
close OUT;
