#!/usr/bin/perl

#Modified by Miaoxin 16-05-27

# for small RNA solexa sequencing
use strict;
use FindBin '$Bin';
use lib "$Bin/lib";
use Trim;

# options
use Getopt::Std;
use vars qw($opt_i $opt_d $opt_a $opt_b $opt_m $opt_p $opt_c $opt_e $opt_s $opt_r $opt_h $opt_z $opt_l);
getopts('i:d:a:b:m:p:c:e:s:r:z:l:h');
my $fq_file        = $opt_i;
my $work_dir       = $opt_d ? $opt_d : "./";
my $a5             = $opt_a ? $opt_a : "GTTCAGAGTTCTACAGTCCGACGATC";
my $a3             = $opt_b ? $opt_b : "CTGTAGGCACCATCA";
my $min_tag_length = defined $opt_m ? $opt_m : 18;
my $phread         = $opt_p ? $opt_p : 33;
my $low_cutoff     = defined $opt_c ? $opt_c : 1;
my $er             = $opt_e ? $opt_e : 0.01;
my $strict         = $opt_s ? $opt_s : 0;
my $rmbad          = $opt_r ? 1 : 0;
my $help           = $opt_h ? 1 : 0;
my $random_adaper  = defined $opt_z ? $opt_z : 'no';

my $seq_len ||= 10;

if ($help) { &usage(); exit; }

unless (-e $fq_file) { &usage(); exit; }

# work dir
$work_dir.="/" unless ($work_dir=~/\/$/);
mkdir $work_dir unless (-e $work_dir);

### remove low quality reads
my %seq_number;
open STDERR,">$opt_l" or die $! if(defined $opt_l);
&get_unique_seq($fq_file,\%seq_number,$phread,$strict,$work_dir);

### cal total reads number
my $tot_read_n=0;
foreach my $seq (keys %seq_number) {
	$tot_read_n += $seq_number{$seq};
}
die "[error] $tot_read_n = $tot_read_n\n" unless ($tot_read_n);
#my $count=sprintf "%.4f",$tot_read_n/$tot_read_n*100;
#print STDERR "high_quality\t$tot_read_n\t$count\%\n";

### trim a3
my %tag_number;
my ($na3_read_n,$nin_read_n) = trimA3(\%seq_number,$a3,\%tag_number,$work_dir,$random_adaper);
my $na3_read_p = sprintf "%.4f", $na3_read_n/$tot_read_n*100;
my $nin_read_p = sprintf "%.4f", $nin_read_n/$tot_read_n*100;
print STDERR "3'adapter_null\t$na3_read_n\t$na3_read_p%\n";
print STDERR "insert_null\t$nin_read_n\t$nin_read_p%\n";

### remove a5 contaminants
my $a5_read_n = rmA5(\%tag_number,$a5,$work_dir);
my $a5_read_p = 0;
$a5_read_p = sprintf "%.4f", $a5_read_n/$tot_read_n*100 if ($tot_read_n);
print STDERR "5'adapter_contaminants\t$a5_read_n\t$a5_read_p%\n";

### remove random adapter 
my %tag_number_new;
if($random_adaper eq 'yes'){
    for my $tag(keys %tag_number){
        my $tag_new = substr($tag,4,length($tag)-8);
#        print "$tag\n$tag_new\n";die;
        $tag_number_new{$tag_new} += $tag_number{$tag};
    }
    %tag_number = %tag_number_new;
}

### plot lengthVSnumber
my %length_number;
foreach my $tag (keys %tag_number) {
	my $length = length $tag;
	my $number = $tag_number{$tag};
	$length_number{$length} += $number;
}
&plot_lengthVSnumber(\%length_number,$work_dir);

### remove too small and polyA
my %clean_number;
my $small_n=0;  # small size
my $polyA_n=0;
my $clean_n=0;  # ideal size
my $small_file=$work_dir."small.txt";
open SIN, ">$small_file" || die $!;
foreach my $tag (sort {$tag_number{$b} <=> $tag_number{$a}} keys %tag_number){
	my $length = length $tag;
	my $number = $tag_number{$tag};
	if ($length < $min_tag_length){ # small insert
		print SIN "$length\t$number\t$tag\n";
		$small_n += $number;
	} else { # ideal insert
		if (isPolyA($tag)){ # polyA
			$polyA_n += $number;
		} else { # clean
			$clean_number{$tag}=$number;
			$clean_n += $number;
		}
	}
}
close SIN;


# write clean
my $clean_file = $work_dir."clean.txt";
&print_clean(\%clean_number,$clean_file);
my $small_p = sprintf "%.4f", $small_n/$tot_read_n*100 if ($tot_read_n);
my $polyA_p = sprintf "%.4f", $polyA_n/$tot_read_n*100 if ($tot_read_n);
my $clean_c = scalar keys %clean_number;
my $clean_p = sprintf "%.4f", $clean_n/$tot_read_n*100 if($tot_read_n);
print STDERR "smaller_than_$min_tag_length"."nt\t$small_n\t$small_p%\n" if($min_tag_length > 0);
print STDERR "polyA\t$polyA_n\t$polyA_p%\n";
#print STDERR "clean_reads_all\t$clean_n\t$clean_p%\n";

### remove bad
if ($rmbad) {
	print STDERR "remove bad tags from clean.txt:\n";
	my ($bad_c,$bad_n) = rmbad(\%clean_number,$er,$work_dir);
	$clean_c -= $bad_c;
	$clean_n -= $bad_n;
	# write clean nobad
	my $nobad_clean_file = $work_dir."clean_nobad.txt";
	&print_clean(\%clean_number,$nobad_clean_file);
	print STDERR "bad tags: $bad_n ($bad_c)\n";
	print STDERR "reads number[clean_nobad.txt]: $clean_n ($clean_c)\n\n";
}

### remove redundant tags
#print STDERR "remove redundant tags from clean_nobad.txt:\n";
#my ($redundant_c,$redundant_n)=rmredundant(\%clean_number,$work_dir);
#$clean_c-=$redundant_c;
#$clean_n-=$redundant_n;
## write clean nobad nr
#my $nr_clean_file=$work_dir."clean_nobad_nr.txt";
#&print_clean(\%clean_number,$nr_clean_file);
#print STDERR "redundant tags: $redundant_n ($redundant_c)\n";
#print STDERR "reads number[clean_nobad_nr.txt]: $clean_n ($clean_c)\n\n";

### filter low abundant
if ($low_cutoff) {
	my $low_c;
	my $low_n;
	my $low_file=$work_dir."clean_lteq"."$low_cutoff".".txt";
	open LT, ">$low_file" || die $!;
	foreach my $tag (sort {$clean_number{$a} <=> $clean_number{$b}} keys %clean_number) {
		if ($clean_number{$tag} <= $low_cutoff) {
			++$low_c;
			$low_n += $clean_number{$tag};
			my $len = length $tag;
			print LT join ("\t",$len,$clean_number{$tag},$tag),"\n";
			delete $clean_number{$tag};
		} else { last; }
	}
	close LT;
	$clean_c -= $low_c;
	$clean_n -= $low_n;
	# output clean
	my $clean_file_new = $work_dir."clean_nolow.txt";
	&print_clean(\%clean_number,$clean_file_new);
#	print STDERR "low abundant (<=$low_cutoff) tags $low_n ($low_c)\n";
	my $low_pp = sprintf "%.4f", $low_n/$tot_read_n*100 if($low_cutoff);
#	my $clean_pp = sprintf "%.4f", $clean_n/$tot_read_n*100 if($tot_read_n);
	print STDERR "low cutoff\t$low_n\t$low_pp%\n";
#	print STDERR "clean_tags\t$clean_n\t$clean_pp%\n";
	`mv -f $work_dir"clean_nolow.txt" $work_dir"clean.txt" `;
}
my $clean_pp = sprintf "%.4f", $clean_n/$tot_read_n*100 if($tot_read_n);
print STDERR "clean_tags\t$clean_n\t$clean_pp%\n";
#
sub get_unique_seq{
	my $file = shift;
	my $seq_number = shift;
	my $phread = shift;
	my $strict = shift;
	my $work_dir = shift;
	
	my $total_read_n = 0;
	my $lowq_read_n = 0;
	my $lowq_file = $work_dir."low_quality.txt";
	
	if ($file =~ /\.gz$/) {open (IN,"gzip -cd <$file|") || die "Can't open the gz file $file:$!";}
	else{open IN,'<',$file || die "Can't open $file:$!";}

	open LQ, ">$lowq_file" || die $!;
	while (<IN>) {
		my $seq = <IN>; chomp $seq; <IN>; my $qua = <IN>; chomp $qua;
		die "[error] fastq format infile is forced. Break at $.\n" if ($seq =~ /[^ATCGN]/);
		++$total_read_n;

		my $lowq = islowq($seq,$qua,$phread,$strict);
		if ($lowq) {
			++$lowq_read_n;
			print LQ "$seq\t$qua\n";
		} else {
			$seq_number->{$seq}++;
		}
	}
	close IN;
	close LQ;
	
	my $high_read_n = $total_read_n-$lowq_read_n;
	my $high_read_p = sprintf "%.4f", $high_read_n/$total_read_n*100;
	
	print STDERR "type\tcount\t%\n";
	print STDERR "clean_reads\t$total_read_n\t100%\n";
	print STDERR "high_quality\t$high_read_n\t$high_read_p\%\n";
}

close STDERR if(defined $opt_l);

# for filter low qulity reads
# 1. N number > 0
# 2. 1-10: [qv < 20] > 1, strict
# 3. 1-10: [qv < 10] > 0, strict
sub islowq {
	my($seq,$qua,$phread,$strict)=@_;
	my $lowq=0;

	# N number
	my $n = $seq =~ tr/N//;
	if ($n > 0) { $lowq = 1; return $lowq; }

	my $seq_now;
	if(length $seq >= $seq_len) { $seq_now = substr($seq,0,$seq_len); } else { $seq_now = $seq; }

	my $lt10_n = 0;
	my $lt20_n = 0;
	foreach my $i (0 .. (length($seq_now)-1)) {
		my $qv = ord(substr($qua,$i,1))-$phread;
		if ($qv < 10) { ++$lt10_n; }
		if ($qv < 20) { ++$lt20_n; }
	}

	if ($strict) {
		if ($lt10_n > 0) { $lowq = 1; return $lowq; }
		if ($lt20_n > 1) { $lowq = 1; return $lowq; }
	} else {
		if ($lt10_n > 1 || $lt20_n > 2) { $lowq = 1; return $lowq; }
	}

	return $lowq;
}

#
sub usage{
	my $usage = << "USAGE";
Program: clean-fastq
Version: 1.0
Contact: chgao\@genedenovo.com

Description:
  Filter Low Quality -> clean reads
  Trim 3\' adaptor
  Remove 5\' adaptor contaminants
  Remove too small and polyA
  Remove error tags
  Stats small RNA frequence
  Out Put -> clean tag

Usage: $0 [options] 2> log
Options:
  -i <file>  fastq format reads file
  -d <dir>   result directory, default="./";
  -a <str>   5\' adaptor sequence, default="GTTCAGAGTTCTACAGTCCGACGATC";
  -b <str>   3\' adaptor sequence, default="CTGTAGGCACCATCA";
  -m <int>   minimum small RNA length, default=18
  -p <int>   quality phread ,default=33
  -e <float> average single base error rate, default=0.01
  -c <int>   frequence cutoff, default=0
  -s         filter low quality reads strictly, default=loose
             strict: base number [1-10, qv<10] <= 0 & base number [1-10, qv<20] <= 1
             loose: base number [1-10, qv<10] <= 1 & base number [1-10, qv<20] <= 2

  -r         remove bad tags (with sequencing error)
  -z <str>   yes/[no] remove random adapter 4bp + tag + 4bp
  -h help

USAGE

	print $usage;
}
