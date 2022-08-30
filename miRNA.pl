#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: miRNA.pl
#
#        USAGE: ./miRNA.pl  
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
#      CREATED: 06/04/2020 10:31:32 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use FindBin qw/$RealBin $Script/;
#use File::Basename qw/basename dirname/;
use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;
use lib "$RealBin/../";
use lib "/Bio/Bin/pipeline//System_Programs/Small_pipe/v1.7/../../hsm_perl_lib/bin/";
use Support_Program;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 
            heatmap 
            enrich 
            pca 
            mirna_mature 
            exp 
            diff 
            target 
            length_distribution 
            align_stat 
            mirna_stat 
            mirna_cluster 
            mirna_edit 
            mirna_family 
            mirdeep2

    e.g: perl $0 <> <>

USAGE

die "$usage" if(@ARGV<1); 

my %opts = (run=>1,"cpu"=>10,"multi_type"=>"qsub-sge");
GetOptions(\%opts,"run=s","pipe_info=s","outdir=s","outfile=s","dis=s",
            "fa=s","pos_txt=s","genome=s","type=s","sp_class=s",
            "exist_target=s","target_mrna=s","shortname=s","taxid=s",
            "ref=s","clean_txt=s","annotation=s","sample_list=s",
            "org=s","mirbase=s","annot=s","max_precursors=s",
            "run_reactome=s","run_do=s",
            "norm=s","para=s","cps=s","multi_type=s",
            );
my $type = shift;
$0 = "$0 $type";

my %apps;
my $soft_conf = "$RealBin/../Small_pipe.soft.conf";
my $soft_path = "$RealBin/../";
&Support_Program::Support_Program_conf("$soft_path",\%apps,$soft_conf);
AddSoft(\%apps,"$RealBin");

InitOpts(run=>$opts{run},'multi_type'=>$opts{multi_type},max_job=>$opts{cpu},shell_dir=>$opts{outdir},
        qsub_sge=>"$apps{perl} $apps{qsub_sge}",
        qsub_pbs=>"$apps{perl} $apps{qsub_pbs}",
        );
my %pipe_info = ReadPipe($opts{pipe_info}) if(exists $opts{pipe_info});
my @mirnas = ReadmiRNA("$pipe_info{conf_file}/mirna.list") if(exists $opts{pipe_info});
MakeDir($opts{outdir}) if(exists $opts{outdir});

my %pipe = (
    heatmap => \&Heatmap,
    enrich => \&Enrich,
    pca => \&PCA,
    mirna_mature => \&miRNAMature,
    exp => \&Exp,
    diff => \&Diff,
    target => \&Target,
    length_distribution => \&LengthDistribution,
    align_stat => \&AlignStat,
    mirna_stat => \&miRNAStat,
    mirna_cluster => \&miRNACluster,
    mirna_edit => \&miRNAEdit,
    mirna_family => \&miRNAFamily,
    mirdeep2 => \&miRdeep2,
);
Error2("Unknown type [$type]") if(!exists $pipe{$type});
&{$pipe{$type}}(\%pipe_info,%opts);

sub Target{
    my $pipe_info_ref = shift;
    my %opts = @_;

    AddOpts('multi_type'=>"threads");

    my $dir = $opts{outdir};
    ## exist_mirna 
    my @targets;
    my $cmd;
    if(exists $pipe_info_ref->{result_exist_mirna_txt}){
        $cmd .= "$apps{perl} $apps{fish} -ff fasta $pipe_info_ref->{result_exist_mirna_txt} $opts{exist_target} | gzip >$dir/exist_target.aln.gz\n";
        push @targets,"$dir/exist_target.aln.gz";
    }

    ## known_mirna 
    $cmd .= TargetPrediction($pipe_info_ref,$dir,%opts,name=>"known");
    push @targets,"$dir/known_target.aln.gz";
    my @class1_targets = @targets;

    ## novel_mirna  
    $cmd .= TargetPrediction($pipe_info_ref,$dir,%opts,name=>"novel");
    push @targets,"$dir/novel_target.aln.gz";

#    $cmd = "echo 1\necho 2\necho 3\necho 4\necho 5 && exit 1\necho 6";
    ## run exist known novel together
    MultiExcu($cmd,$dir,"target",cpu=>3);

    ## class1_mirna 
    Excu2("cat @class1_targets >$dir/class1_target.aln.gz");
    
    ## all_mirna 
    Excu2("cat @targets >$dir/all_target.aln.gz");

    ## aln2table 
    Excu2("$apps{perl} $apps{aln2table} $dir/all_target.aln.gz $dir/all_mirna.target.aln.xls $opts{sp_class} $opts{annot}");

    ## 
    Excu2("cut -f1,2 $dir/all_mirna.target.aln.xls >$dir/all_mirna.target.xls");

    ## stat align 
    Excu2("$apps{perl} $apps{stat_target_aln} $pipe_info_ref->{result_all_exp} $dir/all_mirna.target.xls $dir/all_target.stat.xls");

    return;
}

sub TargetPrediction{
    my $pipe_info_ref = shift;
    my $dir = shift;
    my %opts = @_;
 
    my $name = $opts{name}; ## known novel
    my $mirna_fa_key = "result_$name\_mirna_fa";
    my $cmd = "";
    if($opts{sp_class} eq 'plant'){# plant
        $cmd .= "$apps{perl} $apps{plant_target_prediction} -m $pipe_info_ref->{$mirna_fa_key} -d $opts{target_mrna} -n $name -gz 1 -o $dir\n";
    }else{ #animal
        $cmd .= "$apps{perl} $apps{animal_target_prediction} -multi_type $opts{multi_type} $pipe_info_ref->{$mirna_fa_key} $opts{target_mrna} $opts{shortname} $opts{taxid} $dir/$name\_mirna_target $dir/$name\_target.aln.gz\n";
#        Excu2("ln -sf $name\_mirna_target/all.target.aln $dir/$name\_target.aln");
    }

#    return "$dir/$name\_target.aln";
    return "$cmd";
}

sub Diff{
    my $pipe_info_ref = shift;
    my %opts = @_;

    $opts{norm} //= "yes";
    $opts{para} //= "p-0.05-2";
    my $all_exp_file = "$pipe_info_ref->{result_all_exp}";
    my $conf_dir = "$pipe_info_ref->{conf_file}";
    MakeDir("$opts{outdir}/tmp");
    my @mirna_lists;
    for my $type(qw/exist_mirna known_mirna novel_mirna/){
        my $key = "result_$type\_txt";
        next if(!exists $pipe_info_ref->{$key});
        push @mirna_lists,$pipe_info_ref->{$key};
    }
    my $mirna_list = join(",",@mirna_lists);

    my $word = "sample diff";
    my $group = "-g $conf_dir/group.list";
    my $diff_list = "$conf_dir/group_diff.list";
    my ($pq,$pq_cut,$fc_cut) = split/-/,$opts{para}; ## p-0.05-2
    if($opts{type} eq 'sample_diff'){
        $word = "sample diff";
        $group = "";
        $diff_list = "$conf_dir/sample_diff.list";
    }
    Excu2("$apps{perl} $apps{edgeR} -f $all_exp_file  $group -c $diff_list -m yes -n $opts{norm} -o $opts{outdir}/tmp");
    Excu2("$apps{perl} $apps{edgeR_miRNA_split} -i $opts{outdir}/tmp -g $diff_list -l $mirna_list -p $pq -c $pq_cut -f $fc_cut -t \"$word\" -o $opts{outdir}");

    return;
}

sub Exp{
    my $pipe_info_ref = shift;
    my %opts = @_;

    Excu2("$apps{perl} $apps{get_sample_count_for_exp} $pipe_info_ref->{result_annotation_stat} exist_mirna,exist_mirna_edit,known_mirna,novel_mirna >$opts{outdir}/sample_count.txt");
    my $cmd1;
    my @types = qw/exist_mirna known_mirna novel_mirna/;
    my %exp = map{$_=>"none"} @types;
    for my $type(@types){
        my $key = "result_$type\_txt";
        next if(!exists $pipe_info_ref->{$key});
        Excu2("$apps{perl} $apps{sample_count2tpm} $pipe_info_ref->{$key} $opts{outdir}/sample_count.txt >$opts{outdir}/$type.exp_profile.xls");
        $exp{$type} = "$opts{outdir}/$type.exp_profile.xls";
    }
    if($exp{exist_mirna} ne 'none' || $exp{known_mirna} ne 'none'){
        Excu2("$apps{perl} $apps{combine_exp} $exp{exist_mirna} $exp{known_mirna} >$opts{outdir}/class1_mirna.exp_profile.xls");
    }
    Excu2("$apps{perl} $apps{combine_exp} $exp{exist_mirna} $exp{known_mirna} $exp{novel_mirna} >$opts{outdir}/all_mirna.exp_profile.xls");

    return;
}

sub miRNAMature{
    my $pipe_info_ref = shift;
    my %opts = @_;

    my %mireap = (
            'plant' =>  "$apps{perl} $apps{'mireap-plant'}",
            'animal' => "$apps{perl} $apps{'mireap-animal'}",
            'fungi' => "$apps{perl} $apps{'mireap-animal'}",
            'model_species' => "$apps{perl} $apps{'mireap-animal'}",
            );
    my $sp_class = $opts{sp_class};
    my $outdir = $opts{outdir};
    my $mireap = $mireap{$sp_class};
    MakeDir($outdir);

    Excu2("$apps{perl} $apps{glibc} \"$mireap -i $opts{fa} -m $opts{pos_txt} -r $opts{genome} -o $outdir -t $opts{type}\"");
    if($opts{type} eq 'known'){
        Excu2("$apps{perl} $apps{known_mirna_aln_to_table} $outdir/*.aln >$outdir/$opts{type}_mirna.table.xls");
    }else{
        Excu2("$apps{perl} $apps{novel_mirna_aln_to_table} $outdir/mireap-$opts{type}.aln $outdir/$opts{type}_mirna");
    }
    Excu2("$apps{perl} $apps{exist_known_mirna_create_rnastructure} $outdir/$opts{type}_mirna.table.xls $outdir/structure");

    return;
}

sub PCA{
    my $pipe_info_ref = shift;
    my %opts = @_;

    my $cmd1;
    my $group = "$pipe_info_ref->{conf_file}/group.list";
    for my $mi(@mirnas){
        $cmd1 .= "perl $apps{filter_low_exp_mirna_4pca} $pipe_info_ref->{expression}/$mi.exp_profile.xls 1 >$opts{outdir}/$mi.exp.filter.xls\n";
        $cmd1 .= "perl $apps{R_plot} PCA -scale yes -type PCA,3D -group $group -infile $opts{outdir}/$mi.exp.filter.xls -outpfx $opts{outdir}/$mi\n";
    }
    Excu2("$cmd1",split=>0);
#    print "$cmd1\n";

    return;
}

sub Heatmap{
    my $pipe_info_ref = shift;
    my %opts = @_;

    my $cmd1;
    for my $mi(@mirnas){
        $cmd1 .= "perl $apps{filter_low_exp_mirna_4pca} $pipe_info_ref->{expression}/$mi.exp_profile.xls 1 >$opts{outdir}/$mi.exp.filter.xls\n";
        $cmd1 .= "perl $apps{R_plot} pheatmap -infile $opts{outdir}/$mi.exp.filter.xls -outdir $opts{outdir}/ -outpfx $mi.heatmap\n";
    }
    Excu2($cmd1,split=>0);
#    print "$cmd1\n";

    return;
}

sub Enrich{
    my $pipe_info_ref = shift;
    my %opts = @_;

    my $cmd1;
    my $cmd2;

    $opts{run_do} //= 0;
    $opts{run_reactome} //= 0;

    my $sample_list = "$pipe_info_ref->{conf_file}/sample.list";
    my $ref = $opts{ref};
    my $target_dir = "$pipe_info_ref->{target_prediction}";
    my $all_aln = "$pipe_info_ref->{target_prediction}/all_target.aln.gz";
    $all_aln = $pipe_info_ref->{result_mirna_target} if(exists $pipe_info_ref->{result_mirna_target});
    for my $mi(@mirnas){
        MakeDir("$opts{outdir}/$mi");
    }

    my $enrich_cmd = "perl $apps{enrich} -ref $ref -type nodiff -format -run_do $opts{run_do} -run_reactome $opts{run_reactome}";
    if(-s $sample_list){
        $cmd1 .= "$apps{perl} $apps{aln2glist} $all_aln $pipe_info_ref->{mirna_list} $sample_list $pipe_info_ref->{result_all_exp} samples $opts{outdir}\n";
        for my $mi(@mirnas){
            my $dir = MakeDir("$opts{outdir}/$mi/samples");
#            $cmd1 .= "perl $apps{aln2glist} $sample_list $target_dir $dir\n";
            $cmd2 .= "$enrich_cmd -i $dir\n";
        }
    }

    my $sample_diff_list = "$pipe_info_ref->{conf_file}/sample_diff.list";
    if(-s $sample_diff_list){
        $cmd1 .= "$apps{perl} $apps{aln2glist} $all_aln $pipe_info_ref->{mirna_list} $sample_diff_list $pipe_info_ref->{sample_diff}/all_mirna diff $opts{outdir}\n";
        for my $mi(@mirnas){
            my $dir = MakeDir("$opts{outdir}/$mi/diff");
#            $cmd1 .= "perl $apps{diff_aln2glist} $sample_diff_list $pipe_info_ref->{sample_diff}/$mi $dir\n";
            $cmd2 .= "$enrich_cmd -i $dir\n";
        }
    }

    my $group_diff_list = "$pipe_info_ref->{conf_file}/group_diff.list";
    if(-s $group_diff_list){
        $cmd1 .= "$apps{perl} $apps{aln2glist} $all_aln $pipe_info_ref->{mirna_list} $group_diff_list $pipe_info_ref->{group_diff}/all_mirna group_diff $opts{outdir}\n";
        for my $mi(@mirnas){
            my $dir = MakeDir("$opts{outdir}/$mi/group_diff");
#            $cmd1 .= "perl $apps{diff_aln2glist} $group_diff_list $pipe_info_ref->{group_diff}/$mi $dir\n";
            $cmd2 .= "$enrich_cmd -i $dir\n";
        }
    }
    MultiExcu($cmd1,$opts{outdir},"glist");
    MultiExcu($cmd2,$opts{outdir},"enrich");
#    print "=>$cmd1\n";
#    print "=>$cmd2\n";

    return;
}

sub LengthDistribution{
    my $pipe_info_ref = shift;
    my %opts = @_;

    Excu2("$apps{perl} $apps{length_distribution} $opts{clean_txt} $opts{annotation} $opts{outdir}");
    open my $in_fh,"$opts{sample_list}" or die "$!:$opts{sample_list}\n";
    my $cmd;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        my $outpfx = "$opts{outdir}/$aa[0].length_distribution";
        my $cmd_info = qq(-xlab "Length(nt)" -ylab "Tag Count" -title "Length Distribution" -infile $outpfx.xls -outpfx $outpfx);
        if($opts{annotation} eq 'none'){
            $cmd .= qq($apps{perl} $apps{R_plot} bar -show_text no $cmd_info\n);
        }else{
            $cmd .= qq($apps{perl} $apps{R_plot} fill_bar -type count $cmd_info\n);
        }
    }
    close $in_fh;
    Excu2($cmd);

    return;
}

sub AlignStat{
    my $pipe_info_ref = shift;
    my %opts = @_;

    MakeDir($opts{outdir});
    my $outpfx = "$opts{outdir}/match_$opts{type}.stat";
    Excu2("$apps{perl} $apps{align_stat} $opts{clean_txt} $opts{annotation} $opts{type} $outpfx");
   
    if($opts{type} ne 'genome2'){
        my $cmd = qq($apps{perl} $apps{R_plot} fill_bar -type fill -xlab "Sample" -ylab "Percentage" -title "Tag Match $opts{type} Stat" -rmcol 1);
        Excu2("$cmd -infile $outpfx.xls -outpfx $outpfx");
        Excu2("$cmd -infile $outpfx.unique.xls -outpfx $outpfx.unique");
    }

    return;
}

sub miRNAStat{
    my $pipe_info_ref = shift;
    my %opts = @_;

    my @types = qw/exist_mirna known_mirna novel_mirna/;

    for my $mi(@types){
        my $key = "result_$mi\_txt";
        next if(!exists $pipe_info_ref->{$key});
        my $dir = MakeDir("$opts{outdir}/$mi");
        my $outpfx = "$dir/match_$mi.stat";
        Excu2("$apps{perl} $apps{align_stat} -mirna $pipe_info_ref->{$key} $opts{clean_txt} $opts{annotation} $mi $outpfx");
        Excu2("$apps{perl} $apps{plot_mirna_base} $pipe_info_ref->{$key} $dir $mi");
        open my $in_fh,"$opts{sample_list}" or die "$!:$opts{sample_list}\n";
        my $cmd;
        while(<$in_fh>){
            chomp;
            my @aa = split/\t/,$_;
            my $outpfx = "$dir/$_.$mi\_base_bias.stat";
            $cmd .= qq($apps{perl} $apps{R_plot} fill_bar -type fill -xlab "Position" -ylab "Percentage" -title "Nucleotide Bias at each Position" -infile $outpfx.xls -outpfx $outpfx\n);
            $outpfx = "$dir/$_.$mi\_first_base.stat";
            $cmd .= qq($apps{perl} $apps{R_plot} fill_bar -type fill -xlab "Length(nt)" -ylab "Percentage" -title "First Nuleotide Bias" -infile $outpfx.xls -outpfx $outpfx\n);
        }
        close $in_fh;
        Excu2($cmd);
    }

    if(exists $pipe_info_ref->{exist_mirna_edit}){
        AlignStat($pipe_info_ref,%opts,outdir=>"$opts{outdir}/exist_mirna_edit",type=>"exist_mirna_edit");    
    }
    

    return;
}

sub miRNACluster{
    my $pipe_info_ref = shift;
    my %opts = @_;

    ParaError2(\%opts,"pipe_info,outfile,mirbase","dis[10000]");
    $opts{dis} //= 10000;
    my @types = qw/exist_mirna known_mirna novel_mirna/;
    my @file_lists;
    for my $mi(@types){
        my $key = "result_$mi\_table";
        next if(!exists $pipe_info_ref->{$key});
        push @file_lists,$pipe_info_ref->{$key};
    }
    my $file_list = join(",",@file_lists);

    Excu2("$apps{perl} $apps{mirna_cluster} $file_list $opts{mirbase} $opts{dis} $opts{outfile}");

    return;
}

sub miRNAEdit{
    my $pipe_info_ref = shift;
    my %opts = @_;

    ParaError2(\%opts,"pipe_info,outdir");
    my $exist_mirna_txt = $pipe_info_ref->{result_exist_mirna_txt};
    my $clean_txt = $pipe_info_ref->{result_clean_txt};
    my $edit_aln = $pipe_info_ref->{result_exist_mirna_edit_aln};

    my $outfile = "$opts{outdir}/miRNA_SNP.xls";
    my $stat_pfx = "$opts{outdir}/miRNA_SNP_stat";
    Excu2("$apps{perl} $apps{mirna_edit} $edit_aln $exist_mirna_txt $clean_txt $outfile $stat_pfx.xls");
    Excu2("$apps{perl} $apps{R_plot} fill_bar -type count -title \"\" -xlab \"Sample\" -ylab \"Number of editing sites\" -infile $stat_pfx.xls -outpfx $stat_pfx");

    return;
}

sub miRNAFamily{
    my $pipe_info_ref = shift;
    my %opts = @_;

    ParaError2(\%opts,"pipe_info,outdir,org,mirbase");
    for my $mi(qw/exist_mirna known_mirna/){
        my $key = "result_$mi\_txt";
        next if(!exists $pipe_info_ref->{$key});
        my $outpfx = "$opts{outdir}/$mi.family";
        my $exist_mirna_txt = $pipe_info_ref->{$key};
        Excu2("$apps{perl} $apps{get_mirna_family} $mi $opts{org} $exist_mirna_txt $opts{mirbase} $outpfx.xls $outpfx\_detail.xls");
    }
    
    return;
}

sub miRdeep2{
    my $pipe_info_ref = shift;
    my %opts = @_;

    $opts{max_precursors} //= 50000;
    ParaError2(\%opts,"fa,pos_txt,outdir,genome","max_precursors");
    MakeDir($opts{outdir});

    Excu2("$apps{perl} $apps{pre_mirdeep_fa} $opts{fa} $opts{outdir}/tmp.fa $opts{pos_txt} $opts{outdir}/aln.txt"); 
    Excu2("$apps{perl} $apps{format_genome_4_mirdeep2} $opts{genome} >$opts{outdir}/genome.tmp.fa");
    $opts{genome} = "$opts{outdir}/genome.tmp.fa";
    Excu2("$apps{perl} $apps{miRDeep2} $opts{outdir}/tmp.fa $opts{genome} $opts{outdir}/aln.txt none none none $opts{outdir} -v -g $opts{max_precursors}");
    Excu2("$apps{perl} $apps{mirdeep_result} $opts{outdir}/output.mrd $opts{outdir}/result.csv  $opts{outdir}/pdfs $opts{outdir}");
    Excu2("rm -f $opts{outdir}/genome.tmp.fa");
    Excu2("rm -f $opts{outdir}/tmp.fa $opts{outdir}/aln.txt");
    
    return;
}


sub ReadPipe{
    my $list = shift;

    open DIR,$list or die $!;
    my %hash;
    while(<DIR>){
        chomp;
        my @aa = split/\t/,$_;
        $hash{$aa[0]} = $aa[1];
    }
    close DIR;
    
    my @mirna_lists;
    for my $type(qw/exist_mirna known_mirna novel_mirna/){
        my $key = "result_$type\_txt";
        next if(!exists $hash{$key});
        push @mirna_lists,$hash{$key};
    }
#    my $mirna_list = join(",",@mirna_lists);
    $hash{mirna_list} = join(",",@mirna_lists);

    return %hash;
}

sub AddSoft{
    my $apps_ref = shift;
    my $bin = shift;
    my $src = "$bin/src";
    $apps_ref->{aln2glist} = "$src/aln2glist.pl";
    $apps_ref->{diff_aln2glist} = "$src/diff_aln2glist.pl";
    $apps_ref->{filter_low_exp_mirna} = "$src/filter_low_exp_mirna.pl";
    $apps_ref->{filter_low_exp_mirna_4pca} = "$src/filter_low_exp_mirna_4pca.pl";
    $apps_ref->{known_mirna_aln_to_table} = "$src/known_mirna_aln_to_table.pl";
    $apps_ref->{novel_mirna_aln_to_table} = "$src/novel_mirna_aln_to_table.pl";
    $apps_ref->{get_sample_count_for_exp} = "$src/get_sample_count_for_exp.pl";
    $apps_ref->{sample_count2tpm} = "$src/sample_count2tpm.pl";
    $apps_ref->{stat_target_aln} = "$src/stat_target_aln.pl";
    $apps_ref->{combine_exp} = "$src/combine_exp.pl";
    $apps_ref->{align_stat} = "$src/align_stat.pl";
    $apps_ref->{format_genome_4_mirdeep2} = "$src/format_genome_4_mirdeep2.pl";
    $apps_ref->{mirdeep_result} = "$src/mirdeep_result.pl";
    $apps_ref->{pre_mirdeep_fa} = "$src/pre_mirdeep_fa.pl";
    $apps_ref->{get_mirna_family} = "$src/get_mirna_family.pl";
    $apps_ref->{mirna_edit} = "$src/mirna_edit.pl";
    $apps_ref->{mirna_cluster} = "$src/mirna_cluster.pl";
    $apps_ref->{length_distribution} = "$src/length_distribution.pl";
    $apps_ref->{plot_mirna_base} = "$src/plot_mirna_base.pl";
    $apps_ref->{glibc} = "$RealBin/../glibc-2.14/glibc.pl";
    $apps_ref->{edgeR} = "$RealBin/../edgeR/runedgeR.pl";
    $apps_ref->{miRDeep2} = "$RealBin/../mirdeep2_0_0_7/miRDeep2.pl";
    $apps_ref->{edgeR_miRNA_split} = "$RealBin/../edgeR/edgeR_miRNA_split.pl";
    $apps_ref->{animal_target_prediction} = "$RealBin/../target_prediction/animal_target_prediction.pl";
    $apps_ref->{plant_target_prediction} = "$RealBin/../target_prediction/plant/plant_target_prediction.pl";
    $apps_ref->{aln2table} = "$RealBin/../target_prediction/aln2table.pl";

    return;
}

sub ReadmiRNA{
    my $file = shift;
    
    my @out;
    open my $in_fh,"$file" or die "$!:$file";
    while(<$in_fh>){
        chomp;
        push @out,$_;
    }
    close $in_fh;

    return @out;
}
