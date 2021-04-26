#!/usr/bin/perl -w
use strict;
use Cwd qw(cwd getcwd);
#use PBS::Client;
#TopHat/2.1.0
#Bowtie2/2.2.6
#fastqc/0.11.3
#StringTie/1.3.3
#my $directory="/user/scratch/gent/gvo000/gvo00027/data/RNA_atlas/results";
my $directory="data/fastq";

#my $index="/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/transcriptome_data_38/known";
my $index="data/transcriptome_data_38/known";

#my $index2 = "/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/hg38_index/hg38";
my $index2 = "data/hg38_index/hg38";

my $ensembl_transcriptome="data/Ensemblv86.gtf"
opendir(DIR, $directory) or die $!;

my @fastq_folders=();
while (my $folder = readdir(DIR)) {
        push(@fastq_folders, $folder);
		
}

@fastq_folders=grep(!/^\./, @fastq_folders);
closedir(DIR);

my @commands=();
push (@commands, "#!/bin/sh");
push (@commands, "cd $directory");
foreach (@fastq_folders){
	chomp;
	my $fastq_folder=$_;
	print $fastq_folder."\t";
	opendir(DIR, "$directory/$fastq_folder") or die $!;
	(my $samplename) = $fastq_folder =~ /^([A-Za-z0-9_-]*)/;
	print $samplename."\t";
	push (@commands, "cd $fastq_folder");
	push (@commands, "gunzip *.fastq.gz");
    push (@commands, "cat *R1*.fastq > ".$samplename."_1.fastq");
    push (@commands, "cat *R2*.fastq > ".$samplename."_2.fastq");
	#push (@commands, "rm *R1*.fastq");
	#push (@commands, "rm *R2*.fastq");
	#push (@commands, "module load fastqc/0.11.3");
	push (@commands, "fastqc ".$samplename."_1.fastq");
	push (@commands, "fastqc ".$samplename."_2.fastq");
	#push (@commands, "module load TopHat/2.1.0-intel-2015b");
    #push (@commands, "module load Bowtie2/2.2.6-intel-2015b");
    push (@commands, "tophat -p 8 --no-coverage-search -o ".$fastq_folder."_thout --library-type=fr-firststrand --transcriptome-index=".$index." ".$index2." ".$samplename."_1.fastq ".$samplename."_2.fastq");
    push (@commands, "gzip *.fastq");
    #push (@commands, "module load StringTie/1.3.3-intel-2017a");
    push (@commands, "stringtie $fastq_folder"."_thout/accepted_hits.bam -p 8 -G $ensembl_transcriptome -o $fastq_folder"."_st_assembly.gtf -C $fastq_folder"."_cov_ref.gtf -A $samplename"."_gen_abund.tab");
 	my $cdir = getcwd;
	push (@commands, "echo $cdir/$fastq_folder"."_st_assembly.gtf >> ../../list_gtf_files.txt");
    push (@commands, "cd ..");

	
}

open(my $fh, '>', 'qc_align_assemble.sh');
print $fh join("\n",@commands);
close $fh;
system("chmod 577 ./qc_align_assemble.sh")
system("./qc_align_assemble.sh")

