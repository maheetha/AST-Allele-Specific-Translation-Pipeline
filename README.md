# AST-Allele-Specific-Translation-Pipeline


Pipeline for Allele Specific Translation Using Ribosome Profiling


This directory contains scripts and test files that utilizes the WASP program recently developed (refer to their README for more information) to develop a sensitive pipeline that can take fastq files and develop starter position files (SPFs). These SPFs can then be used in a bootstrapping algorithm to detect allele-specific translation. This README will the general steps in the AST pipeline workflow, as well as a description of the starter position file generator python script.

REQUIRED PROGRAMS:

PYTHON VERSION 3.x Version, preferably 3.4.5
TOPHAT 2.x version, preferably 2.1.1
SAMTOOLS 1.3.1 or later



General Steps in the Allele Specific Translation Workflow 

PHASE 1: 

Extract SNPs from VCFs

The extract_vcf_snps.sh is borrowed from the WASP program and takes in a directory of clearly labeled vcfs and imputes SNPs from them. The columns have been edited to take in input beyond the 2nd, 4th, and 5th column. The extra columns do not affect the other steps of the WASP program, but will be crucial for output of the starter position files. 


PHASE 2: 

Alignment using tophat. Paired ended versus single ended options are available.

PHASE 3: Only if paired ended

Editing bam in file, only if paired ended.

The editing.sh script takes in the bam and edits the names of the reads to make pairs the same. This is only an issue if the reads are paired-ended. 

PHASE 4:

The next step is to find the reads that overlap heterozygous snps and realign them. We use WASP’s find_intersecting_snps.py. A list of changes that have been made to the script itself are appended at the end of this readme. 

PHASE 5:

Filtering and finalizing bam files. Filtering_reads.py, rmdup.py, and rmdup_pe.py are borrowed from the WASP pipeline, and any modifications made have been appended to the end of this ReadMe.  


USAGE and OPTIONS OF AST_mapping_workflow

This program has many mandatory and optional arguments:
      
        Mandatory Arguments:

        -w : the WASP directory with all of the necessary WASP scripts
        
        -i : the full path to the reference genome

        -a : at least one fastq file 

        -n : please provide a distinct, recognizable, base name for all your bam files that come out of this program.

        Conditionally Dependent Mandatory Arguments:

        -v : a directory to the VCF files. Note that this is mandatory if you do not provide a snp directory with valid files

        -s : a directory to the SNP files. We will create one in your current directory if you don't provide this. But this means you have to provide a VCF directory. If you have valid snps files, please put them in an snp directory. 

        Optional Arguments:
        
       -h : this help message

        -p : signals whether you are doing paired-ended or single-ended alignment. no argument need follow this.

        -b : if you have selected paired ended, a second fastq file is necessary.
        

Example Usage:
bash AST_mapping_workflow.sh -w /wasp/directory -n basename -i /path/to/hg19 -s /path/to/snp/directory -p -a fastq1 -b fastq2


Working with Starter Position Files 

A Starter Position File (SPF) is derived from pileup files. This file includes a set of positions with the number of reads whose 5’ end begins at each position, separated by haplotype. These files will be used in the bootstrapping algorithm later. 

The starter_position_file.py takes in a pileup file and a file with snps and haplotype information for ONE chromosome in ONE sample.
 
Both the pileup file and the snps file should be for ONE chromosome. 

The pileup file should be in this format:

CHR	POS	SNP	COUNT	ALLELE	QUALITY
chr1	12345	N	1	C	I
chr1	12346	N	1	A	I
chr1	12347	N	1	A	F
chr1	12348	N	1	C	H
chr1	12349	N	1	C	@
chr1	12350	N	1	T	G

The format of the pileup file above should be the result of the samtools command in the last two lines of the AST_workflow script. The positions in the pileup up file should be in order, and samtools should automatically do this. If not, please sort the pileup file. In addition, please make sure that you have gotten rid of reads that do not match the quality you need. This should also be taken care of by the samtools mpileup command, where you can specify mapping and minimum base pair quality. 

The snp files you’ve created from extract_vcf_snps.sh should look like this:

10177	A	AC	1|0
10352	T	TA	0|1
13273	G	C	1|0
14464	A	T	0|1
14930	A	G	0|1
14933	G	A	0|1
15211	T	G	0|1
15274	A	G,T	1|2
15820	G	T	0|1
15903	G	GC	1|0
54490	G	A	0|1
54712	T	TTTTC	0|1

The fourth column is essential to help us keep track of the haplotype. The starter_positions_file.py program will take alleles that are SNPs, but no insertions and deletions. This program will also take multiple alleles, up to three (the max). 

Usage

NOTE: please use full paths for everything that is not in the directory where you're running the program. For directories, do not put a "/" at the end

python starter_position_file.py --pileup_file pileup --snp_file snp --read_length 75 --unit chr1 --min_reads 20 --output_dir out_dir

OUTPUT

POS	REF	ALT	REF_COUNT	ALT_1_COUNT	ALT_2_COUNT	ALT_3_COUNT	HAP	
22418721	G	A	1	2	0	0	0|1
22418720	G	A	0	0	0	0	0|1
22418719	G	A	0	0	0	0	0|1
22418718	G	A	1	0	0	0	0|1
22418717	G	A	0	0	0	0	0|1
22418716	G	A	0	0	0	0	0|1
22418715	G	A	1	0	0	0	0|1
22418714	G	A	2	1	0	0	0|1
22418713	G	A	0	2	0	0	0|1
22418712	G	A	3	9	0	0	0|1

The output looks like this. With position, reference and alternative allele(s), and the haplotype. 

Testing this pipeline

I've included "accepted_hits_chr1.sorted.bam" and chr1.snps.txt as rather large test files. With these test files, you can start with find_intersecting_snps.py and work through the pipeline, and you can even try starter_position_files.py.

Modifications to WASP files

Due to some difficulties with python versions, all of the WASP files have been modified to work well with python version 3.4.5. You can always download the original files from the WASP github page, and compare with our modified files. All the modifications have been appended to the end of this readme.

Melted Bedgraph for Mappability Filtering

We have provided a bed file formatted file that provides a score for each position on each chromosome. It should be one if the position is uniquely mappable. You can use this file to filter any positions from your starter position files. 

Changes made from original WASP scripts to Work with Python 3.4.5

find_intersecting_snps.py:

1.	line 510 added
	
	for y in range(len(ref_alleles)):
        	item = ref_alleles[y].decode("utf-8")
        	ref_alleles[y] = item

    	for z in range(len(alt_alleles)):
        	item = alt_alleles[z].decode("utf-8")
        	alt_alleles[z] = item

	
2.	line 510 of original has become line 518
	line 518 and 519 have been modified

	 ref_read = read_seq[:idx] + ref_alleles[i].decode("utf-8") + read_seq[idx+1:]
   	 alt_read = read_seq[:idx] + alt_alleles[i].decode("utf-8") + read_seq[idx+1:]

3. 	line 559 & 561 in original has been modified to
	fastq_file1.write(bytes("@%s\n%s\n+%s\n%s\n" %
                          (name, pair[0], name, orig_read1.qual) ,'UTF-8'))


       	 fastq_file2.write(bytes("@%s\n%s\n+%s\n%s\n" %
                          (name, rev_seq, name, orig_read2.qual), 'UTF-8'))

4. 	after line 818 in the original we added
	ref_alleles = ref_alleles.decode('UTF-8')
	
	after line 819 in original we added
	alt_alleles = alt.alleles.decode('UTF-8')


snp_table.py:

1. 	around lines 180 these two lines were added:
	allele2 = allele2.decode('UTF-8')
       	 allele1 = allele1.decode('UTF-8')


util.py

1. 	line 14 was modified to
	DNA_COMP = str.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")

extract_vcfs_snp.sh

Variable $COLUMS added to allow haplotype information to be included in extracted snps. 

