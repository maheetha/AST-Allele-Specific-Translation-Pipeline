#!/bin/bash

module load samtools
module load python/3.4.5
module load bowtie/2.2.9
module load tophat

# Set these environment vars to point to
# your local installation of WASP
while getopts hw:d:s:v:n:i:pa:b: opt; do
  case $opt in
  w)
      WASP=$OPTARG
      ;;
  d)
      DATA_DIR=$OPTARG
      ;;
  s)
      SNP_DIR=$OPTARG
      ;;
  v)
      VCF_DIR=$OPTARG
      ;;
  n)
      NAME=$OPTARG
      ;;
  i)
      INDEX=$OPTARG
      ;;
  p)
      OPTION="paired_ended"
      ;;
  a)
      FASTQ1=$OPTARG
      ;;
  b) 
      FASTQ2=$OPTARG
      ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  h)
      echo "This program has many mandatory and optional arguments:
      
	Mandatory Arguments:

	-h : this help message

	-w : the WASP directory with all of the necessary WASP scripts
	
	-i : the full path to the reference genome

	-a : at least one fastq file 

	-n : please provide a distinct, recognizable, base name for all your bam files that come out of this program.

	Conditionally Dependent Mandatory Arguments:

	-v : a directory to the VCF files. Note that this is mandatory if you do not provide a snp directory with valid files

	-s : a directory to the SNP files. We will create one in your current directory if you don't provide this. But this means you have to provide a VCF directory. If you have valid snps files, please put them in an snp directory. 

	Optional Arguments:

	-p : signals whether it's paired ended. no argument need follow this.

	-b : if you have selected paired ended, a second fastq file is necessary.
	
	"
	exit 1
  esac
done

shift $((OPTIND - 1))


if [ -z ${NAME+x} ]
then
	echo "please provide a name"
	exit 1
fi

if [ -z ${DATA_DIR+x} ]
then
	DATA_DIR=$PWD
fi

echo $DATA_DIR

if [ -z ${VCF_DIR+x} ] && [ -z ${SNP_DIR+x} ]
then
        echo "Please provide either a VCF directory or SNP directory"
        exit 1
fi

if [ -z ${SNP_DIR+x} ]
then
        echo "you haven't provided an SNP directory, so we are assuming you want to extract SNPs"
        mkdir $DATA_DIR/snp_dir
        SNP_DIR=$DATA_DIR/snp_dir
fi

if [ -z ${INDEX+x} ]
then
	echo "You need a reference genome"
	exit 1
fi

if [ -z ${FASTQ1+x} ]
then
	echo "You need to input at least one FASTQ file"
	exit 1
fi

if [ ! -z ${OPTION+x} ] && [ -z ${FASTQ2+x} ]
then
	echo "You have selected paired_end but haven't provided the second fastq file"
	exit 1
fi

if [ -z ${OPTION+x} ] && [ ! -z ${FASTQ2+x} ]
then
        echo "You haven't selected paired_end but you have provided a second fastq file. What do you plan to do with the second one? Just know that the second fastq file will not be used."
       
fi


if [ ! -z ${OPTION+x} ]
then
	OPTION= "single_end"
fi


# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2

#Extract SNPs from VCF files (optional)
$WASP/extract_vcf_snps.sh $VCF_DIR $SNP_DIR 

# Map reads using bowtie2 (or another mapping tool of your choice)
if [[ $OPTION == "paired_end" ]]
then
	tophat $INDEX $DATA_DIR/$FASTQ1 $DATA_DIR/$FASTQ2 
else
	tophat $INDEX $DATA_DIR/$FASTQ1
fi

# Make an output directory
mkdir $DATA_DIR/output_directory

OUT_DIR=$DATA_DIR/output_directory

## Move tophat alignments into that directory
mv $DATA_DIR/tophat_out/accepted_hits.bam $OUT_DIR/$1.bam

## Sort the Bam
samtools sort -o $OUT_DIR/$1.sorted.bam $OUT_DIR/$1.bam

## Edit Bam Names
bash edit_name.sh $OUT_DIR/$1.sorted.bam

## Find which reads to re-align
if [[ $OPTION == "paired_end" ]]
then
	python $WASP/find_intersecting_snps.py -p --is_sorted $OUT_DIR/$1.sorted.bam --snp_dir $SNP_DIR 
	tophat $INDEX $OUT_DIR/$1.sorted.remap.fq1.gz $OUT_DIR/$1.sorted.remap.fq2.gz
else
	python $WASP/find_intersecting_snps.py --is_sorted $OUT_DIR/$1.sorted.bam --snp_dir $SNP_DIR
	tophat $INDEX $OUT_DIR/$1.sorted.remap.fq.gz 
fi

## re-alignment
#tophat $INDEX $OUT_DIR/$1.sorted.remap.fq1.gz $OUT_DIR/$1.sorted.remap.fq2.gz 

## moving out of tophat directory 
mv $OUT_DIR/tophat_out/accepted_hits.bam $OUT_DIR/$1.remapped.bam

## filter the remapped reads for those that don't map uniquely
python $WASP/filter_remapped_reads.py $OUT_DIR/$1.to.remap.bam  \
       $OUT_DIR/$1.remapped.bam $OUT_DIR/$1.remap.keep.bam  

## merge all reads, sort, and index
samtools merge $OUT_DIR/$1.keep.merged.bam \
        $OUT_DIR/$1.keep.bam $OUT_DIR/$1.remap.keep.bam

samtools sort -o $OUT_DIR/$1.keep.merged.sorted.bam $OUT_DIR/$1.keep.merged.bam 

samtools index $OUT_DIR/$1.keep.merged.sorted.bam

## filter reads for duplicates
if [[ $OPTION == "paired_end" ]]
then
	python $WASP/rmdup_pe.py $OUT_DIR/$1.keep.merged.sorted.bam $OUT_DIR/$1.keep.rmdup.merged.sorted.bam
else
	python $WASP/rmdup.py $OUT_DIR/$1.keep.merged.sorted.bam $OUT_DIR/$1.keep.rmdup.merged.sorted.bam
fi

## sort and pileup the final bam. Split the pileup files based on chromosome ID
samtools sort -o $OUT_DIR/$1.keep.rmdup.merged.resorted.bam $OUT_DIR/$1.keep.rmdup.merged.sorted.bam 

samtools mpileup -Q 0 -q 2 $OUT_DIR/$1.keep.rmdup.merged.resorted.bam > $OUT_DIR/$1.pileup

awk -F'\t' '{print > $1".pileup"}' $DATA_DIR/$1.pileup



