#!/bin/bash
#
# This script takes VCF files and generates files that can be used
# by the find_intersecting_snps.py script. 
# The script takes an INPUT directory and an OUTPUT directory.
#
# The INPUT directory is expected to contain files ending with .vcf.gz
# and to contain the name of the chromosome (like chr22 or chr1).
#
# OUTPUT files are named $CHR.snps.txt.gz (where $CHR is the name of the
# chromosome). OUTPUT files contain <position> <allele1> <allele2> on
# each line.
#
while getopts i:o:c: opt; do
  case $opt in
  i)
      INPUT_DIR=$OPTARG
      ;;
  o)
      OUTPUT_DIR=$OPTARG
      ;;
  c)
      COLS= $OPTARG
      ;;
esac

if [[ -z ${INPUT_DIR+x} ]] || [[ -z ${OUTPUT_DIR+x} ]]
then
	echo "Input or Output Directory Missing"
fi

if [[ ! -z ${COLS+x} ]]
then
	COLUMS='$2,$4,$5'$COLS
else
	COLUMS='$2,$4,$5'
fi

if [ ! $INPUT_DIR ]; then
    echo "usage: extract_vcf_snps.sh <input_dir> <output_dir>" >&2
    exit 2
fi

if [ ! $OUTPUT_DIR ]; then
    echo "usage: extract_vcf_snps.sh <input_dir> <output_dir>" >&2
    exit 2
fi

mkdir -p $OUTPUT_DIR

vcf_files=$INPUT_DIR/*vcf
for FILE in $vcf_files; do
     echo $FILE >&2
     CHR=`echo $FILE | sed -n 's/^.*\(chr[0-9A-Z]*\).*.vcf$/\1/p'`
     echo $CHR >&2
     OUTPUT_FILE=$OUTPUT_DIR/$CHR.snps.txt
     egrep -v "^#" $FILE | awk "{print $COLUMS }" > $OUTPUT_FILE
done

