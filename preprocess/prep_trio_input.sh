#!/usr/bin/bash
#Sample script to generate a joint called trio phased vcf from a proband and 2 parent GVCF and BAMS
#useful for testing duoNovo's output against trio based de novo calling


if [[ -z $(which bcftools 2> /dev/null) ]]; then
   echo Dependency not met:  bcftools not found in path.  Please install bcftools or ensure bcftools is in executable PATH.
   exit
fi

if [[ -z $(which glnexus_cli 2> /dev/null) ]]; then
   echo Dependency not met:  glnexus_cli not found in path.  Please install glnexus or ensure glnexus_cli is in executable PATH.
   exit
fi

if [[ -z $(which hiphase 2> /dev/null) ]]; then
   echo Dependency not met:  hiphase not found in path.  Please install hiphase or ensure hiphase is in executable PATH.
   exit
fi


if [[ -z "$8" ]]; then #allows a task ID to be passed to command line if not run as an array
  echo Usage: prep_trio_input.sh PROBAND_GVCF PROBAND_BAM PARENT_1_GVCF PARENT_1_BAM PARENT_2_GVCF PARENT_2_BAM REFERENCE_FASTA OUTPUT_VCF [THREADS:default nproc]
  echo Requires: bcftools, glnexus_cli, hiphase
  echo Inputs:
  echo PROBAND_GVCF: g.vcf file from the proband \(can be g.vcf or g.vcf.gz\)
  echo PROBAND_BAM: bam \(or cram\) file from the proband
  echo PARENT1_GVCF: g.vcf file from the parent1 \(can be g.vcg or g.vcf.gz\)
  echo PARENT1_BAM: bam \(or cram\) file from the parent1
  echo PARENT2_GVCF: g.vcf file from the parent1 \(can be g.vcg or g.vcf.gz\)
  echo PARENT2_BAM: bam \(or cram\) file from the parent1
  echo REFERENCE_FASTA: reference fasta file for the aligned sequences
  echo OUTPUT_VCF: vcf file to be created, a joint called \(glnexus DeepVariant_unfiltered\) duo vcf with both samples with phasing added from hiphase
  echo optional  THREADS: number of thread to use for glnexus and hiphase.  Default is is nproc - $(nproc)
  exit
fi

PROBANDGVCF="$1"
PROBANDBAM="$2"
PARENT1GVCF="$3"
PARENT1BAM="$4"
PARENT2GVCF="$5"
PARENT2BAM="$6"
REF="$7"
OUTPUT="$8"
THREADS=${9:-$(nproc)}


if [ ! -f "$PROBANDGVCF" ]; then
      echo "File not found! : $PROBANDGVCF"
      exit
fi

if [ ! -f "$PROBANDBAM" ]; then
      echo "File not found! : $PROBANDBAM"
      exit
fi

if [ ! -f "$PARENT1GVCF" ]; then
      echo "File not found! : $PARENT1GVCF"
      exit
fi

if [ ! -f "$PARENT1BAM" ]; then
      echo "File not found! : $PARENT1BAM"
      exit
fi

if [ ! -f "$PARENT2GVCF" ]; then
      echo "File not found! : $PARENT2GVCF"
      exit
fi

if [ ! -f "$PARENT2BAM" ]; then
      echo "File not found! : $PARENT2BAM"
      exit
fi

#detecting Proband ID and Parent ID from gVCF header
PROBAND=$(bcftools view -h $PROBANDGVCF | tail -n1 | cut -f10)
PARENT1=$(bcftools view -h $PARENT1GVCF  | tail -n1 | cut -f10)
PARENT2=$(bcftools view -h $PARENT2GVCF  | tail -n1 | cut -f10)

if [[ -z "$PROBAND" ]]; then
      echo FAILED TO identify PROBAND ID from GVCF: $PROBANDGVCF
      exit
fi

if [[ -z "$PARENT1" ]]; then
      echo FAILED TO identify PARENT1 ID from GVCF: $PARENT1GVCF
      exit
fi

if [[ -z "$PARENT2" ]]; then
      echo FAILED TO identify PARENT2 ID from GVCF: $PARENT2GVCF
      exit
fi


TMPDIR=$(dirname $OUTPUT)/tmp_trio_$PROBAND\_$PARENT1\_$PARENT2\_$$

mkdir -p $TMPDIR
glnexus_cli --config DeepVariant_unfiltered --dir $TMPDIR/tmp_GLNEXUS_$PROBAND\_$PARENT1\_$PARENT2_$$ --threads $THREADS $PROBANDGVCF $PARENT1GVCF $PARENT2GVCF | \
                         bcftools view --write-index -Oz -o $TMPDIR/unphased_tmpTRIO_$PROBAND\_$PARENT1\_$PARENT2\_$$.vcf.gz

# bcftools index /scratch/sberger/pmgrc_lr_data/inputs/$PROBAND/duo_mother.vcf.gz  # NOT needed in newer versions of bcftools which accept --write-index parameter

hiphase --bam $PROBANDBAM --bam $PARENT1BAM  --bam $PARENT2BAM \
            --sample-name $PROBAND \
            --sample-name $PARENT1 \
            --sample-name $PARENT2 \
            --threads $THREADS \
            --vcf $TMPDIR/unphased_tmpTRIO_$PROBAND\_$PARENT1\_$PARENT2\_$$.vcf.gz \
            --output-vcf $OUTPUT \
            --reference $REF

rm -rf $TMPDIR

