#!/usr/bin/bash
#example script for preparing duoNovo input vcf
#Inputs: g.vcf files from proband and parent such as outputted by deepvariant, bam or cram files for proband and parent and reference fasta
#Uses bcftools, glnexus_cli, and hiphase

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


if [[ -z "$6" ]]; then #allows a task ID to be passed to command line if not run as an array
  echo Usage: prep_input.sh PROBAND_GVCF PROBAND_BAM PARENT_GVCF PARENT_BAM REFERENCE_FASTA OUTPUT_VCF [THREADS:default nproc]
  echo Requires: bcftools, glnexus_cli, hiphase
  echo Inputs:
  echo PROBAND_GVCF: g.vcf file from the proband (can be g.vcf or g.vcf.gz)
  echo PROBAND_BAM: bam (or cram) file from the proband
  echo PARENT_GVCF: g.vcf file from the parent (can be g.vcg or g.vcf.gz)
  echo PARENT_BAM: bam (or cram) file from the proband
  echo REFERENCE_FASTA: reference fasta file for the aligned sequences
  echo OUTPUT_VCF: vcf file to be created, a joint called (glnexus DeepVariant_unfiltered) duo vcf with both samples with phasing added from hiphase
  echo optional  THREADS: number of thread to use for glnexus and hiphase.  Default is is nproc - $(nproc)
  echo optional  TMP_DIRECTORY: Folder to store glnexus temporary files.  This this be deleted before and after the script.
else

  PROBANDGVCF="$1"
  PROBANDBAM="$2"
  PARENTGVCF="$3"
  PARENTBAM="$4"
  REF="$5"
  OUTPUT="$6"
  THREADS=${7:-$(nproc)}
  TMPDIR=$(dirname $OUTPUT)


   if [ ! -f "$PROBANDGVCF" ]; then
      echo "File not found! : $PROBANDGVCF"
      exit
   fi

   if [ ! -f "$PROBANDBAM" ]; then
      echo "File not found! : $PROBANDBAM"
      exit
   fi

   if [ ! -f "$PARENTGVCF" ]; then
      echo "File not found! : $PARENTGVCF"
      exit
   fi

   if [ ! -f "$PARENTBAM" ]; then
      echo "File not found! : $PARENTBAM"
      exit
   fi

  #detecting Proband ID and Parent ID from gVCF header
  PROBAND=$(bcftools view -h $PROBANDGVCF | tail -n1 | cut -f10)
   PARENT=$(bcftools view -h $PARENTGVCF  | tail -n1 | cut -f10)


  if [[ -z "$PROBAND" ]]; then
      echo FAILED TO identify PROBAND ID from GVCF: $PROBANDGVCF
      exit
  fi

  if [[ -z "$PARENT" ]]; then
      echo FAILED TO identify PARENT ID from GVCF: $PARENTGVCF
      exit
  fi


   rm -rf $TMPDIR/tmp_$PROBAND\_tmp_$$

   mkdir -p $TMPDIR
   glnexus_cli --config DeepVariant_unfiltered --dir $TMPDIR/tmp_$PROBAND\_tmp_$$ --threads $THREADS $PROBANDGVCF $PARENTGVCF | \
                         bcftools view --write-index -Oz -o $OUTPUT\_tmpDUO_$$\_.vcf.gz

   rm -rf $TMPDIR/tmp_$PROBAND\_tmp_$$
    # bcftools index /scratch/sberger/pmgrc_lr_data/inputs/$PROBAND/duo_mother.vcf.gz  NOT needed in newer versions of bcftools which accept --write-index parameter

    hiphase --bam $PROBANDBAM --bam $PARENTBAM  \
             --sample-name $PROBAND \
             --sample-name $PARENT \
             --threads $THREADS
             --vcf  $OUTPUT\_tmpDUO_$$\_.vcf.gz \
             --output-vcf $OUTPUT \
             --reference $REF

   rm -rf $OUTPUT\_tmpDUO_$$\_.vcf

fi
