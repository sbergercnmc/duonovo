#!/bin/bash -i
source ~/.bashrc

TRIO_VCF=$1

#detect proband ID from filename
PROBAND_ID=$( echo $TRIO_VCF | rev | cut -d/ -f1 | rev | cut -d. -f1 )

if [[ -z "$PROBAND_ID" ]]; then
   echo NOT A VALID TRIO
   exit
fi

#confirm existence of trio called vcf
if [[ ! -f "$TRIO_VCF" ]]; then
    echo Can not find TRIO vcf
    exit
fi

echo PROBAND:$PROBAND_ID

if [[ ! -f $AVOUT.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.vcf.gz ]] && [[ ! -f $AVOUT.dnm2.vcf.gz ]]; then
bcftools view $TRIO_VCF > $(readlink -f $TRIO_VCF | sed 's/.PFM.vcf.gz$/.PFM.vcf/')
fi

AVOUT=$(readlink -f $TRIO_VCF | sed 's/.PFM.vcf.gz$/.PFM.annovar/')

echo $(date) :: Starting Annovar Annotation to $AVOUT

 if [[ ! -f $AVOUT.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.vcf.gz ]] && [[ ! -f $AVOUT.dnm2.vcf.gz ]]; then
   table_annovar.pl $(readlink -f $TRIO_VCF | sed 's/.PFM.vcf.gz$/.PFM.vcf/') $AVDB -buildver hg38 -out $AVOUT -protocol refGeneWithVer,bed,bed,gnomad41_genome -operation g,r,r,f -bedfile hg38_cpg.txt,hg38_problematic.txt -remove -polish -vcfinput -nastring . -onetranscript  -thread 20 > $AVOUT.avout 2> $AVOUT.averr
 fi

echo $(date) :: Fixing Annovar VCF Annotation to $AVOUT

   if [[ -f $AVOUT.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.vcf ]] ; then
      cat $AVOUT.hg38_multianno.vcf |sed 's/Name\\x3dNA/True/g' | sed 's/bed2/problematic_region/g' | sed 's/bed/cpg/g' > $AVOUT.clean.vcf
   fi


numX=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $5}'  | sed 's/0/1/')X
father=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $3}')
mother=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $4}')


echo $(date) :: Indexing VCF and reordering to $AVOUT.clean.vcf.gz
if [[ ! -f $AVOUT.clean.vcf.gz ]] && [[ -f $AVOUT.clean.vcf ]]; then
bcftools view -s $PROBAND_ID,$father,$mother -W -Oz -o $AVOUT.clean.vcf.gz $AVOUT.clean.vcf
fi

echo $(date) :: Annotating de novos

if [[ -f $AVOUT.clean.vcf.gz ]] && [[ ! -f  $AVOUT.dnm2.vcf.gz ]]; then
#echo bcftools +trio-dnm2 -p $numX:$PROBAND_ID,$father,$mother --use-NAIVE -Oz -o  $AVOUT.addedParent.dnm2.vcf.gz $AVOUT.with_other_parent.vcf
bcftools +trio-dnm2 -p $numX:$PROBAND_ID,$father,$mother --use-NAIVE -Oz -o  $AVOUT.dnm2.vcf.gz $AVOUT.clean.vcf.gz
fi

