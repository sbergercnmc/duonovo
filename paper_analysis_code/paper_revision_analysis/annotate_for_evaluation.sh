#!/bin/bash -i
source ~/.bashrc

DN_VCF=$1

#detect duoNovo proband ID from vcf header
PROBAND_ID=$( bcftools view -h $DN_VCF | grep duoNovoPARAM | grep proband_column_identifier | cut -d= -f4 | cut -d$">" -f1 )

if [[ -z "$PROBAND_ID" ]]; then
   echo NOT A VALID DUO NOVO VCF
   exit
fi

#detect sequenced parent ID from vcf column headers
PARENT_ID=$( bcftools view -h $DN_VCF | sed 's/\t/\n/g' | tail -n2 | grep -v "$PROBAND_ID" )

#confirm existence of trio called vcf
TRIO_VCFGZ=$(dirname "$DN_VCF")/$PROBAND_ID.PFM.vcf.gz

if [[ ! -f "$TRIO_VCFGZ" ]]; then
    echo Can not find TRIO vcf
    exit
fi

echo $TRIO_VCFGZ
OTHER_PARENT_ID=$(bcftools view -h $TRIO_VCFGZ | sed 's/\t/\n/g' | tail -n3 | grep -v $PROBAND_ID | grep -v $PARENT_ID)

if [[ -z "$PARENT_ID" ]]; then
   echo NOT A VALID DUO VCF
   exit
fi

if [[ -z "$OTHER_PARENT_ID" ]]; then
   echo NOT A VALID TRIO VCF
   exit
fi


echo PROBAND:$PROBAND_ID
echo PARENT:$PARENT_ID
echo OTHER_PARENT:$OTHER_PARENT_ID


AVOUT=$(readlink -f $DN_VCF | sed 's/.duonovo.vcf$/.duonovo.annovar/')

echo $(date) :: Starting Annovar Annotation to $AVOUT

 if [[ ! -f $AVOUT.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.vcf.gz ]] && [[ ! -f $AVOUT.addedParent.dnm2.vcf.gz ]]; then
   table_annovar.pl $DN_VCF $AVDB -buildver hg38 -out $AVOUT -protocol refGeneWithVer,bed,bed,gnomad41_genome -operation g,r,r,f -bedfile hg38_cpg.txt,hg38_problematic.txt -remove -polish -vcfinput -nastring . -onetranscript  -thread 20 > $AVOUT.avout 2> $AVOUT.averr
 fi

echo $(date) :: Fixing Annovar VCF Annotation to $AVOUT

   if [[ -f $AVOUT.hg38_multianno.vcf ]] && [[ ! -f $AVOUT.clean.vcf ]]; then
      cat $AVOUT.hg38_multianno.vcf |sed 's/Name\\x3dNA/True/g' | sed 's/bed2/problematic_region/g' | sed 's/bed/cpg/g' > $AVOUT.clean.vcf
   fi

echo $(date) :: Indexing VCF to $AVOUT.clean.vcf.gz
if [[ ! -f $AVOUT.clean.vcf.gz ]] && [[ -f $AVOUT.clean.vcf ]]; then
bcftools view -W -Oz -o $AVOUT.clean.vcf.gz $AVOUT.clean.vcf
fi

echo $(date) :: Extracting other parent

if [[ ! -f $AVOUT.with_other_parent.vcf ]] && [[ -f $AVOUT.clean.vcf.gz ]]  && [[ -f $TRIO_VCFGZ ]] && [[ ! -f $AVOUT.addParent.dnm2.vcf.gz ]]; then
bcftools view -s $OTHER_PARENT_ID -W -Oz -o $AVOUT.otherparent_tmp.$OTHER_PARENT_ID.vcf.gz $TRIO_VCFGZ
echo $(date) :: Merging other parent
bcftools merge -m none $AVOUT.clean.vcf.gz $AVOUT.otherparent_tmp.$OTHER_PARENT_ID.vcf.gz | \
bcftools view -s $PROBAND_ID,$PARENT_ID,$OTHER_PARENT_ID -i 'INFO/phasing_proband="0|1" || INFO/phasing_proband="1|0"' -o $AVOUT.with_other_parent.vcf
fi

echo $(date) :: Annotating de novos

if [[ -f $AVOUT.with_other_parent.vcf ]] && [[ ! -f $AVOUT.addedParent.dnm2.vcf.gz ]]; then
numX=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $5}'  | sed 's/0/1/')X
father=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $3}')
mother=$(cat $GREGOR/data_model/somalier_output.samples.tsv | awk -v P=$PROBAND_ID '$2==P {print $4}')
#echo bcftools +trio-dnm2 -p $numX:$PROBAND_ID,$father,$mother --use-NAIVE -Oz -o  $AVOUT.addedParent.dnm2.vcf.gz $AVOUT.with_other_parent.vcf
bcftools +trio-dnm2 -p $numX:$PROBAND_ID,$father,$mother --use-NAIVE -Oz -o  $AVOUT.addedParent.dnm2.vcf.gz $AVOUT.with_other_parent.vcf
fi

echo $(date) :: Writing data table

if [[ -f $AVOUT.addedParent.dnm2.vcf.gz ]]; then
 rm -f $AVOUT.avinput
 rm -f $AVOUT.hg38_multianno.* $AVOUT.clean.vcf* $AVOUT.otherparent_tmp.$OTHER_PARENT_ID.vcf.gz*
 rm -f $AVOUT.with_other_parent.vcf
bcftools query -H -f "%CHROM\t%POS\t%REF\t%ALT\t%Func.refGeneWithVer\t[%DNM\t%VA\t]%duoNovo_classification\t%n_de_novo_right_orientation_same_PS\t%n_de_novo_left_orientation_same_PS\t%QC_fail_step\t%cpg\t%problematic_region\t%gnomad41_genome_AF\t[\t%GT\t%GQ\t%DP]"  $AVOUT.addedParent.dnm2.vcf.gz > $AVOUT.datatable.txt
fi
