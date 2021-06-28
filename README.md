# GenePy-1.4
Latest GenePy version using Ensembl VEP annotation, gnomAD random forest flag. Incorporation of tri-allelic variants *in progress*

GenePy v1.4 implements the following improvements from v1.3
* VCF annotation of allele frequencies and CADD is now completed through Ensembl-VEP, which is maintained and updated more consistently than ANNOVAR.
* CADD v.1.6 (also known as CADD-Splice), which improves the scoring of splicing variants in comparison to previous release CADD v.1.5
* Implementation of gnomAD v.3.1.1 Random Forest (RF) flag as an additional variant quality filter.

## GenePy Requirements
* A (multi)sample VCF file (can accept compressed vcf.gz)
* List of genes to generate GenePy scores (gene.list)
* List of variant types to include (variant_type.list)
* Ensembl-VEP installed
* Ensembl-VEP CADD Plugin
* GnomAD v.2.1.1 (exomes) sites vcf
* CADD v.1.6 Installed, with tsv files: whole_genome_SNVs and gnomad genomes indels
* Optional: ANNOVAR with gnomAD RF flag annotation files
* Vcftools
* Python 2.7.x

## Pre-processing for GenePy Scores

Prior to running the GenePy script, we need to annotate our VCF file, and filter accordingly to create the appropriate input file (ALL_genepy.meta)

### Bi-allelic variants only

We want to filter to include only bi-allelic variants

```
# Keep bi-allelic variants
vcftools --gzvcf GENOTYPED_ALL.vcf.gz --min-alleles 2 --max-alleles 2 --recode --remove-filtered-all --out FINAL
```
### Annotation

To run Ensembl-VEP, remember to softlink the location of the VEP software to your home directory

```
ln -s /software/location/VEP /your/home/directory/.vep
```
Annotate with VEP and ANNOVAR (optional)

```
# The VEP annotation script will perform the annotation of chromosomes in parallel to speed-up this process 
grep "^#" FINAL.recode.vcf > vcf_header
grep -v "^#" FINAL.recoded.vcf > FINAL_noheader.recode.vcf

sbatch --array=1-22 vep.sh FINAL_noheader.recode.vcf
sbatch vep_x.sh FINAL_noheader.recode.vcf

# Reassemble the annotated vcf
seq 1 22 > chr
grep "^#" ANNO_CHR1.vcf > vcf_header
while read l; do grep -v "^#" ANNO_CHR${l}.vcf >> temp_vcf; done < chr
cat vcf_header temp_vcf > FINAL_ANNOTATION.vcf

# Annovar annotation
sbatch annotation.sh
```
Our FINAL_ANNOTATION.vcf requires some manipulation to convert it into the correct input for GenePy. In addition, VEP annotation leaves some variants without a CADD score, so these need to be annotated online (https://cadd.gs.washington.edu/score), and the scores inserted into our annotated vcf. 

