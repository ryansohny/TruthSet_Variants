#!/bin/bash
# mhsohny
# 09-23-2024

#SBATCH --job-name=pickvar
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=15:00:00
#SBATCH --mem=5G
#SBATCH -o logs/slurm.%N.%j.out
#SBATCH -e logs/slurm.%N.%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mhsohny@uw.edu
#SBATCH --export=ALL

source activate base
conda activate mcg

dir=/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1

for sample in COLO829BL COLO829T_PassageA_DSA COLO829T_PassageB_DSA
do

:<<'END'
echo "Split Multialleles for ${sample}"
bcftools norm \
--write-index=tbi \
-m -any \
${dir}/${sample}/deepvariant/${sample}.deepvariant.vcf.gz \
-o ${dir}/${sample}/deepvariant/${sample}.deepvariant.split.vcf.gz

echo "Pick only the variant-type-of-interest for ${sample}"
bcftools view \
-v snps \
--write-index=tbi \
${dir}/${sample}/deepvariant/${sample}.deepvariant.split.vcf.gz \
-o ${dir}/${sample}/deepvariant/${sample}.deepvariant.split.snv.vcf.gz
END

zcat ${dir}/${sample}/deepvariant/${sample}.deepvariant.split.snv.vcf.gz | \
awk 'BEGIN {OFS="\t"} { if (/^#/) { print } else { $4 = substr($4, 1, 1); $5 = substr($5, 1, 1); print } }' | \
bgzip > ${dir}/${sample}/deepvariant/${sample}.deepvariant.split.snv.modified.vcf.gz

tabix -p vcf ${dir}/${sample}/deepvariant/${sample}.deepvariant.split.snv.modified.vcf.gz

#Check out the things below (check to see if these are present. If they are there, something went wrong):
#zcat ${sample}.deepvariant.split.snv.vcf.gz | awk '!/^#/ && length($4)>=2 && substr($4, 2) != substr($5, 2)' - | less -S

done
