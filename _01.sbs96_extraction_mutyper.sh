#!/bin/bash
# mhsohny
# 10-20-2024

#SBATCH --job-name=mutyper
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1:00:00
#SBATCH --mem=40G
#SBATCH -o logs/slurm.%N.%j.out
#SBATCH -e logs/slurm.%N.%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mhsohny@uw.edu
#SBATCH --export=ALL

source activate base
conda activate mutyper

DSA=/mmfs1/gscratch/stergachislab/assemblies/DSA_COLO829BL_v3.0.0.fasta
dir=/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1/Mutational_Spectrum

for sample in COLO829T_PassageB_DSA COLO829T_PassageA_DSA; do
  mutyper variants \
  --k 3 \
  "${DSA}" \
  "${dir}/${sample}.deepvariant.split.snv.modified.final.vcf" \
  | mutyper spectra - \
  > "${dir}/01.SBS/${sample}.SBS96"

  mutyper variants \
  --k 3 \
  "${DSA}" \
  "${dir}/${sample}.deepvariant.split.snv.modified.final.tba.vcf" \
  | mutyper spectra - \
  > "${dir}/01.SBS/${sample}.tba.SBS96"

  if [ "$sample" == "COLO829T_PassageB_DSA" ]; then
    mutyper variants \
    --k 3 \
    ${DSA} \
    ${dir}/${sample}.deepvariant.split.snv.modified.final.onlytb.vcf \
    | mutyper spectra - \
    > ${dir}/01.SBS/${sample}_onlytb.SBS96

  elif [ "$sample" == "COLO829T_PassageA_DSA" ]; then
    mutyper variants \
    --k 3 \
    ${DSA} \
    ${dir}/${sample}.deepvariant.split.snv.modified.final.onlyta.vcf \
    | mutyper spectra - \
    > ${dir}/01.SBS/${sample}_onlyta.SBS96
  fi
done