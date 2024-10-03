#!/bin/bash

# Check if both arguments are provided
if [ -z "$1" ] || [ -z "$2" ]; then
	echo "Usage   : ./$0 <input cram> <region>"
	echo "example : ./$0 ./COLO829BL/COLO829BL_DSA_resetmapq.cram haplotype1-0000001:94816619"
	exit 1
fi

inputbcram=$1
region=$2
outputbcram=$(dirname "${inputbcram}")/ROI/$(echo "${region}" | sed 's/:/_/g')_$(basename "${inputbcram}")

echo "${outputbcram}"
echo "${inputbcram}"

mkdir -p $(dirname "${inputbcram}")/ROI

samtools view \
-@ 10 \
-O CRAM \
--output-fmt-option embed_ref=1 \
--write-index \
-o ${outputbcram} \
${inputbcram} \
${region}
