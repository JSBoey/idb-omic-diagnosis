#!/bin/bash -e

for i in 16S ITS; do
  printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" \
    > 1.sequence_processing/${i}.manifest
  fqdir=$(pwd -P)/1.sequence_processing/${i}/0.rawdata
  for j in ${fqdir}/*.1.fq.gz; do
    sample=$(basename ${j} .1.fq.gz)
    printf "%s\t%s\t%s\n" ${sample} ${fqdir}/$(basename ${j}) ${fqdir}/${sample}.2.fq.gz \
      >> 1.sequence_processing/${i}.manifest
  done
done