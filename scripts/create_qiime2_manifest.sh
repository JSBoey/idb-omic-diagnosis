#!/bin/bash -e

for i in 16S ITS; do

    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" \
        > 2.asv/${i}.manifest

    fqdir=$PWD/0.fq/2.primer_trimmed

    for j in ${fqdir}/*${i}_pt.1.fq.gz; do

        printf "%s\t%s\t%s\n" \
            $(basename ${j} _${i}_pt.1.fq.gz) \
            ${j} \
            ${j/\.1\./\.2\.} \
            >> 2.asv/${i}.manifest
            
    done

done