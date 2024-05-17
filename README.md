# Information

## Data source

Project directory: `uoo03475`

- `Data from Azenta/R1_16S/Rawdata` Raw data for 16S amplicons
- `Data from Azenta/R2_ITS/Rawdata` Raw data for 18S amplicons
- `db` SILVA and UNITE databases
- `nobackup_space` -> `nobackup` mirror of project

## Primers

| Gene | Forward | Reverse | Region | Estimated size (bp) |
| --- | --- | --- | --- | --- |
| 16S | CCTACGGRRBGCASCAGKVRVGAAT | GGACTACNVGGGTWTCTAATCC | V3-V4 | ~ 600 |
| ITS | GTGAATCATCGARTC | TCCTCCGCTTATTGAT | NA | ~ 400 |

There does not seem to be a citable source for the primer sequences. 
[Ford et al. (2023)](https://newzealandecology.org/nzje/3545.pdf) has attributed
the sequences as proprietary to GENEWIZ by Azenta Life Sciences.

## Adapters

P7 (read 1) - AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
P5 (read 2) - AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# 0. Set up

## 1. Directories

In nobackup space

```bash
cd /nesi/nobackup/uoo03475

mkdir -p boey/{db,scripts,0.fq}
```

In project space

```bash
cd /nesi/project/uoo03475

mkdir -p db
```

## 2. Databases

Working directory: `/nesi/project/uoo03475/db`

Links to databases are in `db/db_links.txt`

```bash
while IFS= read -r line; do
    wget "$line"
done <db_links.txt
```

Move files tor respective directories

```bash
mkdir -p {SILVA_NR99_138.1,UNITE_10.0,Greengenes2_2022.10}

mv SILVA* SILVA_NR99_138.1
mv 2022.10* Greengenes2_2022.10
mv db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz UNITE_10.0
```

Check integrity of SILVA

```bash
cd SILVA_NR99_138.1
md5sum -c SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz.md5
cd ../
```

Rename UNITE archive

```bash
cd UNITE_10.0
mv db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz sh_qiime_release_04.04.2024.tar.gz
cd ../
```

Copy databases into nobackup space

```bash
cp -r * /nesi/nobackup/uoo03475/boey/db
```

Extract archives

```bash
cd /nesi/nobackup/uoo03475/boey/db

cd Greengenes2_2022.10
unzip -p \
    2022.10.backbone.tax.qza \
    c16a953c-f24d-4d14-927c-40d90ced395e/data/taxonomy.tsv \
    > Greengenes2_taxonomy.tsv

unzip -p \
    2022.10.backbone.full-length.fna.qza \
    a53d9300-5c5c-4774-a2e8-a5e23904f1ae/data/dna-sequences.fasta \
    > Greengenes2_full_length.fna

cd ../UNITE_10.0
tar -xvzf sh_qiime_release_04.04.2024.tar.gz

cd ../SILVA_NR99_138.1/
gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
```

## 3. Copy files

File integrity check

```bash
cd /nesi/project/uoo03475/Data\ from\ Azenta/

for i in R*/Rawdata; do
  cd ${i}
  md5sum -c md5.md5
  cd ../../
done
```

All OK!

Copy and append gene name as suffix

```bash
cd /nesi/nobackup/uoo03475/boey/0.fq

for i in 16S ITS; do
    mkdir -p ${i}
    cp --verbose \
        /nesi/project/uoo03475/Data\ from\ Azenta/*${i}/Rawdata/*.fq.gz \
        ${i}
done

for i in */*.gz; do
    dn=$(dirname ${i})
    mv ${i} ./"$(echo ${i} | sed "s/_L1/_${dn}_L1/g")"
done

mkdir -p 0.raw
mv {16,IT}S/* 0.raw/
rm -rf 16S/ ITS/
```

# 1. Quality control

Check raw sequences

```bash
cd /nesi/nobackup/uoo03475/boey

mkdir -p 1.qc
mkdir -p {0.fq,1.qc}/{1.adapter,2.primer}_trimmed

module load \
    FastQC/0.12.1 \
    MultiQC/1.13-gimkl-2022a-Python-3.10.5

fastqc --outdir 1.qc/0.raw --threads 6 0.fq/0.raw/*.fq.gz
multiqc --interactive --title 0.raw_fq --outdir 1.qc 1.qc/0.raw/
```

Quality problems (3' end, < Q30 after 210bp):

- UC-2_ITS_L1.1 
- UC-8_ITS_L1.1

Quality problems (3' end, < Q30 after 200bp):

- UC-12_ITS_L1.1
- UC-1_ITS_L1.1
- UC-5_ITS_L1.1
- UC-6_ITS_L1.1
- UC-9_ITS_L1.1
- control10_ITS_L1.1
- control12_ITS_L1.1
- control5_ITS_L1.1
- control7_ITS_L1.1

Huh... quality problems are confined to read 1. That's a little unexpected.

Moderate adapter read through:

- UC-4_ITS_L1.1
- UC-4_ITS_L1.2
- control3_ITS_L1.1
- control3_ITS_L1.2

Serious adapter read through:

- UC-10_ITS_L1.1
- UC-10_ITS_L1.2
- UC-3_ITS_L1.1
- UC-3_ITS_L1.2
- control6_ITS_L1.1
- control6_ITS_L1.2
- control9_ITS_L1.1
- control9_ITS_L1.2

ITS data has a bunch of adapter read through and poly G and poly A.

Those with adapter problems don't overlap with those with quality problems.

Trimming non-biological sequences

`scripts/read_trimming.sh`

```bash
#!/bin/bash -e

module purge
module load cutadapt/4.4-gimkl-2022a-Python-3.11.3

# Adapter trimming
# Also trim poly-G tails
echo "Trimming adapters"

A1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
A2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'

for r1 in 0.fq/0.raw/*_L1.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _L1.1.fq.gz)
    od=0.fq/1.adapter_trimmed

    cutadapt \
        -j 8 \
        --nextseq-trim=20 \
        --json 1.qc/1.adapter_trimmed/${bn}_at.cutadapt.json \
        -a ${A1} \
        -A ${A2} \
        -o ${od}/${bn}_at.1.fq.gz \
        -p ${od}/${bn}_at.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/1.adapter_trimmed/${bn}_at.cutadapt.log
done


# Primer trimming (16S)
echo "Trimming primers: 16S"

F16S='^CCTACGGRRBGCASCAGKVRVGAAT'
R16S='^GGACTACNVGGGTWTCTAATCC'

for r1 in 0.fq/1.adapter_trimmed/*16S_at.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _at.1.fq.gz)
    od=0.fq/2.primer_trimmed

    cutadapt \
        -j 8 --minimum-length 1 \
        --json 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.json \
        -g ${F16S} \
        -G ${R16S} \
        -o ${od}/${bn}_pt.1.fq.gz \
        -p ${od}/${bn}_pt.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.log
done


# Primer trimming (ITS)
echo "Trimming primers: ITS"

FITS='^GTGAATCATCGARTC'
RITS='^TCCTCCGCTTATTGAT'

for r1 in 0.fq/1.adapter_trimmed/*ITS_at.1.fq.gz; do
    r2=${r1/\.1\./\.2\.}
    bn=$(basename ${r1} _at.1.fq.gz)
    od=0.fq/2.primer_trimmed

    cutadapt \
        -j 8 --minimum-length 1 \
        --json 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.json \
        -g ${FITS} \
        -G ${RITS} \
        -o ${od}/${bn}_pt.1.fq.gz \
        -p ${od}/${bn}_pt.2.fq.gz \
        ${r1} ${r2} \
        > 1.qc/2.primer_trimmed/${bn}_pt.cutadapt.log
done
```

Check adapter trimming stats

```bash
module purge
module load \
    FastQC/0.12.1 \
    MultiQC/1.13-gimkl-2022a-Python-3.10.5

for i in 1.adapter_trimmed 2.primer_trimmed; do
    fastqc \
        --threads 8 \
        --outdir 1.qc/${i} \
        0.fq/${i}/*
    
    multiqc \
        --interactive \
        --title ${i} \
        --outdir 1.qc \
        1.qc/${i}
done
```

Reads seem okay. Some ITS samples were marked with poly-A tails in read 2. Not 
sure if this is normal but will continue forward.

# 2. Construct ASV table

Here, QIIME2's implementation of DADA2 will be used.

## 1. Import data into QIIME2

```bash
mkdir -p 2.asv
```

Create sample manifest for import.

`scripts/create_qiime2_manifest.sh`

```bash
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
```

Import

```bash
module purge
module load QIIME2/2024.2-shotgun

for i in 2.asv/*.manifest; do
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${i} \
        --output-path 2.asv/$(basename ${i} .manifest)_demux.qza \
        --input-format PairedEndFastqManifestPhred33V2
done
```

Check quality for trimming

```bash
for i in 2.asv/*.qza; do
    qiime demux summarize \
        --i-data ${i} \
        --o-visualization ${i/qza/qzv}
done
```

Doesn't look like 16S needs trimming. ITS may benefit from a truncation at 212.

NeSI doesn't have `dada2` plugin for the `2024.2-shotgun` version. Reverting to 
`2023.5` for all downstream processes.

`scripts/dada2_16S.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=dada2_paired_16S
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise 16S"

qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs 16S_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads 8 \
    --o-table 16S_feature_table.qza \
    --o-representative-sequences 16S_rep_seq.qza \
    --o-denoising-stats 16S_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data 16S_rep_seq.qza \
    --o-visualization 16S_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file 16S_denoise_stats.qza \
    --o-visualization 16S_denoise_stats.qzv
```

Job ID: 46294759

`scripts/dada2_ITS.sh`

```bash
#!/bin/bash -e
#SBATCH --job-name=dada2_paired_ITS
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise ITS"

qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs ITS_demux.qza \
    --p-trunc-len-f 212 \
    --p-trunc-len-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ITS_feature_table.qza \
    --o-representative-sequences ITS_rep_seq.qza \
    --o-denoising-stats ITS_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data ITS_rep_seq.qza \
    --o-visualization ITS_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file ITS_denoise_stats.qza \
    --o-visualization ITS_denoise_stats.qzv
```

Job ID: 46295010

# 2a. Dealing with ITS

The ITS sequences did very poorly on the above method (16-45% retained). There is a different processing method developed at USDA using HMM models called ITSxpress to identify ITS regions prior to DADA2.

NOTE: NeSI's QIIME2 does not have the latest version of ITSxpress. Will try using version 2 from [GitHub](https://github.com/USDA-ARS-GBRU/itsxpress).

## 0. Install ITSxpress

```bash
mkdir -p conda_env/

module purge
module load Mamba/23.1.0-1

mamba create -p conda_env/itsxpress -c bioconda -c conda-forge itsxpress

# Initialise conda environment in .bashrc
conda init
```

## 1. Trim ITS FASTQ

ITSxpress depends on identification of conserved regions within ITS. Some primer sequences overlap with those HMM models and trimming them off (as above) can lead to problems(see [here](https://github.com/USDA-ARS-GBRU/itsxpress/issues/40) and [here](https://github.com/USDA-ARS-GBRU/itsxpress/issues/24)). Therefore, the inputs must be adapter trimmed sequences prior to primer trimming.

`scripts/itsxpress.sh`

```bash
#!/bin/bash -e
#SBATCH --job-name=ITSxpress
#SBATCH --account=uoo03475
#SBATCH --time=30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --array=0-25
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load Mamba/23.1.0-1

a=(0.fq/1.adapter_trimmed/*ITS_at.1.fq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/.1./.2.}
b=$(basename ${r1} _at.1.fq.gz)

o1=0.fq/2.ITSxpress_trimmed/${b}_it.1.fq.gz
o2=${o1/.1./.2.}

id=1.0
taxa=Fungi
region=ITS2
log=1.qc/2.ITSxpress_trimmed/${b}.itsxpress.log

mkdir -p tmp/$SLURM_JOB_ID

source ~/.bashrc
conda activate conda_env/itsxpress

echo "ITSxpress trimming sample ${b}"

itsxpress \
    --fastq ${r1} \
    --fastq2 ${r2} \
    --outfile ${o1} \
    --outfile2 ${o2} \
    --region ${region} \
    --taxa ${taxa} \
    --cluster_id ${id} \
    --log ${log} \
    --tempdir tmp/$SLURM_JOB_ID \
    --threads $SLURM_CPUS_PER_TASK

mamba deactivate

```

JOBID: 46335253

## 2. DADA2 ASV construction

Import ITSxpress trimmed sequences into QIIME2

```bash
module purge
module load QIIME2/2023.5

# Create manifest
printf "%s\t%s\t%s\n" \
    "sample-id" \
    "forward-absolute-filepath" \
    "reverse-absolute-filepath" \
    > 2.asv/ITS_ITSxpress.manifest

for r1 in 0.fq/2.ITSxpress_trimmed/*.1.fq.gz; do
    r2=${r1/.1./.2.}
    b=$(basename ${r1} _ITS_it.1.fq.gz)

    printf "%s\t%s\t%s\n" \
        ${b} \
        '$PWD'/${r1} \
        '$PWD'/${r2} \
        >> 2.asv/ITS_ITSxpress.manifest
done

# Import manifest
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path 2.asv/ITS_ITSxpress.manifest \
    --output-path 2.asv/ITS_ITSxpress_demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

# Check sequence quality
qiime demux summarize \
    --i-data 2.asv/ITS_ITSxpress_demux.qza \
    --o-visualization 2.asv/ITS_ITSxpress_demux.qzv
```

There does not seem to be a need to truncate sequences for DADA2 (ref [here](https://forum.qiime2.org/t/dada2-good-quality-data-cut-or-not/27874)). The ITSxpress software has done some quality filtering which helps a lot.

`scripts/dada2_ITSxpress.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=dada2_ITSxpress
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise ITSxpress trimmed ITS reads"

b=ITS_ITSxpress

qiime dada2 denoise-single \
    --verbose \
    --i-demultiplexed-seqs ${b}_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ${b}_feature_table.qza \
    --o-representative-sequences ${b}_rep_seq.qza \
    --o-denoising-stats ${b}_denoise_stats.qza

qiime feature-table tabulate-seqs \
    --i-data ${b}_rep_seq.qza \
    --o-visualization ${b}_rep_seq.qzv

qiime metadata tabulate \
    --m-input-file ${b}_denoise_stats.qza \
    --o-visualization ${b}_denoise_stats.qzv
```

JOBID: 46336663

For some samples, a substantial majority of the reads are lost via merging. Use DADA2's single-ended denoise pipeline.

# 3. Single-ended denoising

Comparing results for forward read only denoising to recover more reads. This is a particularly bad problem with ITS reads as there is too much loss during merging for some samples.

`scripts/dada2_single.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=dada2_single
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

cd $PWD/2.asv
echo "Denoise single-end forward reads"

a=(16S ITS_ITSxpress)
b=${a[$SLURM_ARRAY_TASK_ID]}

qiime dada2 denoise-single \
    --verbose \
    --i-demultiplexed-seqs ${b}_demux.qza \
    --p-trunc-len 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --o-table ${b}_feature_table.se.qza \
    --o-representative-sequences ${b}_rep_seq.se.qza \
    --o-denoising-stats ${b}_denoise_stats.se.qza

qiime feature-table tabulate-seqs \
    --i-data ${b}_rep_seq.se.qza \
    --o-visualization ${b}_rep_seq.se.qzv

qiime metadata tabulate \
    --m-input-file ${b}_denoise_stats.se.qza \
    --o-visualization ${b}_denoise_stats.se.qzv
```

JOBID: 46341318

## Conclusion

Use merged 16S and single-end forward ITS. 


# 3. Feature classifier

Due to the non-standard primers used for this study, feature classifiers must be
trained from sequences.

NOTE: For UNITE ITS, do not trim references to primer sites! Train the 
classifier on full reference sequences using the `developer` sequences provided
with the QIIME release! See [here](https://docs.qiime2.org/2024.2/tutorials/feature-classifier/#classification-of-fungal-its-sequences)

## 1. Import database

QIIME2 has a helpful RESCRIPt plugin for obtaining pre-formatted SILVA objects.

```bash
module purge
module load QIIME2/2023.5

qiime rescript get-silva-data \
    --verbose \
    --p-version 138.1 \
    --p-target SSURef_NR99 \
    --p-include-species-labels \
    --o-silva-sequences db/SILVA_NR99_138.1/SILVA_NR99_138.1_rna_refseq.qza \
    --o-silva-taxonomy db/SILVA_NR99_138.1/SILVA_NR99_138.1_taxonomy.qza

# Reverse transcribe reference sequences
qiime rescript reverse-transcribe \
    --i-rna-sequences db/SILVA_NR99_138.1/SILVA_NR99_138.1_rna_refseq.qza \
    --o-dna-sequences db/SILVA_NR99_138.1/SILVA_NR99_138.1_refseq.qza
```

Import UNITE developer sequences clustered at 99% sequences similarity.

```bash
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path db/UNITE_10.0/developer/sh_refs_qiime_ver10_99_04.04.2024_dev.fasta \
    --output-path db/UNITE_10.0/UNITE_dev99_10.0_refseq.qza

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --input-path db/UNITE_10.0/developer/sh_taxonomy_qiime_ver10_99_04.04.2024_dev.txt \
    --output-path db/UNITE_10.0/UNITE_dev99_10.0_taxonomy.qza
```

Rename Greengenes2 database (simplifies array declarations)

```bash
mv db/Greengenes2_2022.10/2022.10.backbone.full-length.fna.qza db/Greengenes2_2022.10/Greengenes2_2022.10_refseq.qza
mv db/Greengenes2_2022.10/2022.10.backbone.tax.qza db/Greengenes2_2022.10/Greengenes2_2022.10_taxonomy.qza
```

## 2. Curate SILVA database

SILVA requires curation as per [this document](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/)

`scripts/curate_silva.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=curate_silva
#SBATCH --account=uoo03475
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

silva=db/SILVA_NR99_138.1/SILVA_NR99_138.1_refseq.qza
taxa=db/SILVA_NR99_138.1/SILVA_NR99_138.1_taxonomy.qza

# Remove low-quality sequences
cull=${silva/_refseq/_cull_refseq}

qiime rescript cull-seqs \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-sequences ${silva} \
    --o-clean-sequences ${cull}

# Remove short sequences
filt=${cull/_cull/_filt}
discard=${cull/_cull/_discard}

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences ${cull} \
    --i-taxonomy ${taxa} \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs ${filt} \
    --o-discarded-seqs ${discard}

# Dereplicate
derep_seq=${filt/_filt/_derep}
derep_tax=${taxa/_taxonomy/_derep_taxonomy}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${filt}  \
    --i-taxa ${taxa} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_seq} \
    --o-dereplicated-taxa ${derep_tax}

# In-silico PCR
segment=${derep_seq/_derep/_segment}
F16S='CCTACGGRRBGCASCAGKVRVGAAT'
R16S='GGACTACNVGGGTWTCTAATCC'

qiime feature-classifier extract-reads \
    --i-sequences ${derep_seq} \
    --p-f-primer ${F16S} \
    --p-r-primer ${R16S} \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --p-read-orientation 'forward' \
    --o-reads ${segment}

# Dereplicate segments
derep_segment_seq=${segment/_segment/_derep_segment}
derep_segment_tax=${derep_tax/_derep/_derep_segment}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${segment}  \
    --i-taxa ${derep_tax} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_segment_seq} \
    --o-dereplicated-taxa ${derep_segment_tax}
```

JOBID: 46359074

## 3. Curate Greengenes2

This is a shortened curation that involves only in-silico PCR and a dereplication step.

`scripts/curate_greengenes2.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=curate_gg2
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

gg=db/Greengenes2_2022.10/Greengenes2_2022.10_refseq.qza
tax=db/Greengenes2_2022.10/Greengenes2_2022.10_taxonomy.qza

F16S='CCTACGGRRBGCASCAGKVRVGAAT'
R16S='GGACTACNVGGGTWTCTAATCC'

# In-silico PCR
segment=${gg/_refseq/_segment_refseq}

qiime feature-classifier extract-reads \
    --i-sequences ${gg} \
    --p-f-primer ${F16S} \
    --p-r-primer ${R16S} \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --p-read-orientation 'forward' \
    --o-reads ${segment}

# Dereplicate segments
derep_segment_seq=${segment/_segment/_derep_segment}
derep_segment_tax=${tax/_taxonomy/_derep_segment_taxonomy}

qiime rescript dereplicate \
    --verbose \
    --p-threads $SLURM_CPUS_PER_TASK \
    --i-sequences ${segment}  \
    --i-taxa ${tax} \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${derep_segment_seq} \
    --o-dereplicated-taxa ${derep_segment_tax}
```
JOBID: 46361387

## 4. Train classifiers

Using Naive Bayes.

`scripts/train_classifier.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=train_classifier
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=milan
#SBATCH --array=0-2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(SILVA_NR99_138.1_derep_segment \
   Greengenes2_2022.10_derep_segment \
   UNITE_dev99_10.0)
b=${a[$SLURM_ARRAY_TASK_ID]}
db=$(echo ${b%%_derep_segment} | sed 's/dev99_//g')

seq=db/${db}/${b}_refseq.qza
tax=db/${db}/${b}_taxonomy.qza
cls=db/${db}_classifier.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --p-verbose \
    --i-reference-reads ${seq} \
    --i-reference-taxonomy ${tax} \
    --o-classifier ${cls}
```

JOBID: 46361899

# 4. Phylogeny

```bash
mdkir -p 3.phylogeny
```

A tree is required for UniFrac analyses (will use `rbiom` for this).

```bash
#!/bin/bash -e
#SBATCH --job-name=phylogeny
#SBATCH --account=uoo03475
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(16S_rep_seq.qza \
   ITS_ITSxpress_rep_seq.se.qza)
b=${a[$SLURM_ARRAY_TASK_ID]}
gene=3.phylogeny/${b%%_*}

qiime phylogeny align-to-tree-mafft-fasttree \
    --verbose \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --i-sequences 2.asv/${b} \
    --o-alignment ${gene}.alignment.qza \
    --o-masked-alignment ${gene}.masked_alignment.qza \
    --o-tree ${gene}.unrooted_tree.qza \
    --o-rooted-tree ${gene}.rooted_tree.qza
```

JOBID: 46361961

# 5. Assign taxonomy

```bash
mkdir -p 4.taxonomy
```

`scripts/16S_taxonomy.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=16S_taxonomy
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load QIIME2/2023.5

a=(Greengenes2_2022.10 SILVA_NR99_138.1)
cls=db/${a[$SLURM_ARRAY_TASK_ID]}_classifier.qza
db=${a[$SLURM_ARRAY_TASK_ID]%%_*}

reads=16S_rep_seq.qza

qiime feature-classifier classify-sklearn \
    --verbose \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-reads 2.asv/${reads} \
    --i-classifier ${cls} \
    --o-classification 4.taxonomy/${reads%%_*}.${db}.qza
```

JOBID: 46363815

`scripts/ITS_taxonomy.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=ITS_taxonomy
#SBATCH --account=uoo03475
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err

module purge
module load QIIME2/2023.5

cls=db/UNITE_10.0_classifier.qza
db=$(basename ${cls%%_*})

reads=ITS_ITSxpress_rep_seq.se.qza

qiime feature-classifier classify-sklearn \
    --verbose \
    --p-reads-per-batch 10000 \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --i-reads 2.asv/${reads} \
    --i-classifier ${cls} \
    --o-classification 4.taxonomy/${reads%%_*}.${db}.qza
```

JOBID: 46363835

# 6. Export results

```bash
mkdir -p 5.export

module purge
module load QIIME2/2023.5

# ASV table
for i in 16S_feature_table ITS_ITSxpress_feature_table.se; do
    unzip -p 2.asv/${i}.qza \
      */data/feature-table.biom \
      > 5.export/${i%%_*}.feature_table.biom
done

# Phylogeny
for i in 3.phylogeny/*.rooted_tree.qza; do
    unzip -p ${i} \
        */data/tree.nwk \
        > 5.export/$(basename ${i} .qza).nwk
done

# Taxonomy
for i in 4.taxonomy/*; do
    unzip -p ${i} \
        */data/taxonomy.tsv \
        > 5.export/$(basename ${i} .qza)_taxonomy.tsv
done
```

# Backup

```bash
# QC information
tar -cvf ../project_space/1.qc.tar 1.qc/

# Read-based outputs
tar -cvf ../project_space/2.asv.tar \
    2.asv/ITS_ITSxpress.manifest \
    2.asv/ITS_ITSxpress_demux.qza \
    2.asv/ITS_ITSxpress_denoise_stats.se.qza \
    2.asv/ITS_ITSxpress_rep_seq.se.qza \
    2.asv/ITS_ITSxpress_feature_table.se.qza \
    2.asv/16S.manifest \
    2.asv/16S_demux.qza \
    2.asv/16S_denoise_stats.qza \
    2.asv/16S_rep_seq.qzv \
    2.asv/16S_feature_table.qza

tar -cvf ../project_space/3.phylogeny.tar 3.phylogeny/*

# Database and classifiers
tar -cvf ../project_space/db.tar db/*

# Taxonomies
tar -cvf ../project_space/4.taxonomy.tar 4.taxonomy/*

# Exported tables
tar -cvzf ../project_space/5.export.tar.gz 5.export/*
```
