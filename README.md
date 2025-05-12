# Snakemake Workflow – Exercise 2 (Applied Sequence Analysis)

This repository contains a modular and configurable Snakemake workflow for quality control and mapping steps in NGS preprocessing. It was developed collaboratively as part of Exercise 3 in the "Applied Sequence Analysis" course at Freie Universität Berlin (SS2025).

## 💡 What this workflow does

For each paired-end FASTQ sample, it performs:

1. **Indexing** of the reference genome with `bowtie2-build`
2. **Mapping** with `bowtie2`
3. **Conversion** from SAM to BAM (`samtools view`)
4. **Sorting** of BAM files (`samtools sort`)
5. **Indexing** of sorted BAM files (`samtools index`)
6. **Read count statistics** with `samtools idxstats`
7. **Aggregation** of per-sample read counts into a summary table

## 📁 Repository structure

```
.
├── config/             # Contains config.yaml with paths and parameters
│   └── samples.csv     # Sample information (paths to FASTQ files)
├── resources/          # Contains reference genome (reference.fa)
├── results/            # Output data: sam, bam, sorted bam, idxstats
├── workflow/
│   ├── envs/           # Conda environment YAMLs (bowtie2.yaml, samtools.yaml)
│   ├── rules/          # Modular Snakemake rule files
│   └── Snakefile       # Main Snakemake workflow
├── .gitignore
├── LICENSE.md
└── README.md
```

## ⚙️ How to run

First, activate conda and run:

```bash
conda config --set channel_priority strict
snakemake --use-conda --cores 4
```

## ✅ Output

Final summary file (after complete run):

```
results/summary/idxstats_summary.tsv
```

## 🔬 Developed by

- Dr. Sandro Andreotti – Lecturer, Freie Universität Berlin
- Michael Riethmüller – Master's Student in Bioinformatics, Freie Universität Berlin
- Sajjad Khawari – Master's Student in Bioinformatics, Freie Universität Berlin