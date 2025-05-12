# Snakefile
import pandas as pd
import os

configfile: "config/config.yaml"

# Lade Sample-Tabelle und extrahiere Namen
samples = pd.read_csv(config["samples"], sep="\t", index_col="sample")
SAMPLES = list(samples.index)

include: "rules/qc.smk"
include: "rules/bowtie.smk"
include: "rules/samtools.smk"

#rule all:
#    input:
#        expand("results/stats/{sample}.idxstats.txt", sample=SAMPLES),
 #       "results/stats/aggregated_idxstats.tsv",
  #      expand("results/fastqc_raw/{sample}_1_fastqc.html", sample=SAMPLES),
   #     expand("results/fastqc_trimmed/{sample}_1.trimmed_fastqc.html", sample=SAMPLES),
    #    expand("results/qualimap/{sample}/qualimapReport.html", sample=SAMPLES),
     #   "results/multiqc/multiqc_report.html"


rule all:
    input:
        expand("results/stats/{sample}.idxstats.txt", sample=SAMPLES),
        "results/stats/aggregated_idxstats.tsv",
        expand("results/fastqc_raw/{sample}_1_fastqc.html", sample=SAMPLES) if not config["skip_fastqc"] else [],
        expand("results/fastqc_trimmed/{sample}_1.trimmed_fastqc.html", sample=SAMPLES) if not config["skip_fastqc"] and not config["skip_trimming"] else [],
        expand("results/qualimap/{sample}/qualimapReport.html", sample=SAMPLES) if not config["skip_qualimap"] else [],
        "results/multiqc/multiqc_report.html" if not all([config["skip_fastqc"], config["skip_qualimap"]]) else []