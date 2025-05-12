# rules/qc.smk

rule fastqc_raw:
    input:
        lambda wc: [] if config["skip_fastqc"] else [
            samples.at[wc.sample, "fq1"],
            samples.at[wc.sample, "fq2"]
        ]
    output:
        html1 = "results/fastqc_raw/{sample}_1_fastqc.html",
        html2 = "results/fastqc_raw/{sample}_2_fastqc.html",
        zip1 = "results/fastqc_raw/{sample}_1_fastqc.zip",
        zip2 = "results/fastqc_raw/{sample}_2_fastqc.zip"
    log:
        "logs/qc/fastqc_raw/{sample}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc {input[0]} {input[1]} \
            --outdir results/fastqc_raw &> {log}
        """

rule trimmomatic:
    input:
        fq1 = lambda wc: samples.at[wc.sample, "fq1"],
        fq2 = lambda wc: samples.at[wc.sample, "fq2"],
        adapter = lambda wc: config["adapter_file"] if not config["skip_trimming"] else []
    output:
        fq1 = "results/trimmed/{sample}_1.trimmed.fastq.gz",
        fq2 = "results/trimmed/{sample}_2.trimmed.fastq.gz"
    log:
        "logs/qc/trimmomatic/{sample}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.fq1} {input.fq2} \
            {output.fq1} /dev/null \
            {output.fq2} /dev/null \
            ILLUMINACLIP:{input.adapter}:{config[trimmomatic_params]} &> {log}
        """

rule fastqc_trimmed:
    input:
        lambda wc: [] if config["skip_fastqc"] else [
            f"results/trimmed/{wc.sample}_1.trimmed.fastq.gz",
            f"results/trimmed/{wc.sample}_2.trimmed.fastq.gz"
        ]
    output:
        html1 = "results/fastqc_trimmed/{sample}_1.trimmed_fastqc.html",
        html2 = "results/fastqc_trimmed/{sample}_2.trimmed_fastqc.html",
        zip1 = "results/fastqc_trimmed/{sample}_1.trimmed_fastqc.zip",
        zip2 = "results/fastqc_trimmed/{sample}_2.trimmed_fastqc.zip"

    log:
        "logs/qc/fastqc_trimmed/{sample}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc {input[0]} {input[1]} \
            --outdir results/fastqc_trimmed &> {log}
        """

rule qualimap:
    input:
        lambda wc: [] if config["skip_qualimap"] else [f"results/bam_sorted/{wc.sample}.sorted.bam"]
    output:
        html = "results/qualimap/{sample}/qualimapReport.html"
    log:
        "logs/qc/qualimap/{sample}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        qualimap bamqc \
            -bam {input[0]} \
            -outdir results/qualimap/{wildcards.sample} &> {log}
        """

rule multiqc:
    input:
        expand("results/fastqc_trimmed/{sample}_1.trimmed_fastqc.zip", sample=SAMPLES)
    output:
        html = "results/multiqc/multiqc_report.html"
    log:
        "logs/qc/multiqc.log"
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc results/ -o results/multiqc --filename multiqc_report.html &> {log}
        """