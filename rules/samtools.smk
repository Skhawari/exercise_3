rule sam_to_bam:
    input:
        "results/sam/{sample}.sam"
    output:
        "results/bam/{sample}.bam"
    log:
        "logs/samtools/sam_to_bam/{sample}.log"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -Sb {input} > {output} 2> {log}
        """

rule sort_bam:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/bam_sorted/{sample}.sorted.bam"
    log:
        "logs/samtools/sort_bam/{sample}.log"
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2> {log}
        """

rule index_bam:
    input:
        "results/bam_sorted/{sample}.sorted.bam"
    output:
        "results/bam_sorted/{sample}.sorted.bam.bai"
    log:
        "logs/samtools/index_bam/{sample}.log"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input} &> {log}
        """

rule idxstats:
    input:
        bam = "results/bam_sorted/{sample}.sorted.bam",
        bai = "results/bam_sorted/{sample}.sorted.bam.bai"
    output:
        "results/stats/{sample}.idxstats.txt"
    log:
        "logs/samtools/idxstats/{sample}.log"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output} 2> {log}
        """

rule aggregate_idxstats:
    input:
        expand("results/stats/{sample}.idxstats.txt", sample=SAMPLES)
    output:
        "results/stats/aggregated_idxstats.tsv"
    log:
        "logs/samtools/aggregate_idxstats.log"
    run:
        import pandas as pd
        import os

        aggregated = None
        for infile in input:
            sample = os.path.basename(infile).split(".")[0]
            df = pd.read_csv(infile, sep="\t", header=None,
                             names=["seqname", "seqlen", "mapped", "unmapped"])
            df = df[["seqname", "seqlen", "mapped"]].rename(columns={"mapped": sample})
            if aggregated is None:
                aggregated = df
            else:
                aggregated = pd.merge(aggregated, df[["seqname", sample]], on="seqname")

        # Reorder columns (seqname, seqlen, then samples)
        cols = ["seqname", "seqlen"] + [col for col in aggregated.columns if col not in ["seqname", "seqlen"]]
        aggregated = aggregated[cols]

        aggregated.to_csv(output[0], sep="\t", index=False)