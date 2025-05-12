# rules/bowtie.smk

rule build_bowtie2_index:
    input:
        fasta = config["reference"]
    output:
        expand("results/bowtie2_index/index.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ])
    log:
        "logs/bowtie2/index_build.log"
    threads: 1
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        bowtie2-build {input.fasta} results/bowtie2_index/index \
        &> {log}
        """

rule bowtie2_map:
    input:
        r1 = lambda wc: samples.at[wc.sample, "fq1"],
        r2 = lambda wc: samples.at[wc.sample, "fq2"],
        index = expand("results/bowtie2_index/index.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ])
    output:
        "results/sam/{sample}.sam"
    log:
        "logs/bowtie2/{sample}.log"
    threads: 4
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        bowtie2 \
            --{config[bowtie2][sensitivity]} \
            -N {config[bowtie2][seed_mismatch]} \
            -x results/bowtie2_index/index \
            -1 {input.r1} -2 {input.r2} \
            -S {output} -p {threads} \
            &> {log}
        """