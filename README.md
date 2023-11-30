# ngs-snakemake-workflow
Workflow for next generation sequencing using snakemake for small dataset.

SAMPLES = ['A','B','C']
TRIMMED_INPUTS = expand("data/samples/{sample}_trimmed.fastq", sample=SAMPLES)
FASTQC_INPUTS = expand("data/samples/{sample}.fastq", sample=SAMPLES)
rule all:
   input:
          "qc/multiqc_report.html",
          TRIMMED_INPUTS,
          expand("mapped_reads/{sample}.bam", sample=SAMPLES),
          expand("sorted_reads/{sample}.bam", sample=SAMPLES),
          expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES),
          "calls/all.vcf"
rule fastqc:
    input:
        FASTQC_INPUTS
    output:
        "qc/fastqc/{sample}_fastqc.html",
        "qc/fastqc/{sample}_fastqc.zip"
    conda:
        "environment.yaml"
    shell:
        "fastqc {input} --outdir qc/fastqc"
rule multiqc:
    input:
        expand("qc/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "qc/multiqc_report.html"
    conda:
        "environment.yaml"
    shell:
        "multiqc qc/fastqc --outdir qc"
rule trim:
    input:
        "data/samples/{sample}.fastq"
    output:
        "data/samples/{sample}_trimmed.fastq"
    conda:
        "environment.yaml"
    shell:
        "trimmomatic SE -threads 4 -phred33 "
        "{input} {output} "
        "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 "
        "SLIDINGWINDOW:4:15 MINLEN:36"
rule bwa_map:
    input:
        'data/genome.fa',
        'data/samples/{sample}.fastq',
    output:
        'mapped_reads/{sample}.bam',
    shell:
        'bwa mem {input} | samtools view -Sb - > {output}'

rule samtools_sort:
    input:
        'mapped_reads/{sample}.bam'
    output:
         'sorted_reads/{sample}.bam'
    shell:
         'samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}'
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "environment.yaml"
    shell:
        "samtools index {input}"
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    conda:
        "environment.yaml"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"



