
rule bwa_mem2_index:
    """
    Index genome
    """
    input:
        custom_genome('results/custom-genome/combined.fasta')
    output:
        multiext("results/idx/idx",".0123",".amb",".ann",".bwt.2bit.64",".bwt.8bit.32",".pac")
    log:
        "results/logs/bwa-mem2_index/transposons.log"
    resources:
        time=60,
        mem=20000,
        cpus=4
    params:
        prefix="results/idx/idx"
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/bwa-mem2/index"

rule bwa_mem2_mem:
    input:
        reads = rules.trim_qual.output,
        idx = rules.bwa_mem2_index.output,
    output:
        temp("results/mapped/{sample}-{subsample}.bam")
    log:
        "results/logs/bwa_mem2/{sample}-{subsample}.log"
    resources:
        time=120,
        mem=32000,
        cpus=24
    params:
        index="results/idx/idx",
        #extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        extra=lambda wc: r"-R '@RG\tID:{r}\tSM:{s}\tLB:{l}'".format(s=wc.sample, r=pep.get_sample(wc.sample).rgid, l=wc.subsample),
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 24
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/bwa-mem2/mem"

rule samtools_fixmate:
    input:
        "results/mapped/{sample}-{subsample}.bam"
    output:
        temp("results/fixed/{sample}-{subsample}.bam")
    threads:
        4
    resources:
        cpus=4
    params:
        extra = ""
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.70.0/bio/samtools/fixmate/"

rule samtools_sort:
    input:
        "results/fixed/{sample}-{subsample}.bam"
    output:
        temp("results/sorted/{sample}-{subsample}.bam")
    resources:
        cpus=8
    params:
        extra = "-m 4G",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.70.0/bio/samtools/sort"

rule picard_mark_duplicates:
    input:
        "results/sorted/{sample}-{subsample}.bam"
    output:
        bam=temp("results/marked/{sample}-{subsample}.bam"),
        metrics="results/marked/{sample}-{subsample}.metrics.txt"
    log:
        "results/logs/picard_mark_duplicates/{sample}-{subsample}.log"
    params:
        "REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT",
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.70.0/bio/picard/markduplicates"

rule filter_bwa_reads:
    input:
        "results/marked/{sample}-{subsample}.bam"
    output:
        "results/filt/{sample}-{subsample}.bam"
    params:
        "-O BAM -F 256 -q {mapq}".format(mapq = config.get('MIN_BM2_MAPQ'))
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/samtools/view"

rule samtools_merge:
    input:
        lambda wc: expand("results/filt/{s}-{sub}.bam",s=wc.sample,sub=pep.get_sample(wc.sample).subsample_name)
    output:
        "results/merged/{sample}.bam"
    resources:
        cpus=8
    params:
        "" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.70.0/bio/samtools/merge"
