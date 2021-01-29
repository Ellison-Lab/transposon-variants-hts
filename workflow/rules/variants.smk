def get_fastqs(wc):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]
    tmp2 = tmp[tmp['subsample_name'] == wc.subsample]
    return [tmp2.get('fastq_r1')[0], tmp2.get('fastq_r2')[0]]

rule bwa_mem2_index:
    """
    Index consensus transposon sequences for use with bwa-mem2
    """
    input:
        config.get('CONSENSUS_TE_FASTA')
    output:
        multiext("results/idx/transposons",".0123",".amb",".ann",".bwt.2bit.64",".bwt.8bit.32",".pac")
    log:
        "results/logs/bwa-mem2_index/transposons.log"
    params:
        prefix="results/idx/transposons"
    wrapper:
        "0.70.0/bio/bwa-mem2/index"

rule bwa_mem2_mem:
    input:
        reads = get_fastqs,
        idx = rules.bwa_mem2_index.output,
    output:
        temp("results/mapped/{sample}-{subsample}.bam")
    log:
        "results/logs/bwa_mem2/{sample}-{subsample}.log"
    params:
        index="results/idx/transposons",
        #extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        extra=lambda wc: r"-R '@RG\tID:{r}\tSM:{s}\tLB:{l}'".format(s=wc.sample, r=pep.get_sample(wc.sample).rgid, l=wc.subsample),
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.70.0/bio/bwa-mem2/mem"

rule samtools_fixmate:
    input:
        "results/mapped/{sample}-{subsample}.bam"
    output:
        temp("results/fixed/{sample}-{subsample}.bam")
    threads:
        4
    params:
        extra = ""
    wrapper:
        "0.70.0/bio/samtools/fixmate/"

rule samtools_sort:
    input:
        "results/fixed/{sample}-{subsample}.bam"
    output:
        temp("results/sorted/{sample}-{subsample}.bam")
    params:
        extra = "-m 4G",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.70.0/bio/samtools/sort"

rule picard_mark_duplicates:
    input:
        "results/sorted/{sample}-{subsample}.bam"
    output:
        bam=temp("results/marked/{sample}-{subsample}.bam"),
        metrics="results/marked/{sample}-{subsample}.metrics.txt"
    log:
        "results/logs/picard_mark_duplicates/{sample}-{subsample}.log"
    params:
        "REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT",
    wrapper:
        "0.70.0/bio/picard/markduplicates"

rule samtools_merge:
    input:
        lambda wc: expand("results/marked/{s}-{sub}.bam",s=wc.sample,sub=pep.get_sample(wc.sample).subsample_name)
    output:
        "results/merged/{sample}.bam"
    params:
        "" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.70.0/bio/samtools/merge"

rule samtools_index:
    input:
        "results/merged/{sample}.bam"
    output:
        "results/merged/{sample}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.70.0/bio/samtools/index"

rule freebayes:
    input:
        ref=config.get('CONSENSUS_TE_FASTA'),
        # you can have a list of samples here
        samples=expand("results/merged/{s}.bam",s=SAMPLES),
        # the matching BAI indexes have to present for freebayes
        indexes=expand("results/merged/{s}.bam.bai",s=SAMPLES),
    output:
        "results/freebayes/calls.bcf"  # either .vcf or .bcf
    log:
        "results/logs/freebayes.log"
    params:
        extra="--pooled-continuous -F 0.001 -C 3 -p 1 --min-base-quality 20",         # freq-based calling
        chunksize=100000,
        normalize=False,  # flag to use bcftools norm to normalize indels
    threads: 8
    wrapper:
        "0.70.0/bio/freebayes"

rule bcftools_filter:
    """
    filter and make calls not-multiallelic
    """
    input:
        rules.freebayes.output
    output:
        vcf="results/freebayes/calls.filt.vcf.gz",
        tbi="results/freebayes/calls.filt.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools norm -m-any {input} | \
        bcftools filter -Oz -o {output.vcf} \
            -i '%QUAL>10 && SAP>10 && SRP>10 && EPP>10' {input} &&
        tabix {output.vcf}
        """
