def get_fastqs(wc, r="r1"):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]
    tmp2 = tmp[tmp['subsample_name'] == wc.subsample]
    return [tmp2.get('fastq_'+r)[0]]

rule get_fqs:
    input:
        r1 = lambda wc: get_fastqs(wc, "r1"),
        r2 = lambda wc: get_fastqs(wc, "r2"),
    output:
        r1 = "results/fastq/{sample}/{subsample}_r1.fastq.gz"
        r2 = "results/fastq/{sample}/{subsample}_r2.fastq.gz"
    shell:
        """
        wget -O {output.r1} {input.r1} &&
        wget -O {output.r2} {input.r2}
        """

rule trim_qual:
    input:
        r1 = rules.get_fqs.output.r1
        r2 = rules.get_fqs.output.r2
    output:
        temp("results/fastq-trim-qual/{sample}/{subsample}_r1.fastq"),
        temp("results/fastq-trim-qual/{sample}/{subsample}_r2.fastq"),
    params:
        q = config.get('MIN_FQ_QUAL', 20),
        min_read_len = config.get('MIN_READ_LEN', 35)
    resources:
        time=60,
        mem=5000,
        cpus=8
    conda:
        "../envs/cutadapt.yaml"
    threads:
        8
    shell:
        "cutadapt -q {params.q} -m {params.min_read_len} -j {threads} -o {output[0]} -p {output[1]} {input.r1} {input.r2}"
