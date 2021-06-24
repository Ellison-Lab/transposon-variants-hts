rule trim_qual:
    input:
        r1 = lambda wc: AUTO.remote(pep.get_sample(wc.sample).fastq_r1),
        r2 = lambda wc: AUTO.remote(pep.get_sample(wc.sample).fastq_r2),
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
