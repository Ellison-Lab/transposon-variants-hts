def get_fastqs(wc):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]
    tmp2 = tmp[tmp['subsample_name'] == wc.subsample]
    return [tmp2.get('fastq_r1')[0], tmp2.get('fastq_r2')[0]]

rule trim_qual:
    input:
        reads = get_fastqs
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
        "cutadapt -q {params.q} -m {params.min_read_len} -j {threads} -o {output[0]} -p {output[1]} {input.reads}"
