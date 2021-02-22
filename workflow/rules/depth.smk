rule mosdepth:
    input:
        bam = "results/merged/{sample}.bam",
        bai = "results/merged/{sample}.bam.bai",
    output:
        multiext('results/depth/{sample}/{sample}','.mosdepth.global.dist.txt','.mosdepth.region.dist.txt', ".mosdepth.summary.txt", ".regions.bed.gz",".regions.bed.gz.csi")
    threads:
        12
    resources:
        time=20,
        mem=10000,
        cpus=12
    params:
        pfx = 'results/depth/{sample}/{sample}',
        ws = config.get('MOSDEPTH_WINDOW_SIZE'),
    conda:
        '../envs/mosdepth.yaml'
    shell:
        """
        mosdepth -n -t {threads} --by {params.ws} --use-median {params.pfx} {input.bam}
        """

rule copies:
    input:
        cov = 'results/depth/{sample}/{sample}.mosdepth.summary.txt',
        fasta = config.get('CONSENSUS_TE_FASTA')
    output:
        tsv = 'results/copies/{sample}.tsv'
    resources:
        time=30,
        mem=20000,
        cpus=2
    conda:
        "../envs/bioc-general.yaml"
    script:
        "../scripts/copies.R"
