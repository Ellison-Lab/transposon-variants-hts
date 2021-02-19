rule get_pileups:
    input:
        bam = "results/merged/{sample}.bam",
        fasta = config.get('CONSENSUS_TE_FASTA')
    output:
        tsv = "results/pileups/{sample}.pileups.tsv.gz"
    params:
        bq = config.get('MIN_SNP_BQ')
    conda:
        "../envs/bioc-general.yaml"
    script:
        '../scripts/pileups.R'


rule get_snps:
    input:
        pileups = expand("results/pileups/{s}.pileups.tsv.gz",s=SAMPLES),
        fasta = config.get('CONSENSUS_TE_FASTA')
    output:
        bed = "results/snps/snps.bed",
        tsv = "results/snps/snps.tsv.gz",
    params:
        min_snp_depth =  config.get('MIN_SNP_SUPPORT')
    conda:
        "../envs/bioc-general.yaml"
    script:
        '../scripts/snps.R'
