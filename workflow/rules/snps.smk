rule get_pileups:
    input:
        bam = "results/merged/{sample}.bam",
        bai = "results/merged/{sample}.bam.bai",
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
        vcfs = expand("results/snps/{s}-snps.vcf",s=SAMPLES),
    params:
        min_snp_depth =  config.get('MIN_SNP_SUPPORT')
    conda:
        "../envs/bioc-general.yaml"
    script:
        '../scripts/snps.R'


# rule bcftools_pileups:
#     input:
#         bams = expand("results/merged/{s}.bam",s=SAMPLES),
#         bai = expand("results/merged/{s}.bam.bai",s=SAMPLES),
#         fa = config.get("COMBINED_GENOME"),
#         te_fa = config.get("CONSENSUS_TE_FASTA")
#     output:
#         vcf = 'results/bcftools/pileups/{te}.vcf.gz',
#     threads:
#         2
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#         bcftools mpileup -p -I -a FORMAT/AD --threads {threads} -Oz -o {output.vcf} -d 1000000 -r {wildcards.te} -f {input.fa} {input.bams}
#         """
#
# rule bcftools_call:
#     input:
#         vcf = rules.bcftools_pileups.output.vcf,
#     output:
#         vcf = 'results/bcftools/call/{te}.vcf.gz',
#     threads:
#         2
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#          bcftools call --ploidy 1 -A -v -m --threads {threads} -Oz -o {output.vcf} {input.vcf}
#         """
#
# rule bcftools_norm:
#     input:
#         vcf = rules.bcftools_call.output.vcf,
#     output:
#         vcf = 'results/bcftools/norm/{te}.vcf.gz',
#     threads:
#         2
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#         bcftools norm -m-any -Oz --threads {threads} -o {output.vcf} {input.vcf}
#         """
#
# rule bcftools_filt:
#     input:
#         vcf = rules.bcftools_norm.output.vcf,
#     output:
#         vcf = 'results/bcftools/filt/{te}.vcf.gz',
#     threads:
#         2
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#         bcftools filter --threads {threads} -Oz -o {output.vcf} -i 'FORMAT/AD==0' {input.vcf}
#         """


## bcftools merged
## bcftools isec
        # """
        # REGIONS=$(cat {input.te_fa} | grep '>' | cut -f 1 -d '#' | tr -d '>' | tr '\n' ',')
        # let 'THREADS = {threads} - 3'
        # bcftools mpileup -I --threads $THREADS -Ou -d 1000000 -r $REGIONS -f {input.fa} {input.bams} | \
        #     bcftools call -v -m -Ou |
        #     bcftools norm -m-any -Ou | \
        #     bcftools +fill-tags -t AO -Ou | \
        #     bcftools filter -Oz -o {output.vcf} -i 'FORMAT/AO==0' && \
        # tabix {output.vcf}
        #
        # """
