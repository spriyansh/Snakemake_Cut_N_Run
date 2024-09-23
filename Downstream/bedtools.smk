#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Bed tool convert and Normalization
"""

# Set variables
bedtools_PATH = config["tool_paths"]["bedtools_PATH"]


# Conversion to bed
rule bed_convert:
    input:
        bam=rules.sort_post_process_bam.output.processedBam,
    output:
        bed=f"{config['file_paths']['hg38_bed_dir']}/{{sample}}.bed",
    shell:
        f"""
        {bedtools_PATH} bamtobed -i {{input.bam}} -bedpe > {{output.bed}}
        """


# Bed to Cleanbed
rule clean_bed_convert:
    input:
        bed=rules.bed_convert.output.bed,
    output:
        cleanBed=f"{config['file_paths']['hg38_clean_bed_dir']}/{{sample}}.bed",
    shell:
        """
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input.bed} > {output.cleanBed}
        """


# Clean bed to Fragment
rule fragment_convert:
    input:
        cleanBed=rules.clean_bed_convert.output.cleanBed,
    output:
        fragmentBed=f"{config['file_paths']['hg38_fragment_bed_dir']}/{{sample}}.bed",
    shell:
        f"""
        cut -f 1,2,6 {{input.cleanBed}} | sort -k1,1 -k2,2n -k3,3n > {{output.fragmentBed}}
        """


# bin fragments
rule fragment_binSplit_500:
    input:
        fragmentBed=rules.fragment_convert.output.fragmentBed,
    output:
        fragmentBed500=f"{config['file_paths']['hg38_bin_500_bed_dir']}/{{sample}}.bed",
    shell:
        """
        awk -v w=500 '{{print $1, int(($2 + $3)/(2*w))*w + w/2}}' {input.fragmentBed} | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{{print $2, $3, $1}}' |  sort -k1,1V -k2,2n  > {output.fragmentBed500}
        """


# bin fragments
rule replicate_reproducibility:
    input:
        hg38_logs=expand(
            "{path}/{sample}.bed",
            path=config["file_paths"]["hg38_bin_500_bed_dir"],
            sample=SAMPLES,
        ),
    output:
        corrplot_png="{}/images/03_A_Corrplot.png".format(config["file_paths"]["RData"]),
        corrplot_table="{}/tables/03_Correlation_Matrix.txt".format(
            config["file_paths"]["RData"]
        ),
    shell:
        """
                                                                                                                Rscript Rscripts/03_replicate_reproducibility.R {hg38_path} {img_path} {table_path}
                                                                                                                """.format(
            hg38_path=config["file_paths"]["hg38_bin_500_bed_dir"],
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )
