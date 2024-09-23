#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Global Alignment with bowtie2
"""

# Diabled Dovetail "--dovetail"

hg38_Index = ind = config["index"]["bowtie2_hg_38_index"]
spikeIn_Index = ind = config["index"]["bowtie2_spikeIn_index"]
SAMPLES = config["samples"]


rule bowtie2:
    priority: 2
    input:
        fwd=rules.trimgalore.output.fwd_paired,
        rev=rules.trimgalore.output.rev_paired,
    output:
        sam=f"{config['file_paths']['hg38_bowtie2_sam_dir']}/{{sample}}.sam",
    log:
        log_file=f"{config['file_paths']['hg38_bowtie2_sam_dir']}/{{sample}}.log",
    threads: config["params"]["bowtie2_threads"]
    shell:
        f"""
        {config['tool_paths']['bowtie2_PATH']} --end-to-end --threads {{threads}} --very-sensitive \
        --no-mixed --no-discordant --phred33 -I 25 -X 700 \
        -x {hg38_Index} \
        -1 {{input.fwd}} -2 {{input.rev}} -S {{output.sam}} &> {{log.log_file}}
        """


rule bowtie2_spikeIn:
    priority: 2
    input:
        fwd=rules.trimgalore.output.fwd_paired,
        rev=rules.trimgalore.output.rev_paired,
    output:
        sam=f"{config['file_paths']['spikeIn_bowtie2_sam_dir']}/{{sample}}.sam",
    log:
        log_file=f"{config['file_paths']['spikeIn_bowtie2_sam_dir']}/{{sample}}.log",
    threads: config["params"]["bowtie2_threads"]
    shell:
        f"""
        {config['tool_paths']['bowtie2_PATH']} --end-to-end --threads {{threads}} --very-sensitive \
        --no-mixed --no-discordant --phred33 -I 25 -X 700 \
        -x {spikeIn_Index} \
        -1 {{input.fwd}} -2 {{input.rev}} -S {{output.sam}} &> {{log.log_file}}
        """


# Summarize
rule summarize_alignments:
    input:
        hg38_logs=expand(
            "{path}/{sample}.log",
            path=config["file_paths"]["hg38_bowtie2_sam_dir"],
            sample=SAMPLES,
        ),
        spikeIn_logs=expand(
            "{path}/{sample}.log",
            path=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            sample=SAMPLES,
        ),
    output:
        sequencing_depth_png="{}/images/01_A_Sequencing_Depth.png".format(
            config["file_paths"]["RData"]
        ),
        alignable_fragment_png="{}/images/01_B_Alignable_Fragment_hg38.png".format(
            config["file_paths"]["RData"]
        ),
        alignable_rate_png="{}/images/01_C_Alignment_Rate_hg38.png".format(
            config["file_paths"]["RData"]
        ),
        alignable_rate_spikeIn_png="{}/images/01_D_Alignment_Rate_SpikeIn.png".format(
            config["file_paths"]["RData"]
        ),
        combined_alignment_summary_png="{}/images/01_Combined_Alignment_Summary.png".format(
            config["file_paths"]["RData"]
        ),
        alignment_summary_table="{}/tables/01_Alignment_Summary.txt".format(
            config["file_paths"]["RData"]
        ),
    shell:
        """
                                                                                                                Rscript Rscripts/01_summarize_alignments.R \
                                                                                                                    {hg38_path} \
                                                                                                                    {spikeIn_path} \
                                                                                                                    {img_path} \
                                                                                                                    {table_path}
                                                                                                                """.format(
            hg38_path=config["file_paths"]["hg38_bowtie2_sam_dir"],
            spikeIn_path=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )
