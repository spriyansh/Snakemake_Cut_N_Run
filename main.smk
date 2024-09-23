#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Main file for the snakemake
Description: Snakemake Pipeline to analyze
    data from SKM-1 leukemic cell line generated
    using CUT&RUN (Cleavage Under Targets and Release Using
    Nuclease) protocol, sequenced on Illumina NovaSeq 6000.
"""


# Set configuration File
# configfile: "config.json"

# # Copy and process samples
# dup_samples = []
# copy_samples = []

# # Add Logic
# if config["params"]["remove_duplicates"] == "control":
#     for sample in config["samples"]:
#         if config["params"]["control_antibody"] in sample:
#             dup_samples.append(sample)
#         else:
#             copy_samples.append(sample)
# elif config["params"]["remove_duplicates"] == "all":
#     dup_samples = config["samples"]
# elif config["params"]["remove_duplicates"] == "none":
#     copy_samples = config["samples"]

# def get_bam_path(sample):
#     if sample in dup_samples:
#         return f"{config['file_paths']['hg38_picard_dir']}/{sample}.rem.bam"
#     else:
#         return f"{config['file_paths']['hg38_picard_dir']}/{sample}.mark.bam"


# Include the rules
include: "QC_Prep_Rules/trimgalore.smk"
include: "QC_Prep_Rules/fastqc.smk"
include: "Align_Rules/bowtie2.smk"
include: "Align_Rules/samtools.smk"
include: "Downstream/bedtools.smk"
include: "Downstream/normalize.smk"
include: "Downstream/deeptools.smk"
include: "SEACR/seacr_peak_calling.smk"
include: "SEACR/chip_seekr_peak_annotation.smk"
include: "Differential_Exp/deseq2_diff_exp.smk"


SAMPLES = config["samples"]
SAMPLES_ANTIBODY = [
    sample
    for sample in SAMPLES
    if not any(ctrl in sample for ctrl in config["control"])
]
# Rule to specify which outputs you want to achieve
rule all:
    input:
        # FastQC on all samples before Trimming
        expand(
            "{fastqc_dir}/{sample}_1_fastqc.html",
            fastqc_dir=config["file_paths"]["qc_before_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_2_fastqc.html",
            fastqc_dir=config["file_paths"]["qc_before_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_1_fastqc.zip",
            fastqc_dir=config["file_paths"]["qc_before_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_2_fastqc.zip",
            fastqc_dir=config["file_paths"]["qc_before_trim_dir"],
            sample=config["samples"],
        ),
        # Trimming reads with trimGalore
        expand(
            "{trimmed_read_dir}/{sample}_1.fq.gz",
            trimmed_read_dir=config["file_paths"]["trimmed_fq_dir"],
            sample=config["samples"],
        ),
        expand(
            "{trimmed_read_dir}/{sample}_2.fq.gz",
            trimmed_read_dir=config["file_paths"]["trimmed_fq_dir"],
            sample=config["samples"],
        ),
        expand(
            "{trimmed_read_dir}/{sample}_1.trim_report.txt",
            trimmed_read_dir=config["file_paths"]["trimmed_fq_dir"],
            sample=config["samples"],
        ),
        expand(
            "{trimmed_read_dir}/{sample}_2.trim_report.txt",
            trimmed_read_dir=config["file_paths"]["trimmed_fq_dir"],
            sample=config["samples"],
        ),
        # FastQC on all samples after Trimming
        expand(
            "{fastqc_dir}/{sample}_1_fastqc.html",
            fastqc_dir=config["file_paths"]["qc_after_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_2_fastqc.html",
            fastqc_dir=config["file_paths"]["qc_after_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_1_fastqc.zip",
            fastqc_dir=config["file_paths"]["qc_after_trim_dir"],
            sample=config["samples"],
        ),
        expand(
            "{fastqc_dir}/{sample}_2_fastqc.zip",
            fastqc_dir=config["file_paths"]["qc_after_trim_dir"],
            sample=config["samples"],
        ),
        # SAM Alignments
        expand(
            "{sam_dir}/{sample}.sam",
            sam_dir=config["file_paths"]["hg38_bowtie2_sam_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.sam",
            sam_dir=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            sample=config["samples"],
        ),
        ## Summarize Alignment
        expand(
            "{sam_dir}/{sample}.seqDepth",
            sam_dir=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bam",
            sam_dir=config["file_paths"]["hg38_bowtie2_bam_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bam.bai",
            sam_dir=config["file_paths"]["hg38_bowtie2_bam_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.txt",
            sam_dir=config["file_paths"]["hg38_fragmentLen_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bam",
            sam_dir=config["file_paths"]["hg38_processed_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bed",
            sam_dir=config["file_paths"]["hg38_bed_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bed",
            sam_dir=config["file_paths"]["hg38_clean_bed_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bed",
            sam_dir=config["file_paths"]["hg38_fragment_bed_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bed",
            sam_dir=config["file_paths"]["hg38_bin_500_bed_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bedgraph",
            sam_dir=config["file_paths"]["hg38_normal_bedgraph_dir"],
            sample=config["samples"],
        ),
        expand(
            "{sam_dir}/{sample}.bed",
            sam_dir=config["file_paths"]["hg38_normal_bed_dir"],
            sample=config["samples"],
        ),
        # expand(
        #     "{bw_dir}/{sample}.bw",
        #     bw_dir=config["file_paths"]["hg38_deeptools_heatmaps_dir"],
        #     sample=config["samples"],
        # ),
        expand(
            f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/{{sample}}.bw",
            sample=SAMPLES_Antibody,
        ),
        # Compute Matrix and Make Heatmaps
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_3K_5K.mat.gz",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_2K_3K.mat.gz",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_1K_2K.mat.gz",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_5H_1K.mat.gz",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_25_5H.mat.gz",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B5k_L3K.png",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B3k_L2K.png",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B2k_L1K.png",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B1k_L5H.png",
        f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B5H_L250.png",
        ## SEACR
        expand(
            "{dir}/{sample}.control.bed",
            dir=hg38_SEACR_PEAKS_dir,
            sample=all_files_map.keys(),
        ),
        expand(
            "{dir}/{sample}.top.bed",
            dir=hg38_SEACR_PEAKS_dir,
            sample=all_files_map.keys(),
        ),
        expand(
            "{dir}/{sample}_summitRegion.bed",
            dir=config['file_paths']['hg38_deeptools_heatmaps_centers_dir'],
            sample=SAMPLES_Antibody,
        ),
        expand(
            "{dir}/{sample}.centered.mat.gz",
            dir=config['file_paths']['hg38_deeptools_heatmaps_centers_dir'],
            sample=SAMPLES_Antibody,
        ),
        expand(
            "{dir}/{sample}.png",
            dir=config['file_paths']['hg38_deeptools_heatmaps_centers_dir'],
            sample=SAMPLES_Antibody,
        ),
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
        fragment_length_violin_png="{}/images/02_A_Fragment_Length_Violin.png".format(
            config["file_paths"]["RData"]
        ),
        fragment_length_png="{}/images/02_B_Fragment_Length.png".format(
            config["file_paths"]["RData"]
        ),
        combined_fragment_length_png="{}/images/02_Combined_Fragment_Length.png".format(
            config["file_paths"]["RData"]
        ),
        fragment_length_summary_table="{}/tables/02_Fragment_Length.txt".format(
            config["file_paths"]["RData"]
        ),
        corrplot_png="{}/images/03_A_Corrplot.png".format(config["file_paths"]["RData"]),
        corrplot_table="{}/tables/03_Correlation_Matrix.txt".format(
            config["file_paths"]["RData"]
        ),
        spikeIn_scaling_factor_png="{}/images/04_A_SpikeIn_Scaling_Factor.png".format(
            config["file_paths"]["RData"]
        ),
        normalization_fragment_count_png="{}/images/04_B_Normalization_Fragment_Count.png".format(
            config["file_paths"]["RData"]
        ),
        scale_normal_png="{}/images/04_Combined_Scaling_Normalization.png".format(
            config["file_paths"]["RData"]
        ),
        new_alignSummary_table="{}/tables/04_Alignment_Summary_Scaling.txt".format(
            config["file_paths"]["RData"]
        ),
        scale_table="{}/tables/04_Scaling_factors.txt".format(
            config["file_paths"]["RData"]
        ),
        peak_reproduction_png="{}/images/05_A_Number_of_Peaks.png".format(
            config["file_paths"]["RData"]
        ),
        peak_numbers_png="{}/images/05_B_Peak_Reproduction.png".format(
            config["file_paths"]["RData"]
        ),
        combined_peak_summary_png="{}/images/05_Combined_Peak_Calling.png".format(
            config["file_paths"]["RData"]
        ),
        seace_peak_summary_table="{}/tables/05_SEACR_peak_summary.txt".format(
            config["file_paths"]["RData"]
        ),
        peak_count_table="{}/tables/06_A_Count_Table.txt".format(
            config["file_paths"]["RData"]
        ),
        peak_metadata_table="{}/tables/06_B_Peak_Metadata_Table.txt".format(
            config["file_paths"]["RData"]
        ),
        sample_metadata_table="{}/tables/07_Sample_Metadata.txt".format(
            config["file_paths"]["RData"]
        ),
