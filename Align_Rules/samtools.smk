#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Samtools based rules
"""


# Set variables
samtools_PATH = config["tool_paths"]["samtools_PATH"]
SAMPLES = config["samples"]
SAMPLES_Antibody = [sample for sample in config["samples"] if "IgG" not in sample]


def get_samples_without_IgG(wildcards):
    return [
        f"{config['file_paths']['hg38_processed_dir']}/{sample}.bam"
        for sample in config["samples"]
        if "IgG" not in sample
    ]


rule calc_spikeIn_SeqDep:
    priority: 3
    input:
        spikeInSAM=rules.bowtie2_spikeIn.output.sam,  #f"{config['file_paths']['spikeIn_bowtie2_sam_dir']}/{{sample}}.sam",
    output:
        spikeInSeqDepth=f"{config['file_paths']['spikeIn_bowtie2_sam_dir']}/{{sample}}.seqDepth",
    threads: config["params"]["samtools_threads"]
    shell:
        """
        seqDepthRaw=$({samtools_PATH} view -@ {threads} -F 0x04 {input.spikeInSAM} | wc -l)
        seqDepth=$(($seqDepthRaw / 2))
        echo $seqDepth > {output.spikeInSeqDepth}
        """


rule samtools_conversion:
    input:
        sam=rules.bowtie2.output.sam,  #f"{config['file_paths']['hg38_bowtie2_sam_dir']}/{{sample}}.sam",
    output:
        bam=f"{config['file_paths']['hg38_bowtie2_bam_dir']}/{{sample}}.bam",
        bai=f"{config['file_paths']['hg38_bowtie2_bam_dir']}/{{sample}}.bam.bai",
    threads: config["params"]["samtools_threads"]
    shell:
        f"""
        {samtools_PATH} view -Sb {{input.sam}} -o {{output.bam}}.temp
        {samtools_PATH} sort -@ {{threads}} {{output.bam}}.temp -o {{output.bam}}
        rm {{output.bam}}.temp
        {samtools_PATH} index {{output.bam}} {{output.bai}}
        """


rule fragment_size_estimation:
    input:
        bam=rules.samtools_conversion.output.bam,
    output:
        fragmentLen=f"{config['file_paths']['hg38_fragmentLen_dir']}/{{sample}}.txt",
    threads: config["params"]["samtools_threads"]
    shell:
        r"""
        {samtools_PATH} view -F 0x04 {input.bam} \
        | awk  -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' \
        | sort | uniq -c | awk -v OFS="\t" '{{print $2, $1/2}}' > {output.fragmentLen}
        """


rule summarize_fragment_size:
    input:
        hg38_logs=expand(
            "{path}/{sample}.txt",
            path=config["file_paths"]["hg38_fragmentLen_dir"],
            sample=SAMPLES,
        ),
    output:
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
    shell:
        """
                                                                                Rscript Rscripts/02_fragment_size_comparison.R {hg38_path} {img_path} {table_path}
                                                                                """.format(
            hg38_path=config["file_paths"]["hg38_fragmentLen_dir"],
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )


rule sort_post_process_bam:
    input:
        bam=rules.samtools_conversion.output.bam,
    output:
        processedBam=f"{config['file_paths']['hg38_processed_dir']}/{{sample}}.bam",
    params:
        minQualityScore=config["params"]["minQualityScoreAlignment"],
    threads: config["params"]["samtools_threads"]
    shell:
        f"""
        {samtools_PATH} sort -@ {{threads}} -n {{input.bam}} -O BAM -o temp_sorted_{{wildcards.sample}}.bam
        {samtools_PATH} view -bS -F 0x04 temp_sorted_{{wildcards.sample}}.bam > {{output.processedBam}}
        rm temp_sorted_{{wildcards.sample}}.bam
        """


rule deeptools_heatmaps_prepare:
    input:
        processedBam=lambda wildcards: f"{config['file_paths']['hg38_processed_dir']}/{wildcards.sample}.bam",
    output:
        raw_bw=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/{{sample}}.bw",
        run_log=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/{{sample}}.run.log",
    params:
        minQualityScore=config["params"]["minQualityScoreAlignment"],
    threads: config["params"]["samtools_threads"]
    shell:
        f"""
        {samtools_PATH} sort -@ {{threads}} {{input.processedBam}} -O BAM -o temp_sorted_{{wildcards.sample}}.bam
        {samtools_PATH} index temp_sorted_{{wildcards.sample}}.bam temp_sorted_{{wildcards.sample}}.bam.bai
        bamCoverage -b temp_sorted_{{wildcards.sample}}.bam -o {{output.raw_bw}} &> {{output.run_log}}
        """
