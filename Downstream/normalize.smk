#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Normalzie
"""

# Set variables
bedtools_PATH = config["tool_paths"]["bedtools_PATH"]
hg38_chromSize = config["index"]["hg38_chromSize"]
hg38_genome_bed = config["index"]["hg38_genome_bed"]


# bin fragments
rule calc_norm_factors:
    priority: 8
    input:
        hg38_logs=expand(
            "{path}/{sample}.seqDepth",
            path=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            sample=SAMPLES,
        ),
        align_summary=rules.summarize_alignments.output.alignment_summary_table,
    output:
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
    shell:
        """
                                                                                            Rscript Rscripts/04_calculate_scaling_factor.R {hg38_path} {img_path} {table_path}
                                                                                                                """.format(
            hg38_path=config["file_paths"]["spikeIn_bowtie2_sam_dir"],
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )


# Normalize
rule normalize_hg38:
    input:
        hg38CleanFragmentBed=rules.fragment_convert.output.fragmentBed,
        scale_factor_file=rules.calc_norm_factors.output.scale_table,
    output:
        normalizedBedgraph=f"{config['file_paths']['hg38_normal_bedgraph_dir']}/{{sample}}.bedgraph",
    run:
        # Read the tab-separated file into a dictionary
        scale_factors = {}
        # Genralize this step
        with open(input.scale_factor_file) as file:
            next(file)  # skip header line
            for line in file:
                parts = line.strip().split("\t")
                scale_factors[parts[0]] = float(parts[1])
        sample_name = os.path.basename(input.hg38CleanFragmentBed)
        scale_factor = scale_factors.get(sample_name, None)
        if scale_factor is None:
            raise ValueError(f"No scale factor found for sample: {sample_name}")
        shell(
            f"""
            bedtools genomecov -bg -scale {scale_factor} -i {input.hg38CleanFragmentBed} -g {hg38_chromSize} > {output.normalizedBedgraph}
            """
        )


# Normlaized Bed
rule convert_to_normalized_bed:
    priority: 9
    input:
        hg38_Normalized_Bedgraph=rules.normalize_hg38.output.normalizedBedgraph,
    output:
        hg38_Normalized_Bed=f"{config['file_paths']['hg38_normal_bed_dir']}/{{sample}}.bed",
    shell:
        f"""
        {bedtools_PATH} intersect -a {{input.hg38_Normalized_Bedgraph}} -b {hg38_genome_bed} > {{output.hg38_Normalized_Bed}}
        """
