#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Deeptools
"""

# Set variables
deeptools_PATH = config["tool_paths"]["deeptools_PATH"]
SAMPLES = config["samples"]
SAMPLES_ANTIBODY = [
    sample
    for sample in SAMPLES
    if not any(ctrl in sample for ctrl in config["control"])
]
hg38_Regions = config["index"]["hg38_gtf"]


rule compute_matrix:
    input:
        # List all .bw files here; this is an example
        bw_files=expand(
            "{path}/{sample}.bw",
            path=config["file_paths"]["hg38_deeptools_heatmaps_dir"],
            sample=SAMPLES,
        ),
        hg38_Regions=hg38_Regions,
    output:
        matrix_3K_5K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_3K_5K.mat.gz",
        matrix_2K_3K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_2K_3K.mat.gz",
        matrix_1K_2K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_1K_2K.mat.gz",
        matrix_5H_1K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_5H_1K.mat.gz",
        matrix_25_5H=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/matrix_gene_25_5H.mat.gz",
    threads: 24  # Adjust as needed
    shell:
        """
        # 5000 Body Length and 3000 Length on either side of TSS
        computeMatrix scale-regions -S {input.bw_files} \
                                      -R {input.hg38_Regions} \
                                      --beforeRegionStartLength 3000 \
                                      --regionBodyLength 5000 \
                                      --afterRegionStartLength 3000 \
                                      --skipZeros -o {output.matrix_3K_5K} -p {threads}

        # 3000 Body Length and 2000 Length on either side of TSS
        computeMatrix scale-regions -S {input.bw_files} \
                                      -R {input.hg38_Regions} \
                                      --beforeRegionStartLength 2000 \
                                      --regionBodyLength 3000 \
                                      --afterRegionStartLength 2000 \
                                      --skipZeros -o {output.matrix_2K_3K} -p {threads}

        # 2000 Body Length and 1000 Length on either side of TSS
        computeMatrix scale-regions -S {input.bw_files} \
                                      -R {input.hg38_Regions} \
                                      --beforeRegionStartLength 1000 \
                                      --regionBodyLength 2000 \
                                      --afterRegionStartLength 1000 \
                                      --skipZeros -o {output.matrix_1K_2K} -p {threads}
        # 1000 Body Length and 500 Length on either side of TSS
        computeMatrix scale-regions -S {input.bw_files} \
                                      -R {input.hg38_Regions} \
                                      --beforeRegionStartLength 500 \
                                      --regionBodyLength 1000 \
                                      --afterRegionStartLength 500 \
                                      --skipZeros -o {output.matrix_5H_1K} -p {threads}
        # 500 Body Length and 250 Length on either side of TSS
        computeMatrix scale-regions -S {input.bw_files} \
                                      -R {input.hg38_Regions} \
                                      --beforeRegionStartLength 250 \
                                      --regionBodyLength 500 \
                                      --afterRegionStartLength 250 \
                                      --skipZeros -o {output.matrix_25_5H} -p {threads}
        """


rule create_heatmaps:
    input:
        matrix_3K_5K=rules.compute_matrix.output.matrix_3K_5K,
        matrix_2K_3K=rules.compute_matrix.output.matrix_2K_3K,
        matrix_1K_2K=rules.compute_matrix.output.matrix_1K_2K,
        matrix_5H_1K=rules.compute_matrix.output.matrix_5H_1K,
        matrix_25_5H=rules.compute_matrix.output.matrix_25_5H,
    output:
        png_3K_5K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B5k_L3K.png",
        png_2K_3K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B3k_L2K.png",
        png_1K_2K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B2k_L1K.png",
        png_5H_1K=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B1k_L5H.png",
        png_25_5H=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/Histone_Signals_heatmap_B5H_L250.png",
    threads: 24  # Adjust as needed
    shell:
        """
        plotHeatmap -m {input.matrix_3K_5K} -out {output.png_3K_5K} --sortUsing sum
        plotHeatmap -m {input.matrix_2K_3K} -out {output.png_2K_3K} --sortUsing sum
        plotHeatmap -m {input.matrix_1K_2K} -out {output.png_1K_2K} --sortUsing sum
        plotHeatmap -m {input.matrix_5H_1K} -out {output.png_5H_1K} --sortUsing sum
        plotHeatmap -m {input.matrix_25_5H} -out {output.png_25_5H} --sortUsing sum
        """

rule extract_regions:
    input:
        controlPeakBed=f"{config['file_paths']['hg38_seacr_peaks_dir']}/{{SAMPLES_ANTIBODY}}.control.bed",
    output:
        summitRegionBed=f"{config['file_paths']['hg38_deeptools_heatmaps_centers_dir']}/{{SAMPLES_ANTIBODY}}_summitRegion.bed",
    shell:
        """
        awk '{{split($6, summit, ":"); split(summit[2], region, "-"); printf "%s\\t%d\\t%d\\n", summit[1], region[1], region[2]}}' {input.controlPeakBed} > {output.summitRegionBed}
        """

rule compute_matrix_summit:
    input:
        raw_bw=f"{config['file_paths']['hg38_deeptools_heatmaps_dir']}/{{SAMPLES_ANTIBODY}}.bw",
        hg38_Regions=rules.extract_regions.output.summitRegionBed,
    output:
        matrix=f"{config['file_paths']['hg38_deeptools_heatmaps_centers_dir']}/{{SAMPLES_ANTIBODY}}.centered.mat.gz",
    threads: 8,
    shell:
        """
        computeMatrix reference-point -S {input.raw_bw} \
                                      -R {input.hg38_Regions} \
                                      --skipZeros -p {threads} \
                                      --referencePoint center \
                                      -a 3000 -b 3000 \
                                      -o {output.matrix}
        """

rule plot_compute_matrix_summit:
    input:
        matrix=rules.compute_matrix_summit.output.matrix,
    output:
        png=f"{config['file_paths']['hg38_deeptools_heatmaps_centers_dir']}/{{SAMPLES_ANTIBODY}}.png"
    threads: 8
    shell:
        """
        plotHeatmap -m {input.matrix} -out {output.png} \
                    --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" \
                    --xAxisLabel "" --regionsLabel "Peaks"
        """