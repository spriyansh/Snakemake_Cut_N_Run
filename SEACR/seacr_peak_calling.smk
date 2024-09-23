# Use config values
SEACR_PATH = config["tool_paths"]["SEACR_PATH"]
seacr_peakType = config["params"]["peakType"]
seacr_significance = config["params"]["significance"]

all_files_map = {
    item["key"]: item["value"] for item in config["exp_control"]
}

hg38_normalized_bedgraph_dir = config["file_paths"]["hg38_normal_bedgraph_dir"]
hg38_SEACR_PEAKS_dir = config["file_paths"]["hg38_seacr_peaks_dir"]


rule run_SEACR_against_control:
    input:
        target=f"{hg38_normalized_bedgraph_dir}/{{sample}}.bedgraph",
        control=lambda wildcards: f"{hg38_normalized_bedgraph_dir}/{all_files_map[wildcards.sample]}.bedgraph",
    output:
        against_control=f"{hg38_SEACR_PEAKS_dir}/{{sample}}.control.bed",  # Removed .bed here
    threads: 24
    shell:
        f"""
        bash {SEACR_PATH} {{input.target}} {{input.control}} norm stringent {{wildcards.sample}}
        mv {{wildcards.sample}}.stringent.bed {{output.against_control}}
        """


rule run_SEACR_top_peaks:
    input:
        target=f"{hg38_normalized_bedgraph_dir}/{{sample}}.bedgraph",
    output:
        top_peaks=f"{hg38_SEACR_PEAKS_dir}/{{sample}}.top.bed",  # Removed .bed here
    threads: 24
    shell:
        f"""
        bash {SEACR_PATH} {{input.target}} {seacr_significance} non stringent {{wildcards.sample}}
        mv {{wildcards.sample}}.stringent.bed {{output.top_peaks}}
        """


rule summarize_peaks:
    input:
        seacr_results_folder=expand(
            "{path}/{sample}.{peaktype}.bed",
            path=config["file_paths"]["hg38_seacr_peaks_dir"],
            sample=all_files_map.keys(),
            peaktype=["control", "top"],
        ),
    output:
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
    shell:
        """
                                    Rscript Rscripts/05_summarize_peaks.R {peak_path} {img_path} {table_path}
                                """.format(
            peak_path=config["file_paths"]["hg38_seacr_peaks_dir"],
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )
