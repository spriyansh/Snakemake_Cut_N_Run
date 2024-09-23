# Use config values
SEACR_PATH = config["tool_paths"]["SEACR_PATH"]
seacr_peakType = config["params"]["peakType"]
seacr_significance = config["params"]["significance"]

all_files_map = {
    item["key"]: item["value"] for item in config["exp_control"]
}

hg38_normalized_bedgraph_dir = config["file_paths"]["hg38_normal_bedgraph_dir"]
hg38_SEACR_PEAKS_dir = config["file_paths"]["hg38_seacr_peaks_dir"]


rule run_annotate_peaks:
    input:
        seacr_results_folder=expand(
            "{path}/{sample}.{peaktype}.bed",
            path=config["file_paths"]["hg38_seacr_peaks_dir"],
            sample=all_files_map.keys(),
            peaktype=["control", "top"],
        ),
    output:
        peak_count_table="{}/tables/06_A_Count_Table.txt".format(
            config["file_paths"]["RData"]
        ),
        peak_metadata_table="{}/tables/06_B_Peak_Metadata_Table.txt".format(
            config["file_paths"]["RData"]
        ),
    shell:
        """
        Rscript Rscripts/06_Annotate_Peaks.R {peak_path} {processed_bam_path} {table_path} {gtf_path}
                                """.format(
            peak_path=config["file_paths"]["hg38_seacr_peaks_dir"],
            processed_bam_path=config["file_paths"]["hg38_processed_dir"],
            table_path=config["file_paths"]["RData"] + "/tables/",
            gtf_path=config["index"]["hg38_gtf"],
        )
