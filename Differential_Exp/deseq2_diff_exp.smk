# Use config values
SEACR_PATH = config["tool_paths"]["SEACR_PATH"]
seacr_peakType = config["params"]["peakType"]
seacr_significance = config["params"]["significance"]

all_files_map = {
    item["key"]: item["value"] for item in config["exp_control"]
}

hg38_normalized_bedgraph_dir = config["file_paths"]["hg38_normal_bedgraph_dir"]
hg38_SEACR_PEAKS_dir = config["file_paths"]["hg38_seacr_peaks_dir"]


rule run_deseq2_combinations:
    input:
        counts = rules.run_annotate_peaks.output.peak_count_table,
        peak_metadata = rules.run_annotate_peaks.output.peak_metadata_table,
    output:
        sample_metadata="{}/tables/07_Sample_Metadata.txt".format(
            config["file_paths"]["RData"]
        ),
    shell:
        """
        Rscript Rscripts/07_deseq2_differential_expression.R {input_dir} {table_path} {img_path} {table_path}
                                """.format(
            input_dir=config["file_paths"]["RData"] + "/tables/",
            img_path=config["file_paths"]["RData"] + "/images/",
            table_path=config["file_paths"]["RData"] + "/tables/",
        )
