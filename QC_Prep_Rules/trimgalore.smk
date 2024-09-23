#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: Trim the fastq files using Trimgalore
"""


rule trimgalore:
    priority: 1
    input:
        fwd=f"{config['file_paths']['raw_fq_dir']}/{{sample}}_1.fq.gz",
        rev=f"{config['file_paths']['raw_fq_dir']}/{{sample}}_2.fq.gz",
    output:
        fwd_paired=f"{config['file_paths']['trimmed_fq_dir']}/{{sample}}_1.fq.gz",
        rev_paired=f"{config['file_paths']['trimmed_fq_dir']}/{{sample}}_2.fq.gz",
        trimming_report_fwd=f"{config['file_paths']['trimmed_fq_dir']}/{{sample}}_1.trim_report.txt",
        trimming_report_rev=f"{config['file_paths']['trimmed_fq_dir']}/{{sample}}_2.trim_report.txt",
    threads: config["params"]["trimgalore_threads"]
    shell:
        f"""
        mkdir -p logs
        {config['tool_paths']['TRIM_GALORE_PATH']} --cores {{threads}} --paired {{input.fwd}} {{input.rev}} -o {config['file_paths']['trimmed_fq_dir']}
        mv {config['file_paths']['trimmed_fq_dir']}/{{wildcards.sample}}_1.fq.gz_trimming_report.txt {{output.trimming_report_fwd}}
        mv {config['file_paths']['trimmed_fq_dir']}/{{wildcards.sample}}_2.fq.gz_trimming_report.txt {{output.trimming_report_rev}}
        mv {config['file_paths']['trimmed_fq_dir']}/{{wildcards.sample}}_1_val_1.fq.gz {{output.fwd_paired}}
        mv {config['file_paths']['trimmed_fq_dir']}/{{wildcards.sample}}_2_val_2.fq.gz {{output.rev_paired}}
        """
