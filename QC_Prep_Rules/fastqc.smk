#!/usr/bin/env python3

"""
Name: Priyansh Srivastava
Email: spriyansh29@gmail.com
Web: https://www.metapriyansh.com/
Script Title: RUN FASTQC on the Raw FastQ Files and Trimmed FastQ Files
"""


rule fastqc_Raw:
    priority: 1
    input:
        fwd=f"{config['file_paths']['raw_fq_dir']}/{{sample}}_1.fq.gz",
        rev=f"{config['file_paths']['raw_fq_dir']}/{{sample}}_2.fq.gz",
    output:
        fwd_html=f"{config['file_paths']['qc_before_trim_dir']}/{{sample}}_1_fastqc.html",
        rev_html=f"{config['file_paths']['qc_before_trim_dir']}/{{sample}}_2_fastqc.html",
        fwd_zip=f"{config['file_paths']['qc_before_trim_dir']}/{{sample}}_1_fastqc.zip",
        rev_zip=f"{config['file_paths']['qc_before_trim_dir']}/{{sample}}_2_fastqc.zip",
    threads: config["params"]["fastqc_threads"]
    log:
        out=f"{config['file_paths']['qc_before_trim_dir']}/logs/{{sample}}_fastqc.out",
        err=f"{config['file_paths']['qc_before_trim_dir']}/logs/{{sample}}_fastqc.err",
    shell:
        f"{config['tool_paths']['FASTQC_PATH']} -q -t {{threads}} -o {config['file_paths']['qc_before_trim_dir']} {{input.fwd}} {{input.rev}}"


rule fastqc_Trim:
    priority: 2
    input:
        fwd=rules.trimgalore.output.fwd_paired,
        rev=rules.trimgalore.output.rev_paired,
    output:
        fwd_html=f"{config['file_paths']['qc_after_trim_dir']}/{{sample}}_1_fastqc.html",
        rev_html=f"{config['file_paths']['qc_after_trim_dir']}/{{sample}}_2_fastqc.html",
        fwd_zip=f"{config['file_paths']['qc_after_trim_dir']}/{{sample}}_1_fastqc.zip",
        rev_zip=f"{config['file_paths']['qc_after_trim_dir']}/{{sample}}_2_fastqc.zip",
    threads: config["params"]["fastqc_threads"]
    log:
        out=f"{config['file_paths']['qc_after_trim_dir']}/logs/{{sample}}_fastqc.out",
        err=f"{config['file_paths']['qc_after_trim_dir']}/logs/{{sample}}_fastqc.err",
    shell:
        f"{config['tool_paths']['FASTQC_PATH']} -q -t {{threads}} -o {config['file_paths']['qc_after_trim_dir']} {{input.fwd}} {{input.rev}}"
