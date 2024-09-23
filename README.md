# SKM1 CUT & RUN with Snakemake
Snakemake Pipeline to analyze data from SKM-1 leukemic cell line produced using Cleavage Under Targets and Release Using Nuclease protocol sequenced on Illumina NovaSeq 6000. The Pipeline is based on the [CUT&Tag Data Processing and Analysis Tutorial](https://yezhengstat.github.io/CUTTag_tutorial/index.html).

## Docker Image
The docker image to reproduce the results of the analysis is available at
[Docker-Hub](https://hub.docker.com/layers/spriyansh29/cnr_snakemake/dev/images/sha256-09097d493ed10382d234c300740f6cac35d185652887347df241aafbac576bbf?context=explore)

### Docker run
To run the docker, docker should be installed on the system. Information about how to install Docker is available at the [docker website](https://docs.docker.com/engine/install/). After the installation of Docker, execute the following command,
```
docker run -it --rm --name priyansh_cnr_snakemake -v /data/priyansh_data/skm1_cut_n_run/snakemake_pipeline_merged/:/outputFolder -v /data/priyansh_data/skm1_cut_n_run/absolute_raw/lane_merged/:/rawFQ -v /home/priyansh/index_files/:/indexFolder -v /data/priyansh_data/skm1_cut_n_run/SKM1_CUT_N_RUN_Snakemake/:/scripts spriyansh29/cnr_snakemake:snakemake
```

## Run the pipeline
To run the pipeline execute the following command after going to the directory scripts

### Move to scripts
```
cd /scripts
```
### Execute the command
```
snakemake -s main.smk -c 24 --configfile config_merged.json --rerun-incomplete -k --dry-run
```


## Results Diretcory Structure 
---