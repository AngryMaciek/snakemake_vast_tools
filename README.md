# Snakemake pipeline for VAST-TOOLS execution
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[General information about the project]

description here: hsa mmu only

## Snakemake pipeline execution
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires Python 3 and can be most easily installed via the bioconda package from the anaconda cloud service.

### Step 1: Download and installation of Miniconda3
Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  source .bashrc
  ```

macOS:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  source .bashrc
  ```

### Step 2: Pandas and Snakemake installation

To execute the workflow one would require pandas python library and snakemake workflow menager.  
Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest versions:
  ```bash
  conda install -c conda-forge pandas
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify required information (input/output/parameters) in the config.yaml and information about RNA-Seq samples in the design table. The design table is a TSV file with four obligatory columns:
* The first column "sample" would serve as sample unique IDs. Please do not use "." character within the ID.
* Two columns "fq1" and "fq2" contain paths to the RNA-Seq samples. Please do not use "." character within the filename before the file extension. In case of paired-end sequencing data "fq1" should contain forwards read and "fq2" reverse reads. In case of single-end sequencing data please leave "fq2" with empty strings.
* Column "condition" contains either: "treated" or "untreated" mask and refers to the experiment group of a particular sample.

Write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are two scripts to start the pipeline, depending on whether you want to run locally/on a SLURM computational cluster. Subsequent snakemake rules are executed within a singularity container provided by VAST-TOOLS developers (therefore in order to use this pipeline one needs to have `singularity` installed as well). For the cluster execution it might be required to adapt the 'cluster_config.json' before starting the run.
  ```bash
  bash snakemake_local_run_singularity_container.sh
  bash snakemake_cluster_run_singularity_container.sh
  ```

## License

GPL v3.0
