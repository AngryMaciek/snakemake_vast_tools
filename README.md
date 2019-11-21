# Snakemake pipeline for VAST-TOOLS execution
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[General information about the project]

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
Specify all the required information (input/output/parameters) in the config.yaml 

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


---------------------

* design table description, single / paired, forward first?
* hsa mmu only
* description here
* no dots inside fastq sapmle names and sample IDs



Before each run adjust the config.yaml with the right information (input/output/parameters).
The main input to the pipeline is the design table which has to have the following format:

FORMAT
