##############################################################################
#
#   Snakemake pipeline for VAST-TOOLS
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 20-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################

# imports
import sys
import os

# local rules
localrules: create_output_dir, all

# get FASTQ files related to all samples
def get_fastq(wildcards):
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    x = [design_table.at[wildcards.sample,"fq1"], design_table.at[wildcards.sample,"fq2"]]
    x = [i for i in x if str(i)!="nan"]
    return x

# get all samples from the design table
def get_samples():
    design_table = pd.read_csv(config["design_file"], sep="\t")
    return list(design_table["sample"])

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TXT_final_results = expand(directory("{output_dir}/vast_tools_align_{sample}"),output_dir=config["output_dir"], sample=get_samples())

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_random_samples}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Align sequencing reads
##############################################################################

rule VT_align:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        LIST_samples = get_fastq
    output:
        DIR_align_dir = directory("{output_dir}/vast_tools_align_{sample}")
    params:
        species = config["species"],
        db_dir = "VASTDB",
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "VT_align_{sample}.log"),
        queue = "1day",
        time = "23:59:59"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "VT_align_{sample}.log"),
    resources:
        threads = 4,
        mem = 100000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "VT_align_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    singularity:
        "docker://zavolab/vast-tools:2.0.2"
    shell:
        """
        vast-tools align \
        {input.samples} \
        --dbDir {params.dbDir} \
        --output {output.align_dir} \
        --sp {params.species} \
        --useFastq \
        --EEJ_counts \
        --expr \
        -cores {resources.threads} \
        &> {log.local_log}
        """






##############################################################################
### Merge the results
##############################################################################

#rule merge_results:
#    input:
#        TXT_result_files = \
#            lambda wildcards: [os.path.join(wildcards.output_dir,
#                "random_samples", f) for f in config["samples_filenames"]]
#    output:
#        TXT_final_results = os.path.join("{output_dir}", "results.txt")
#    params:
#        LOG_cluster_log = \
#            os.path.join("{output_dir}", "cluster_log/merge_results.log"),
#        queue = "30min",
#        time = "00:05:00"
#    resources:
#        threads = 1,
#        mem = 5000
#    log:
#        LOG_local_log = \
#            os.path.join("{output_dir}", "local_log", "merge_results.log")
#    run:
#        # read all the sampled numbers:
#        numbers = []
#        for i in input.TXT_result_files:
#            with open(i) as f:
#                numbers.append(f.read().splitlines()[0])
#        # save into one file:
#        with open(output.TXT_final_results, "w") as outfile:
#                outfile.write("\n".join(numbers))

