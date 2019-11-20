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
import shutil
import glob
import pandas as pd

# local rules
localrules: all, VASTDB_download, create_output_dir, \
    group_alignments, adjust_inclusion_table

# get FASTQ files related to all samples
def get_fastq(w):
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    x = [design_table.at[w.sample,"fq1"], design_table.at[w.sample,"fq2"]]
    x = [i for i in x if str(i)!="nan"]
    return x

# get all samples from the design table
def get_samples():
    design_table = pd.read_csv(config["design_file"], sep="\t")
    return list(design_table["sample"])

#def get_CONTROL_samples():
#    design_table = pd.read_csv(config["design_file"],sep="\t")
#    return list(design_table[design_table["condition"]=="untreated"]["sample"])

#def get_TREATMENT_samples():
#    design_table = pd.read_csv(config["design_file"],sep="\t")
#    return list(design_table[design_table["condition"]=="treated"]["sample"])

# add flags to force use new annotation versions
def annotation_update(species):
    if species=="Hsa":
        return "-a hg38"
    elif species=="Mmu":
        return "-a mm10"
    else:
        return ""

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TXT_final_results = expand(directory("{output_dir}/INCLUSION_LEVELS_FULL.tsv"),output_dir=config["output_dir"])

##############################################################################
### Before even the analysis starts:
### In case VASTDB is not present - download it
##############################################################################

rule VASTDB_download:
    output:
        DIR_VASTDB = config["VASTDB"]
    shell:
        """
        vast-tools align \
        {input.LIST_samples} \
        --dbDir {params.db_dir} \
        --output {output.DIR_align_dir} \
        --sp {params.species} \
        --useFastq \
        --EEJ_counts \
        --expr \
        -cores {resources.threads} \
        &> {log.LOG_local_log}
        """

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
        "docker://vastgroup/vast-tools:v2.2.2"
    shell:
        """
        vast-tools align \
        {input.LIST_samples} \
        --dbDir {params.db_dir} \
        --output {output.DIR_align_dir} \
        --sp {params.species} \
        --useFastq \
        --EEJ_counts \
        --expr \
        -cores {resources.threads} \
        &> {log.LOG_local_log}
        """

##############################################################################
### Restructure the alignment results for further
### VAST-TOOLS analysis
##############################################################################

rule group_alignments:
    input:
        DIR_align_dir = \
            expand(directory("{output_dir}/vast_tools_align_{sample}"), \
                output_dir=config["output_dir"], sample=get_samples())
    output:
        DIR_grouped_dir = directory(os.path.join("{output_dir}", "grouped"))
    params:
        DIR_outdir = "{output_dir}",
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "group_alignments.log"),
        queue = "30min",
        time = "00:30:00"
    resources:
        threads = 1,
        mem = 5000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "group_alignments.log")
    run:
        # copy aligned dirs into new location
        os.makedirs(output.DIR_grouped_dir)
        os.makedirs(os.path.join(output.DIR_grouped_dir, "to_combine"))
        dirs = [os.path.join(i, "to_combine") for i in input.DIR_align_dir]
        for path in dirs:
            src_files = os.listdir(path)
            for file_name in src_files:
                shutil.copy(os.path.join(path, file_name), \
                    os.path.join(output.DIR_grouped_dir, "to_combine"))

#################################################################################
### Combine aligned results
#################################################################################

rule VT_combine:
    input:
        DIR_grouped_dir = os.path.join("{output_dir}", "grouped")
    output:
        TEMP_combine = temp(os.path.join("{output_dir}", "VT_combine.temp"))
    params:
        species = config["species"],
        annotation_update = annotation_update(config["species"]),
        db_dir = "VASTDB",
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "VT_combine.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "VT_combine.log"),
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", "cluster_log", \
            "VT_combine_benchmark.log")
    singularity:
        "docker://vastgroup/vast-tools:v2.2.2"
    shell:
        """
        vast-tools combine \
        --o {input.DIR_grouped_dir} \
        --sp {params.species} \
        {params.annotation_update} \
        --dbDir {params.db_dir}; \
        touch {output.TEMP_combine} \
        &> {log.LOG_local_log}
        """

#################################################################################
### Normalize the final inclusion table
#################################################################################

rule adjust_inclusion_table:
    input:
        TEMP_combine = os.path.join("{output_dir}", "VT_combine.temp")
    output:
        TSV_inclusion_table = \
            os.path.join("{output_dir}", "INCLUSION_LEVELS_FULL.tsv")
    params:
        DIR_grouped_dir = directory(os.path.join("{output_dir}", "grouped")),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "adjust_inclusion_table.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "adjust_inclusion_table.log")
    resources:
        threads = 1,
        mem = 5000
    run:
        fname_regex = os.path.join(params.DIR_grouped_dir, "INCLUSION_*")
        df = pd.read_csv(glob.glob(fname_regex)[0], sep="\t", header=0)
        #remove_suffix = "_" in df.columns.values[7]
        #if remove_suffix:
        #    df.columns = list(df.columns.values[:6]) + [c.split("_")[0]+c.split("_")[1][1:] for c in df.columns.values[6:]]
        df.to_csv(output.TSV_inclusion_table, sep="\t", index=False)














#################################################################################
### Compare the samples to each other
#################################################################################

#rule VT_compare:
#    input:
#        inclusion_table = "{output_dir}/INCLUSION_LEVELS_FULL.tab"
#    output:
#        outfile = "{output_dir}/compare_output.tsv"
#    params:
#        c_samples = ",".join(get_CONTROL_samples()),
#        t_samples = ",".join(get_TREATMENT_samples()),
#        species = config["species"],
#        min_dpsi = config["min_dPSI"],
#        dbDir = "VASTDB",
#        cluster_log = "{output_dir}/cluster_log/vast_tools_compare.log",
#        queue = "6hours",
#        time = "6:00:00"
#    log:
#        local_log = "{output_dir}/local_log/vast_tools_compare.log",
#    resources:
#        threads = 1,
#        mem = 5000
#    benchmark:
#        "{output_dir}/cluster_log/vast_tools_compare.benchmark.log"
#    singularity:
#        "docker://zavolab/vast-tools:2.0.2"
#    shell:
#        """
#        vast-tools compare \
#        {input.inclusion_table} \
#        --min_range -100 \
#        --min_dPSI {params.min_dpsi} \
#        --dbDir {params.dbDir} \
#        --sp {params.species} \
#        --GO \
#        -a {params.c_samples} \
#        -b {params.t_samples} \
#        --outFile compare_output.tsv \
#        &> {log.local_log}
#        """











#################################################################################
### Differential Splicing Analysis
#################################################################################

#rule VT_diff:
#    input:
#        inclusion_table = "{output_dir}/INCLUSION_LEVELS_FULL.tab"
#    output:
#        outfile = "{output_dir}/diff_output.tsv"
#    params:
#        MV_param = config["MV_param"],
#        outdir = "{output_dir}",
#        c_samples = ",".join(get_CONTROL_samples()),
#        t_samples = ",".join(get_TREATMENT_samples()),
#        min_coverage = config["min_coverage"],
#        cluster_log = "{output_dir}/cluster_log/vast_tools_diff.log",
#        queue = "6hours",
#        time = "6:00:00"
#    log:
#        local_log = "{output_dir}/local_log/vast_tools_diff.log",
#    resources:
#        threads = 8,
#        mem = 5000
#    benchmark:
#        "{output_dir}/cluster_log/vast_tools_diff.benchmark.log"
#    singularity:
#        "docker://zavolab/vast-tools:2.0.2"
#    shell:
#        """
#        vast-tools diff \
#        -o {params.outdir} \
#        -a {params.c_samples} \
#        -b {params.t_samples} \
#        -e {params.min_coverage} \
#        -c {resources.threads} \
#        -r {params.MV_param} \
#        1> {output.outfile} \
#        2> {log.local_log}
#        """
