##############################################################################
#
#   Snakemake pipeline for VAST-TOOLS
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 20-11-2019
#   LICENSE: Apache_2.0
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

# get samples marked as untreated
def get_CONTROL_samples():
    design_table = pd.read_csv(config["design_file"], sep="\t")
    mask = design_table["condition"] == "untreated"
    return list(design_table[mask]["sample"])

# get samples marked as treated
def get_TREATMENT_samples():
    design_table = pd.read_csv(config["design_file"], sep="\t")
    mask = design_table["condition"] == "treated"
    return list(design_table[mask]["sample"])

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
        TSV_compare_outfile = expand( \
            os.path.join("{output_dir}", "VT_compare_output.tsv"), \
            output_dir=config["output_dir"]),
        TSV_diff_outfile = expand( \
            os.path.join("{output_dir}", "INCLUSION_LEVELS_FULL.DIFF.txt"), \
            output_dir=config["output_dir"])

##############################################################################
### Before even the analysis starts:
### In case VASTDB is not present - download it
##############################################################################

rule VASTDB_download:
    output:
        DIR_VASTDB = directory(config["VASTDB"])
    shell:
        """
        (mkdir {output.DIR_VASTDB} && \
        cd {output.DIR_VASTDB} && \
        wget http://vastdb.crg.eu/libs/vastdb.hsa.16.02.18.tar.gz && \
        tar xzvf vastdb.hsa.16.02.18.tar.gz && \
        wget http://vastdb.crg.eu/libs/vastdb.mmu.16.02.18.tar.gz && \
        tar xzvf vastdb.mmu.16.02.18.tar.gz) \
        2>&1
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
        LIST_samples = get_fastq,
        DIR_db_dir = config["VASTDB"]
    output:
        DIR_align_dir = directory("{output_dir}/vast_tools_align_{sample}")
    params:
        species = config["species"],
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
        --dbDir {input.DIR_db_dir} \
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
        DIR_grouped_dir = os.path.join("{output_dir}", "grouped"),
        DIR_db_dir = config["VASTDB"]
    output:
        TEMP_combine = temp(os.path.join("{output_dir}", "VT_combine.temp"))
    params:
        species = config["species"],
        annotation_update = annotation_update(config["species"]),
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
        --dbDir {input.DIR_db_dir}; \
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
        TSV_design_table = config["design_file"],
        DIR_grouped_dir = directory(os.path.join("{output_dir}", "grouped")),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "adjust_inclusion_table.log"),
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "adjust_inclusion_table.log")
    run:
        # VAST_TOOLS takes one of the FASTQ files name as sample name
        # to be the header of a column in the inclusion table.
        # Adjust that for the real sample sames from the design table.
        fname_regex = os.path.join(params.DIR_grouped_dir, "INCLUSION_*")
        df = pd.read_csv(glob.glob(fname_regex)[0], sep="\t", header=0)
        design_table = \
            pd.read_csv(params.TSV_design_table, sep="\t", index_col=0)
        name_dict = {}
        for sample,row in design_table.iterrows():
            old_name = row.fq1.split("/")[-1].split(".")[0]
            name_dict[old_name] = sample
        new_headers = list(df.columns.values[:6])
        for col in df.columns.values[6:]:
            if col in name_dict.keys():
                new_headers.append(name_dict[col])
            elif col in [x+"-Q" for x in name_dict.keys()]:
                new_headers.append(name_dict[col[:-2]]+"-Q")
            else:
                assert False # should not happen
        df.columns = new_headers
        df.to_csv(output.TSV_inclusion_table, sep="\t", index=False)

#################################################################################
### Compare treated/untreated samples to each other
#################################################################################

rule VT_compare:
    input:
        TSV_inclusion_table = \
            os.path.join("{output_dir}", "INCLUSION_LEVELS_FULL.tsv"),
        DIR_db_dir = config["VASTDB"]
    output:
        TSV_compare_outfile = \
            os.path.join("{output_dir}", "VT_compare_output.tsv")
    params:
        TSV_compare_outfile_name = "VT_compare_output.tsv",
        STRING_c_samples = ",".join(get_CONTROL_samples()),
        STRING_t_samples = ",".join(get_TREATMENT_samples()),
        species = config["species"],
        min_dpsi = config["min_dPSI"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "VT_compare.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "VT_compare.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", "cluster_log", \
            "VT_compare_benchmark.log")
    singularity:
        "docker://vastgroup/vast-tools:v2.2.2"
    shell:
        """
        vast-tools compare \
        {input.TSV_inclusion_table} \
        --min_range -100 \
        --min_dPSI {params.min_dpsi} \
        --dbDir {input.DIR_db_dir} \
        --sp {params.species} \
        --GO \
        -a {params.STRING_c_samples} \
        -b {params.STRING_t_samples} \
        --outFile {params.TSV_compare_outfile_name} \
        &> {log.LOG_local_log}
        """

#################################################################################
### Differential Splicing Analysis
#################################################################################

rule VT_diff:
    input:
        TSV_inclusion_table = \
            os.path.join("{output_dir}", "INCLUSION_LEVELS_FULL.tsv")
    output:
        TSV_diff_outfile = \
            os.path.join("{output_dir}", "INCLUSION_LEVELS_FULL.DIFF.txt")
    params:
        MV_param = config["MV_param"],
        min_coverage = config["min_coverage"],
        DIR_outdir = "{output_dir}",
        STRING_c_samples = ",".join(get_CONTROL_samples()),
        STRING_t_samples = ",".join(get_TREATMENT_samples()),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "VT_diff.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "VT_diff.log")
    resources:
        threads = 8,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", "cluster_log", \
            "VT_diff_benchmark.log")
    singularity:
        "docker://vastgroup/vast-tools:v2.2.2"
    shell:
        """
        vast-tools diff \
        -o {params.DIR_outdir} \
        -a {params.STRING_c_samples} \
        -b {params.STRING_t_samples} \
        -e {params.min_coverage} \
        -c {resources.threads} \
        -r {params.MV_param} \
        2> {log.LOG_local_log}
        """
