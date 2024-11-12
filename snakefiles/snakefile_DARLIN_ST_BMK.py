import os
import sys
from pathlib import Path
import pandas as pd
from darlin import help_functions as hf
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

##################
## preprocessing
################## 
cfg_type=config['cfg_type']
template=config['template']

CARLIN_dir=hf.update_CARLIN_dir(CARLIN_dir,config['template'])
print("Updated CARLIN_dir:"+ str(CARLIN_dir))

if template.startswith('cCARLIN'):
    locus = 'CA'
elif template.startswith('Rosa'):
    locus = 'RA'
elif template.startswith('Tigre'):
    locus = 'TA'

SampleList=config['SampleList']
raw_fastq_dir = config['raw_fastq_dir']

kernel=config['python_DARLIN_pipeline']['kernel']


##################
## start the rules
################## 
rule all:
    input:
        expand("BST_output/{sample}_{locus}/all.done", sample=SampleList, locus=locus)

rule extract_DARLIN_seq:
    input:
        r1=raw_fastq_dir+"/{sample}_{locus}_R1.fastq.gz",
        r2=raw_fastq_dir+"/{sample}_{locus}_R2.fastq.gz"
    output:
        r1="cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    conda:
        kernel
    shell:
        "python ../../bin/run_cutadapt.py {template} 4 {raw_fastq_dir} ./ {wildcards.sample}_{wildcards.locus}"

rule write_BMK_config:
    input:
        r1="cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    output:
        config="BST_config/{sample}_{locus}.config.txt"
    shell:
        "python ../../bin/write_BMK_config.py {wildcards.sample} {wildcards.locus}"

rule extract_spatial_barcode:
    input:
        config="BST_config/{sample}_{locus}.config.txt"
    output:
        select_id="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.select_id",
        parsed_barcodes="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.bc_umi_read.tsv.id",
    conda:
        "BST-env"
    shell:
        "~/tools/BSTMatrix_v2.4.f.1/BSTMatrix -c {input.config} -s 1"

rule run_DARLIN_pipeline:
    input:
        select_id="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.select_id",
        parsed_barcodes="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.bc_umi_read.tsv.id",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    output:
        BST_output=directory("BST_output/{sample}_{locus}/02.Umi2Gene"),
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
        umi2gene='BST_output/{sample}_{locus}/02.Umi2Gene/out.umi_gene.tsv',
        notebook="BST_output/{sample}_{locus}/spatial_DARLIN-BMK.ipynb"
    shell:
        "papermill ../../QC/spatial_DARLIN-BMK.ipynb -k {kernel} {output.notebook} "
        "-p select_id_file {input.select_id} "
        "-p parsed_barcodes_file {input.parsed_barcodes} "
        "-p features_file {output.feat} "
        "-p umi2gene_file {output.umi2gene} "
        "-p dar_reads {input.r2} "
        "-p sample {wildcards.sample} "
        "-p output_dir {output.BST_output}"

rule generate_slim_fastq:
    input:
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv'
    output:
        r1="slim_fastq/{sample}_{locus}_R1.fastq.gz",
        r2="slim_fastq/{sample}_{locus}_R2.fastq.gz",
        notebook="BST_output/{sample}_{locus}/extra-gen_BMK_fastq.ipynb"
    shell:
        "papermill ../../QC/extra-gen_BMK_fastq.ipynb -k {kernel} {output.notebook} "
        "-p sample {wildcards.sample} -p locus {wildcards.locus} "
        "-p called_bc_file {input.feat} "
        "-p R1_file {output.r1} "
        "-p R2_file {output.r2} "
        "-p sc10xv3_cell_barcode_whitelist {CARLIN_dir}/cfg/10xV3_barcodes.txt.gz"

rule allele_calling:
    input:
        r1="slim_fastq/{sample}_{locus}_R1.fastq.gz",
        r2="slim_fastq/{sample}_{locus}_R2.fastq.gz"
    output:
        seq="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/Actaul_CARLIN_seq.txt",
        anno="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/AlleleAnnotations.txt"
    run:
        output_dir=config['data_dir']+'/DARLIN/'+config['template']+'_cutoff_override_1'
        input_dir=config['data_dir']+'/slim_fastq/'
        print(input_dir)

        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        read_cutoff_UMI_override=1
        read_cutoff_CB_override=1
        
        file_size = os.path.getsize(f'{input.r1}')*5/1000000000
        print(f"{input.r2}:   FileSize {file_size} G")
        requested_memory=5
        if requested_memory<20:
            requested_memory=20 # at least request 20G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 200G memory
        print(f"{wildcards.sample}_{wildcards.locus}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}_{wildcards.locus}',exist_ok=True)
        
        # do not run sbatch within this command
        command_1=f"bash {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample}_{wildcards.locus} {cfg_type} {template} {read_cutoff_UMI_override} {read_cutoff_CB_override} {requested_memory} {0} {CARLIN_max_run_time}"
        combined_command=f"{command_1}"
        print("Run on terminal directly")
        os.system(combined_command)

rule update_features:
    input:
        seq="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/Actaul_CARLIN_seq.txt",
        anno="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/AlleleAnnotations.txt",
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
    output:
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features_allele.tsv',
        done='BST_output/{sample}_{locus}/update_features.done'
    shell:
        "python ../../bin/update_features.py {input.seq} {input.anno} {input.feat} {output.feat} && "
        "touch {output.done}"

rule generate_level_matrix:
    input:
        config="BST_config/{sample}_{locus}.config.txt",
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
        umi2gene='BST_output/{sample}_{locus}/02.Umi2Gene/out.umi_gene.tsv',
        done='BST_output/{sample}_{locus}/update_features.done'
    output:
        done='BST_output/{sample}_{locus}/generate_matrix.done'
    conda:
        "BST-env"
    shell:
        "~/tools/BSTMatrix_v2.4.f.1/BSTMatrix -c {input.config} -s 3,4,5 && "
        "touch {output.done}"

rule finish:
    input:
        r2="slim_fastq/{sample}_{locus}_R2.fastq.gz",
        seq="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/Actaul_CARLIN_seq.txt",
        anno="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/AlleleAnnotations.txt",
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features_allele.tsv',
        done='BST_output/{sample}_{locus}/generate_matrix.done'
    output:
        done="BST_output/{sample}_{locus}/all.done"
    shell:
        "touch {output.done}"
