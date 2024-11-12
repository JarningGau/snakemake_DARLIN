import os
import sys

sample = sys.argv[1] # 'E14.5_embryo_dox_E10.5'
locus = sys.argv[2] # 'CA'

# This is the config file for BMK pipeline
config = f'''
### Globle parameters
## data 数据路径
FQ1	cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz
FQ2	cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz

## Flu info file 荧光解码文件路径
FLU	../image/{sample}-HE.txt

## AllheStat.py & CellSplit 组织识别及细胞分割
#HE染色图片路径、是否识别组织内部空白(1识别，0忽略，默认1)
HE	../image/{sample}-HE.tif
#INSIDE	1
#GRAY	200
#是否做细胞分割及荧光图片路径、颜色通道（颜色通道需要根据实际使用的进行选择0、1或2）
CellSplit	True
fluorescence	../image/{sample}-FL.tif
fluorescence_channl	   2 
#FLGRAY	15
#cells_npy	/path/to/cells/npyfile
#细胞分割参数
#YAML	/path/to/cell_split/parameter/file
#enhance	1

## Ref genome 参考基因组版本信息、STAR建库目录、gff/gtf文件及features文件路径
#GenomeVer	xxxx
# INDEX	/mnt/d/resource/star/mm10-2020-A/
# GFF		/mnt/d/resource/star/mm10-2020-A/genes.gtf
FEATURE	./BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv

## output    输出结果及输出文件前缀
OUTDIR	./BST_output/{sample}_{locus}
PREFIX	out

### Local parameters
## fastq2BcUmi    barcode版本类型(一般为V2版本)及barcode识别线程数
BCType	V2
BCThreads	8

## Umi2Gene       SART比对参数
#Sjdboverhang	100
#STARThreads	    8

## ENV      python和Rscript的路径，如不提供则使用系统环境中的版本(不提供请注释掉以下参数)
##PYTHON	/path/to/python/dir/
##Rscript	/path/to/Rscript/dir/
'''

if not os.path.exists('BST_config'):
    os.makedirs('BST_config')

if not os.path.exists(f'./BST_output/{sample}_{locus}'):
    os.makedirs(f'./BST_output/{sample}_{locus}')

with open(f'BST_config/{sample}_{locus}.config.txt', 'w') as f:
    f.write(config)