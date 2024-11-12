import os
import sys

templete = sys.argv[1] # {Tigre_2022_v2, Rosa_v2, cCARLIN}, long_primer_set: {Tigre_2022,Rosa,cCARLIN}
cores = int(sys.argv[2]) # 8
raw_fastq_dir = sys.argv[3]
data_dir = sys.argv[4]
sample = sys.argv[5]

primer_len = 15

PRIMERS = {
    "CA" : {
        "5prime": "AGCTGTACAAGTAAGCGGC",
        "3prime": "AGAATTCTAACTAGAGCTCGCTGATCAGCCTCGACTGTGCCTTCT"
    },
    "RA" : {
        "5prime":"GTACAAGTAAAGCGGCC",
        "3prime":"GTCTGCTGTGTGCCTTCTAGTT"
    },
    "TA" : {
        "5prime": "TCGGTACCTCGCGAA",
        "3prime": "GTCTTGTCGGTGCCTTCTAGTT"
    }
}

if templete.startswith('cCARLIN'):
    prime3 = PRIMERS['CA']['3prime']
    prime5 = PRIMERS['CA']['5prime']
elif templete.startswith('Rosa'):
    prime3 = PRIMERS['RA']['3prime']
    prime5 = PRIMERS['RA']['5prime']
elif templete.startswith('Tigre'):
    prime3 = PRIMERS['TA']['3prime']
    prime5 = PRIMERS['TA']['5prime']

prime3 = prime3[:primer_len] if len(prime3) > primer_len else prime3
prime5 = prime5[-primer_len:] if len(prime5) > primer_len else prime5

os.path.exists(f'{data_dir}/cutadapt') or os.makedirs(f'{data_dir}/cutadapt')
in_fq1 = f"{raw_fastq_dir}/{sample}_R1.fastq.gz"
in_fq2 = f"{raw_fastq_dir}/{sample}_R2.fastq.gz"
tmp_fq1 = f"{data_dir}/cutadapt/{sample}_R1.tmp.fastq.gz"
tmp_fq2 = f"{data_dir}/cutadapt/{sample}_R2.tmp.fastq.gz"
out_fq1 = f"{data_dir}/cutadapt/{sample}_R1.trimmed.fastq.gz"
out_fq2 = f"{data_dir}/cutadapt/{sample}_R2.trimmed.fastq.gz"

## Extract DARLIN sequence (Reads2)
## -j threads
## -A Adapter | 3' adapter to be removed from second read in a pair [R2].
## -G Adapter | 5' adapter to be removed from second read in a pair [R2].
## -q cutoff | Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. 
##           | Set to 0 means do not perform this trimming.
## -m cutoff | Minimum length after trimming.
## --discard-untrimmed | [Important!!] Discard reads that do not contain an adapter.
cmd1 = f"cutadapt -j {cores} -A {prime3} -q 0 -m 10 --discard-untrimmed -o {tmp_fq1} -p {tmp_fq2} {in_fq1} {in_fq2}"
cmd2 = f"cutadapt -j {cores} -G {prime5} -q 0 -m 10 --discard-untrimmed -o {out_fq1} -p {out_fq2} {tmp_fq1} {tmp_fq2}"

os.system(cmd1)
if os.system(cmd2) == 0:  # Check if cmd2 completed successfully
    os.system(f"rm {tmp_fq1} {tmp_fq2}")

