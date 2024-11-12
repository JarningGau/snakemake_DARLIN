import sys
import pandas as pd

seq_file = sys.argv[1]
anno_file = sys.argv[2]
feat_file = sys.argv[3]
out_file = sys.argv[4]

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("length of seq1 and seq2 are not equal!")
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def find_best_match(seq_listA, seq_listB):
    best_matches = {}
    for seqB in seq_listB:
        best_match = None
        min_mismatches = float('inf')
        for seqA in seq_listA:
            if len(seqA) == len(seqB):
                mismatches = hamming_distance(seqA, seqB)
                if mismatches < min_mismatches:
                    min_mismatches = mismatches
                    best_match = seqA
        if best_match is not None:
            best_matches[seqB] = best_match
    return best_matches


#### Main ####
feat_df = pd.read_csv(feat_file, sep='\t', header=None, names=['clone_id', 'clone_bc', 'type'])
feat_df = feat_df.drop('clone_bc', axis=1)
clone_barcodes = feat_df['clone_id'].tolist()

## generate new DataFrame with callable information
allele_df = pd.DataFrame({
    "callable_seq": [i.rstrip().replace('-', '') for i in open(seq_file).readlines()],
    "allele": [i.rstrip() for i in open(anno_file).readlines()]
})
correct_seqs = find_best_match(clone_barcodes, allele_df['callable_seq'])
allele_df['clone_id'] = allele_df['callable_seq'].apply(lambda x: correct_seqs[x])

data_final = pd.merge(feat_df, allele_df, how="left", on="clone_id")
data_final.pop('callable_seq')
data_final = data_final[['clone_id', 'allele', 'type']]

data_final.to_csv(out_file, index=0, sep='\t')