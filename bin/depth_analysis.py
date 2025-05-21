import pysam
import csv
import glob
import os
import sys
from collections import defaultdict

# Usage: python depth_annotate_no_fasta.py META_ID
meta_id   = sys.argv[1]
bam_files = glob.glob("*.bam")

# codon->AA table
codon_to_aa = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# will hold: {(bam,ref,pos): {counts...}}
depth_dict = {}

for bam in bam_files:
    with pysam.AlignmentFile(bam, "rb") as bf:
        ref = bf.references[0]
        for col in bf.pileup(contig=ref):
            pos = col.pos + 1
            counts = {'A':0,'T':0,'C':0,'G':0,'N':0}
            for pr in col.pileups:
                if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                    base = pr.alignment.query_sequence[pr.query_position]
                    counts[base] = counts.get(base,0) + 1
                else:
                    counts['N'] += 1
            total = sum(counts.values())
            # ratios
            ratios = {b: (counts[b]/total if total>0 else 0) for b in counts}
            depth_dict[(os.path.basename(bam), ref, pos)] = {
                **counts,
                **{f"{b}_Ratio": ratios[b] for b in ratios},
                "TotalDepth": total
            }

# Now build consensus‐base per position
consensus = {}
for key, vals in depth_dict.items():
    # pick base with max count (ties broken arbitrarily)
    b, r, p = key
    # exclude N when picking consensus?
    pick = max(['A','T','C','G'], key=lambda x: vals.get(x,0))
    consensus[key] = pick

# Group into codons
# {(bam,ref,codon_idx): [bases at pos1,pos2,pos3]}
codon_bases = defaultdict(lambda: [None,None,None])

for (bam,ref,pos), base in consensus.items():
    idx = (pos-1)//3 + 1
    frame = (pos-1)%3             # 0,1,2 → position inside the codon
    codon_bases[(bam,ref,idx)][frame] = base

# Build final CSV rows
out_rows = []
for (bam,ref,pos), vals in depth_dict.items():
    idx = (pos-1)//3 + 1
    bases = codon_bases[(bam,ref,idx)]
    # if any base is missing (e.g. low coverage), fill with 'N'
    codon = "".join(b if b else 'N' for b in bases)
    aa    = codon_to_aa.get(codon, 'X')
    row = [
        meta_id,
        bam,
        ref,
        pos,
        vals["TotalDepth"],
        vals["A"], vals["T"], vals["C"], vals["G"], vals["N"],
        vals["A_Ratio"], vals["T_Ratio"], vals["C_Ratio"], vals["G_Ratio"], vals["N_Ratio"],
        idx,
        codon,
        aa
    ]
    out_rows.append(row)

# Write annotated file
out_fn = f"{meta_id}_annotated_depth.csv"
with open(out_fn, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow([
        'MetaID','BAM','Reference','Position','TotalDepth',
        'A_Count','T_Count','C_Count','G_Count','N_Count',
        'A_Ratio','T_Ratio','C_Ratio','G_Ratio','N_Ratio',
        'CodonIndex','Codon','AA'
    ])
    w.writerows(out_rows)

print(f"✅ Wrote {out_fn} — includes ratios, consensus codon & AA")
