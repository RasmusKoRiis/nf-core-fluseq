import pysam
import csv
import glob
import sys
from collections import defaultdict

# Usage: python depth_analysis.py META_ID
meta_id = sys.argv[1]
bam_files = glob.glob("*.bam")

# Full codon → AA map
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

# Step 1: collect depth & ratios per (bam, ref, pos)
depth = {}
for bam in bam_files:
    with pysam.AlignmentFile(bam, "rb") as bf:
        ref = bf.references[0]
        for col in bf.pileup(contig=ref):
            pos = col.pos + 1
            # counts
            counts = {'A':0,'T':0,'C':0,'G':0,'N':0}
            for pr in col.pileups:
                if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                    b = pr.alignment.query_sequence[pr.query_position]
                    counts[b] = counts.get(b,0) + 1
                else:
                    counts['N'] += 1
            total = sum(counts.values())
            ratios = {b: (counts[b]/total if total>0 else 0) for b in counts}
            depth[(bam,ref,pos)] = (counts, ratios)

# Step 2: determine consensus base and group codons
consensus = {}
codons = defaultdict(lambda: [None,None,None])

for (bam,ref,pos),(counts,ratios) in depth.items():
    # consensus base (highest count among A,T,C,G)
    maj = max(['A','T','C','G'], key=lambda x: counts[x])
    consensus[(bam,ref,pos)] = maj
    # assign into codon frames
    idx   = (pos-1)//3 + 1
    frame = (pos-1)%3
    codons[(bam,ref,idx)][frame] = maj

# Step 3: write long CSV
out_fn = f"{meta_id}_long.csv"
with open(out_fn, 'w', newline='') as fo:
    w = csv.writer(fo)
    w.writerow([
        'MetaID','BAM','Reference','Position',
        'A_Count','T_Count','C_Count','G_Count','N_Count',
        'A_Ratio','T_Ratio','C_Ratio','G_Ratio','N_Ratio',
        'CodonIndex','Frame','Codon','AA','RefPos','Base','Ratio'
    ])
    for (bam,ref,pos),(counts,ratios) in depth.items():
        idx   = (pos-1)//3 + 1
        frame = ((pos-1)%3) + 1
        trip  = "".join(codons[(bam,ref,idx)][i] or 'N' for i in range(3))
        aa    = codon_to_aa.get(trip, 'X')
        refpos= f"{ref}-{pos}"
        # emit one row per base
        for base in ['A','T','C','G']:
            w.writerow([
                meta_id,
                bam,
                ref,
                pos,
                counts['A'], counts['T'], counts['C'], counts['G'], counts['N'],
                ratios['A'], ratios['T'], ratios['C'], ratios['G'], ratios['N'],
                idx,
                frame,
                trip,
                aa,
                refpos,
                base,
                ratios[base]
            ])

print(f"Wrote {out_fn}")
