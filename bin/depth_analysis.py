import pysam, csv, glob, sys
from collections import defaultdict

meta_id = sys.argv[1]            # e.g. INF001
bam_files = glob.glob("*.bam")

# codon → AA
codon_to_aa = { … same codon map … }

depth = {}

for bam in bam_files:
    with pysam.AlignmentFile(bam, "rb") as bf:
        ref = bf.references[0]
        for col in bf.pileup(contig=ref):
            pos = col.pos + 1
            counts = dict.fromkeys(['A','T','C','G','N'], 0)
            for pr in col.pileups:
                if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                    b = pr.alignment.query_sequence[pr.query_position]
                    counts[b] = counts.get(b,0)+1
                else:
                    counts['N'] += 1
            total = sum(counts.values())
            ratios = {b: (counts[b]/total if total>0 else 0) for b in counts}
            depth[(bam,ref,pos)] = (counts, ratios)

# consensus & codon grouping
consensus = {}
codons = defaultdict(lambda: [None,None,None])

for (bam,ref,pos),(counts,ratios) in depth.items():
    # 1) consensus base
    maj = max(['A','T','C','G'], key=lambda b: counts[b])
    consensus[(bam,ref,pos)] = maj
    # 2) codon
    idx   = (pos-1)//3 + 1
    frame = (pos-1)%3
    codons[(bam,ref,idx)][frame] = maj

# write long CSV
out_fn = f"{meta_id}_long.csv"
with open(out_fn,'w',newline='') as fo:
    w=csv.writer(fo)
    w.writerow([
      'MetaID','BAM','Reference','Position','Base','Ratio',
      'CodonIndex','Frame','Codon','AA','RefPos'
    ])
    for (bam,ref,pos),(counts,ratios) in depth.items():
        idx   = (pos-1)//3 + 1
        frame = (pos-1)%3 + 1
        trip  = "".join(codons[(bam,ref,idx)][i] or 'N' for i in range(3))
        aa    = codon_to_aa.get(trip,'X')
        rp    = f"{ref}-{pos}"
        for base in ['A','T','C','G']:
            w.writerow([
              meta_id,
              bam,
              ref,
              pos,
              base,
              ratios[base],
              idx,
              frame,
              trip,
              aa,
              rp
            ])
print("Wrote",out_fn)
