from pathlib import Path
from Bio import Align, SeqIO

records = list(SeqIO.parse("tp53_multi.fa", "fasta"))
id_a, seq_human = records[0].id.split()[0], str(records[0].seq)
id_b, seq_mouse = records[1].id.split()[0], str(records[1].seq)

# 1. Global
aligner_global = Align.PairwiseAligner()
aligner_global.mode = "global"
aligner_global.match_score = 1.0
aligner_global.mismatch_score = 0.0
aligner_global.open_gap_score = -1.0
aligner_global.extend_gap_score = -0.5

# 2. Local
aligner_local = Align.PairwiseAligner()
aligner_local.mode = "local"
aligner_local.match_score = 1.0
aligner_local.mismatch_score = 0.0
aligner_local.open_gap_score = -1.0
aligner_local.extend_gap_score = -0.5

aln_global = aligner_global.align(seq_human, seq_mouse)[0]
aln_local = aligner_local.align(seq_human, seq_mouse)[0]

print(f"Scor Aliniament Global: {aln_global.score:.1f}")
print(f"Scor Aliniament Local:  {aln_local.score:.1f}")

print("\n--- Fragment Aliniament Global (primele 100 bp) ---")
print(format(aln_global)[:200])

print("\n--- Fragment Aliniament Local (primele 100 bp) ---")
print(format(aln_local)[:200])