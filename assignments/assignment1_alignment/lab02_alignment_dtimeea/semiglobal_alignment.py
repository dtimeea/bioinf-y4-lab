from Bio import Align, SeqIO

records = list(SeqIO.parse("tp53_multi.fa", "fasta"))
seq_human = str(records[0].seq)
seq_mouse = str(records[1].seq)

aligner_semi = Align.PairwiseAligner()
aligner_semi.mode = "global"
aligner_semi.match_score = 1.0
aligner_semi.mismatch_score = 0.0
aligner_semi.open_gap_score = -1.0
aligner_semi.extend_gap_score = -0.5

# Fara penalizare la capete:
aligner_semi.target_end_gap_score = 0.0
aligner_semi.query_end_gap_score = 0.0

aln_semi = aligner_semi.align(seq_human, seq_mouse)[0]
print(f"Scor Aliniament Semiglobal: {aln_semi.score:.1f}")
print(format(aln_semi)[:200])