import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO

FASTA_PATH = Path("tp53_multi.fa")
records = list(SeqIO.parse(FASTA_PATH, "fasta"))

ids = [rec.id.split()[0] for rec in records]
seqs = [str(rec.seq) for rec in records]


def p_distance(s1: str, s2: str) -> float:
  min_len = min(len(s1), len(s2))
  mismatches = sum(c1 != c2 for c1, c2 in zip(s1[:min_len], s2[:min_len]))
  return mismatches / min_len


n = len(ids)
dist_matrix = np.full((n, n), np.nan)

for i, j in itertools.combinations(range(n), 2):
  dist_matrix[i, j] = p_distance(seqs[i], seqs[j])

df = pd.DataFrame(dist_matrix, index=ids, columns=ids)
print("\n--- TASK 1: Matrice p-distance (triunghi superior) ---")
print(df.round(4).to_string())

df.round(4).to_csv("distance_matrix.csv")
print("\nSalvat cu succes in distance_matrix.csv!")