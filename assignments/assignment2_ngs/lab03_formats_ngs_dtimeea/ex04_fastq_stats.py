import gzip
import os
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

HANDLE = "dtimeea"

# Cauta fisierul FASTQ: prioritate fisierul tau, apoi fallback pe sample
fastq_candidates = [
    Path(f"data/work/{HANDLE}/lab03/your_reads.fastq.gz"),
    Path(f"data/work/{HANDLE}/lab03/reads.fastq.gz"),
    Path("data/sample/your_reads.fastq.gz"),
    Path("data/sample/sample_reads.fastq.gz"),
    Path("data/sample/reads.fastq.gz"),
    Path("data/sample/tp53_reads.fastq.gz")
]

FASTQ_PATH = None
for p in fastq_candidates:
    if p.exists():
        FASTQ_PATH = p
        break

# Daca nu a fost gasit un fisier gz, cautam varianta necomprimata (.fastq)
if not FASTQ_PATH:
    for p in list(Path("data").rglob("*.fastq*")):
        FASTQ_PATH = p
        break

if not FASTQ_PATH or not FASTQ_PATH.exists():
    raise FileNotFoundError("Nu s-a gasit niciun fisier FASTQ in repo!")

print(f"Fisier FASTQ utilizat: {FASTQ_PATH}")

lengths = []
phred_means = []
n_bases_total = 0
total_bases = 0

open_fn = gzip.open if str(FASTQ_PATH).endswith(".gz") else open

with open_fn(FASTQ_PATH, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq_str = str(record.seq)
        seq_len = len(seq_str)
        lengths.append(seq_len)
        
        n_bases_total += seq_str.upper().count("N")
        total_bases += seq_len
        
        phred_scores = record.letter_annotations["phred_quality"]
        phred_means.append(float(np.mean(phred_scores)))

total_reads = len(lengths)
avg_length = float(np.mean(lengths)) if total_reads else 0.0
prop_n = (n_bases_total / total_bases) if total_bases else 0.0
overall_phred = float(np.mean(phred_means)) if phred_means else 0.0

# 1. Salvare raport QC
report_file = f"qc_report_{HANDLE}.txt"
with open(report_file, "w", encoding="utf-8") as rep:
    rep.write(f"Raport FASTQ QC — Handle: {HANDLE}\n")
    rep.write(f"Fisier analizat: {FASTQ_PATH.name}\n")
    rep.write(f"Numar total citiri: {total_reads}\n")
    rep.write(f"Lungime medie a citirilor: {avg_length:.2f} bp\n")
    rep.write(f"Proportie baze ambigue (N): {prop_n:.6f} ({prop_n*100:.4f}%)\n")
    rep.write(f"Scor Phred mediu per baza: {overall_phred:.2f}\n")

print(f"Raport salvat in {report_file}")

# 2. Generare grafic bonus (qc_plot_dtimeea.png)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.hist(lengths, bins=20, color="#4A90E2", edgecolor="black")
ax1.set_title("Distributia lungimii citirilor")
ax1.set_xlabel("Lungime (bp)")
ax1.set_ylabel("Frecventa")

ax2.hist(phred_means, bins=20, color="#50E3C2", edgecolor="black")
ax2.axvline(20, color="orange", linestyle="--", label="Prag Q20 (99%)")
ax2.axvline(30, color="green", linestyle="--", label="Prag Q30 (99.9%)")
ax2.set_title("Distributia scorurilor Phred medii")
ax2.set_xlabel("Scor Phred mediu")
ax2.set_ylabel("Frecventa")
ax2.legend(fontsize=8)

plt.tight_layout()
plot_file = f"qc_plot_{HANDLE}.png"
plt.savefig(plot_file, dpi=200)
plt.close()
print(f"Grafic salvat in {plot_file}")