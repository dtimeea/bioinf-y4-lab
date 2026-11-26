
"""
Exercițiu 04 — FASTQ QC pe date proprii

Calculează:
- numărul total de citiri
- lungimea medie a citirilor
- proporția bazelor N
- scorul Phred mediu

Salvează raportul în:
labs/03_formats&NGS/submissions/dtimeea/qc_report_dtimeea.txt
"""

import gzip
from pathlib import Path
from Bio import SeqIO

# Username GitHub
handle = "dtimeea"

in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")

out_report = Path(
    f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt"
)
out_report.parent.mkdir(parents=True, exist_ok=True)

# Deschidem fișierul FASTQ existent
if in_fastq_plain.exists():
    reader = SeqIO.parse(str(in_fastq_plain), "fastq")
elif in_fastq_gz.exists():
    reader = SeqIO.parse(gzip.open(in_fastq_gz, "rt"), "fastq")
else:
    raise FileNotFoundError(
        f"Nu am găsit fișierul FASTQ:\n"
        f"  {in_fastq_plain}\n"
        f"  {in_fastq_gz}\n"
        "Rulați mai întâi ex01_fetch_fastq.py."
    )

num_reads = 0
total_length = 0
total_n = 0
total_phred = 0
total_bases = 0

# Calculăm statisticile pentru fiecare citire
try:
    for record in reader:
        seq_str = str(record.seq)
        phred = record.letter_annotations["phred_quality"]

        num_reads += 1
        total_length += len(seq_str)
        total_n += seq_str.upper().count("N")
        total_bases += len(seq_str)
        total_phred += sum(phred)
finally:
    if in_fastq_gz.exists() and not in_fastq_plain.exists():
        reader.stream.close() if hasattr(reader, "stream") else None

# Calculăm mediile și proporția bazelor N
len_mean = total_length / num_reads if num_reads > 0 else 0.0
n_rate = total_n / total_bases if total_bases > 0 else 0.0
phred_mean = total_phred / total_bases if total_bases > 0 else 0.0

# Salvăm raportul
with open(out_report, "w", encoding="utf-8") as out:
    out.write(f"Reads: {num_reads}\n")
    out.write(f"Mean length: {len_mean:.2f}\n")
    out.write(f"N rate: {n_rate:.4f}\n")
    out.write(f"Mean Phred: {phred_mean:.2f}\n")

print(f"[OK] QC report -> {out_report.resolve()}")
