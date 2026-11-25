# Laboratorul 3 - Formate si NGS

## Activitati realizate

- Am rulat demo01_fastq_qc.py.
- Am rulat demo02_mapping_toy.py.
- Am implementat si rulat ex01_fetch_fastq.py.
- Am implementat si rulat ex02_fastq_stats.py.
- Am salvat raportul QC si scripturile in submissions/dtimeea.

## Sursa datelor

- Accession: SRR6109280
- Organism: Homo sapiens
- Experiment: CRISPR-DS Sequencing of TP53
- Platforma: Illumina MiSeq
- Sursa: https://www.ebi.ac.uk/ena/browser/view/SRR6109280
- Fisier utilizat: SRR6109280_1.fastq.gz

Fisierul a fost descarcat prin API-ul ENA si salvat in
data/work/dtimeea/lab03/your_reads.fastq.gz.

## Rezultate QC

| Indicator | Valoare |
|---|---:|
| Numar de citiri | 133195 |
| Lungime medie | 280.00 |
| N-rate | 0.0000 |
| Scor Phred mediu | 35.24 |

Scorul Phred mediu indica o calitate ridicata a citirilor.
Proportia bazelor N este foarte mica, fiind rotunjita la 0.0000.

## Reflectie

Verificarea calitatii datelor este esentiala inainte de
analiza variantelor deoarece erorile de secventiere pot
fi interpretate gresit drept mutatii reale. Citirile de
calitate slaba pot afecta alinierea si identificarea
variantelor genetice. Prin QC putem identifica problemele
datelor si reduce riscul rezultatelor fals pozitive.
