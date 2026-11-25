# Laboratorul 4 - Filogenetica

## Activitati realizate

- Am rulat demo01_distance_matrix.py.
- Am descarcat 7 secvente proteice p53 din UniProt.
- Am realizat alinierea multipla folosind MAFFT v7.505.
- Am calculat matricea de distante cu Biopython.
- Am construit arborele Neighbor-Joining.
- Am salvat arborele in format Newick.

## Dataset utilizat

Sursa: https://www.uniprot.org/

Am utilizat secvente ale proteinei p53 din 7 specii:

| Specie | Accession UniProt | Lungime initiala |
|---|---|---:|
| Human | P04637 | 393 |
| Mouse | P02340 | 390 |
| Rat | P10361 | 391 |
| Dog | Q29537 | 381 |
| Cow | P67939 | 386 |
| Chicken | P10360 | 367 |
| Zebrafish | P79734 | 373 |

Secventele au fost descarcate din intrari UniProt reviewed.

Fisier initial:
data/work/dtimeea/lab04/p53_raw.fasta

Fisier aliniat:
data/work/dtimeea/lab04/your_sequences.fasta

Aliniamentul contine 7 secvente si 419 pozitii.

## Rezultate

Matricea de distante a fost calculata folosind
DistanceCalculator("identity") din Biopython.

Cea mai mica distanta este intre Mouse si Rat (0.121718).
Cea mai mare distanta este intre Chicken si Zebrafish
(0.534606).

Arborele Neighbor-Joining grupeaza Mouse si Rat,
evidentiind similaritatea secventelor lor.

Arborele obtinut este neînrădăcinat si reprezinta
o estimare a relatiilor evolutive pe baza proteinei p53,
nu o reconstructie definitiva a evolutiei speciilor.

Rezultatul a fost salvat in tree_dtimeea.nwk.

## Reflectie

O matrice de distante prezinta diferentele dintre
perechile de secvente, in timp ce arborele filogenetic
ofera o reprezentare a relatiilor evolutive estimate.

Arborele permite identificarea grupurilor de secvente
apropiate, a ramificatiilor si a distantelor evolutive
estimate prin lungimea ramurilor. Astfel, relatiile
dintre secvente sunt mai usor de interpretat decat
prin simpla analiza a valorilor din matrice.
