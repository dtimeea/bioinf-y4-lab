# Laboratorul 2 - Sequence Alignment

## Activitati realizate
- Am rulat demo01_pairwise.py.
- Am rulat demo02_distance_matrix.py.
- Am implementat si rulat Needleman-Wunsch.
- Am implementat si rulat Smith-Waterman.

## Date utilizate
Am folosit secventele TP53 umane descarcate din NCBI
in Lab 1, din data/work/dtimeea/lab01/my_tp53.fa.

| Secventa | Lungime (nt) |
|---|---:|
| OZ582015.1 | 359 |
| OZ582014.1 | 282 |
| OZ582013.1 | 291 |

## Rezultate
Pentru secventele OZ582015.1 si OZ582014.1:
- Needleman-Wunsch: scor -101.
- Smith-Waterman: scor 220.
- Demo global (100 nt): scor -9.
- Demo local (100 nt): scor 10.

Scorurile NW si SW nu sunt direct comparabile,
deoarece parametrii de punctaj sunt diferiti.

## Reflectie: Global vs. Local
Alinierea globala este potrivita pentru secvente
de lungimi similare si inrudite pe intreaga lungime.

Alinierea locala este potrivita pentru identificarea
regiunilor conservate in secvente divergente
sau de lungimi diferite.
