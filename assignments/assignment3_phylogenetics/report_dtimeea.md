# Assignment 3 - Phylogenetics

Autor: Dobrean Timeea
GitHub: dtimeea

## 1. Dataset si matricea de distante

Am utilizat 11 secvente proteice p53 din UniProt,
extinzand datasetul reconstruit in Lab 4 de la 7
la 11 specii.

Accession-urile utilizate sunt:
P04637, P02340, P10361, Q29537, P67939,
P10360, P79734, P41685, Q9TUB2, Q95330, P56424.

Secventele au fost aliniate cu MAFFT v7.505,
folosind strategia L-INS-i.
Alinierea contine 423 de pozitii.

Matricea de distante 11 x 11 a fost calculata
cu Biopython DistanceCalculator("identity")
si salvata in distance_matrix_dtimeea.csv.

Perechea cea mai apropiata este om-macac,
cu distanta 0.0402.

Distanta maxima este 0.5248, observata pentru
perechile soarece-peste zebra si gaina-peste zebra.

## 2. Arbore Neighbor-Joining

Arborele a fost construit cu Biopython si salvat
in format Newick, in tree_dtimeea.nwk.

Principalele grupari identificate sunt:

- Om si macac: primate.
- Soarece si sobolan: rozatoare.
- Pisica si caine: carnivore.
- Porc si vaca: artiodactile.
- Gaina si peste zebra: ramuri divergente
  fata de mamifere.

Gruparile mamiferelor sunt compatibile cu
relatiile taxonomice cunoscute.

Arborele este neinradacinat si reflecta
similaritatea proteinei p53, nu intreaga
istorie evolutiva a speciilor.

## 3. Comparatie cu Clustal Omega online

Acelasi dataset a fost incarcat in Clustal Omega
prin platforma EMBL-EBI.

Alinierea MAFFT contine 423 de pozitii,
iar alinierea Clustal Omega contine 413.

Am verificat ca secventele fara gap-uri
coincid exact cu secventele originale.

Am construit un al doilea arbore NJ folosind
alinierea Clustal Omega si aceeasi metoda
de calcul a distantelor.

Rezultatele comparatiei topologice:

- Bipartitii MAFFT: 8.
- Bipartitii Clustal Omega: 8.
- Bipartitii comune: 8.
- Diferente topologice: 0.

Cele doua metode produc aceeasi topologie
pentru datasetul analizat, in ciuda diferentelor
de lungime dintre alinieri.

## 4. Regiuni conservate

In alinierea Clustal Omega au fost identificate
132 de pozitii complet conservate din 413,
reprezentand aproximativ 31.96%.

Regiunile consecutive de minimum 5 aminoacizi
complet conservati sunt:

| Pozitii MSA | Secventa | Lungime |
|---|---|---:|
| 245-259 | MCNSSCMGGMNRRPI | 15 |
| 278-290 | FEVRVCACPGRDR | 13 |
| 178-188 | EVVRRCPHHER | 11 |
| 201-208 | HLIRVEGN | 8 |
| 124-129 | GTAKSV | 6 |

Pozitiile sunt coordonate ale alinierii,
nu pozitii individuale ale proteinei umane.

Conservarea acestor regiuni sugereaza existenta
unor constrangeri structurale si functionale.
Proteina p53 are un rol important in reglarea
ciclului celular si raspunsul la deteriorarea ADN.

Pozitiile complet conservate sunt marcate cu *
in conserved_alignment_dtimeea.txt.

## 5. Bonus - Vizualizare

Arborele filogenetic a fost reprezentat grafic
cu matplotlib in tree_dtimeea.png.

Comparatia grafica dintre cei doi arbori este
disponibila in tree_comparison_dtimeea.png.

## 6. Limitari si concluzii

Distanta identity include pozitiile cu gap-uri,
inclusiv potrivirile gap-gap, si nu reprezinta
un model evolutiv complex.

Nu a fost efectuata o analiza bootstrap.
Topologia identica nu demonstreaza automat
corectitudinea tuturor relatiilor evolutive.

Pentru datasetul analizat, principalele grupari
sunt stabile fata de alegerea metodei MSA.

## Surse

- https://www.uniprot.org/
- https://mafft.cbrc.jp/alignment/software/
- https://www.ebi.ac.uk/jdispatcher/msa/clustalo
- https://biopython.org/

## Utilizarea AI

ChatGPT a fost utilizat pentru asistenta la
scrierea codului, interpretarea rezultatelor
si structurarea raportului.
Scripturile au fost executate in GitHub Codespaces,
iar alinierea Clustal Omega a fost realizata online.
