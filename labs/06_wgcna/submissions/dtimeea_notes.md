# Lab 06 - Gene Co-expression Networks

## Activitati realizate

- Rulat demo01_corr_threshold.py.
- Descarcat un dataset real RNA-Seq din GEO.
- Pregatit matricea de expresie genica.
- Implementat ex01_gce_networks.py.
- Calculat matricea de corelatie Spearman.
- Construit matricea de adiacenta si graful.
- Detectat modulele prin algoritmul Louvain.
- Generat modules_dtimeea.csv.

## Dataset

- Accession: GSE96058
- Sursa: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
- Date RNA-Seq din cancer mamar.
- Datasetul GSE96058 a fost ales in locul GSE115469,
  deoarece acesta din urma contine date din ficat uman.
- Datele originale sunt deja normalizate si transformate
  logaritmic, deci transformarea nu a fost repetata.

## Preprocesare

- Numar initial de randuri: 30865.
- Selectate primele 100 de coloane de expresie.
- Eliminate valorile lipsa si infinite.
- Agregate identificatoarele duplicate prin medie.
- Prag de varianta: 0.1.
- Gene dupa filtrare: 22989.
- Selectate 500 de gene cu varianta maxima.
- Matrice finala: 500 x 100.
- Subsetul este exploratoriu si nu este neaparat
  reprezentativ pentru intregul studiu.
- Nu a fost verificata prezenta replicilor tehnice
  in subsetul selectat.

## Parametrii retelei

- Metrica: Spearman.
- Prag de corelatie: 0.7.
- Corelatii pozitive.
- Matrice de adiacenta binara.
- Graf neorientat, fara auto-conexiuni.
- Algoritm: Louvain.
- Seed: 42.

## Rezultate

- Gene analizate: 500.
- Gene izolate eliminate: 485.
- Noduri: 15.
- Muchii: 31.
- Module detectate: 4.

| Modul | Numar gene |
|-------|------------|
| 1     | 9          |
| 2     | 2          |
| 3     | 2          |
| 4     | 2          |

Pragul de 0.7 a produs o retea rara, pastrand
doar genele cu asocieri pozitive puternice.

## Reflectie

Spre deosebire de clusteringul clasic, care grupeaza
obiectele pe baza similaritatii sau distantei,
analiza retelelor de co-expresie reprezinta genele
ca noduri si corelatiile dintre acestea ca muchii.

Algoritmul Louvain identifica module pe baza
structurii conexiunilor din graf. Genele din acelasi
modul pot avea comportamente de expresie similare,
dar corelatia nu demonstreaza o relatie cauzala
sau o functie biologica comuna.

## Fisiere rezultate

- ex07_gce_networks.py
- modules_dtimeea.csv

Datele de lucru sunt pastrate local in
data/work/dtimeea/lab06/ si nu sunt incluse in Git.
