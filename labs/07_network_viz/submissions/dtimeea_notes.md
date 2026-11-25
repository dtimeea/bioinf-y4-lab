# Lab 07 - Network Visualization and Hub Genes

## Activitati realizate

- Rulat demo01_viz.py.
- Corectat calea de export din demo.
- Refolosit matricea de expresie din Lab 6.
- Incarcat modules_dtimeea.csv.
- Reconstruit matricea de adiacenta.
- Vizualizat reteaua de co-expresie.
- Identificat genele hub pe baza gradului.
- Exportat imaginea PNG si tabelul CSV.

## Dataset si parametri

- Dataset: GSE96058.
- Matrice de expresie: 500 gene x 100 coloane.
- Corelatie Spearman pozitiva.
- Prag de corelatie: 0.7.
- Graf neorientat si neponderat.
- Module preluate din Lab 6.

## Metoda de vizualizare

Am utilizat spring_layout din NetworkX,
cu seed 42 si 200 de iteratii.

Nodurile sunt colorate in functie de modul,
iar dimensiunea lor este proportionala cu gradul.

Primele 10 gene dupa grad sunt evidentiate
prin contur negru si etichete.

## Rezultate

- Numar noduri: 15.
- Numar muchii: 31.
- Numar module: 4.

Rezultatele coincid cu reteaua din Lab 6.

## Gene hub

| Gene | Degree | Module |
|------|--------|--------|
| AF343406 | 8 | 1 |
| KC823651 | 8 | 1 |
| AF343414 | 7 | 1 |
| HM222546 | 7 | 1 |
| JC037840 | 7 | 1 |
| FJ840426 | 6 | 1 |
| X17117 | 6 | 1 |
| JN576489 | 5 | 1 |
| IGHV1-18 | 2 | 1 |
| AK124272 | 1 | 2 |

Cele mai conectate gene sunt AF343406 si KC823651,
ambele avand gradul 8.

Clasamentul este realizat la nivelul intregului graf.
Gradul ridicat nu demonstreaza automat importanta
biologica sau implicarea cauzala intr-o boala.

## Reflectie

Vizualizarea aduce avantaje fata de analiza numerica
din Lab 6 deoarece permite observarea directa a
structurii retelei, a modulelor si a conexiunilor
dintre gene.

Culorile faciliteaza diferentierea modulelor,
iar dimensiunea nodurilor ajuta la identificarea
genelor cu grad ridicat.

Reprezentarea grafica completeaza matricea de
corelatie si rezultatele numerice, facilitand
interpretarea structurii retelei.

Reteaua este exploratorie si depinde de subsetul
de expresie si de pragul selectat in Lab 6.

## Fisiere rezultate

- ex01_network_viz.py
- network_dtimeea.png
- hubs_dtimeea.csv
