# Lab 08 - Machine Learning

## Activitati realizate

- Rulat demo01_sup_unsup.py.
- Rulat demo02_expr.py.
- Rulat demo03_logreg.py.
- Implementat ex01_sup_pipeline.py.
- Implementat ex02_logreg_vs_ref.py.
- Comparat Random Forest cu Logistic Regression.
- Comparat clusterele KMeans cu etichetele reale.

## Dataset

Am utilizat datasetul GSE96058, folosit si in laboratoarele 6 si 7.

Sursa:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058

Datele reprezinta expresia genica din probe de cancer mamar.

Am utilizat metadatele platformei GPL11154 si caracteristica
clinica ER status pentru clasificare.

Am exclus replicile tehnice si probele fara eticheta valida.

Din cele 2569 probe ER pozitive si 214 probe ER negative
disponibile am selectat aleatoriu cate 200 din fiecare clasa.

Datasetul final contine:

- 400 probe.
- 500 gene selectate aleatoriu.
- 200 probe ER_positive.
- 200 probe ER_negative.
- 0 valori lipsa.

Datele de expresie erau deja transformate.
Nu am aplicat o noua transformare logaritmica.

Fisier local:
data/work/dtimeea/lab08/expression_matrix_dtimeea.csv

## Configurarea modelelor

- Train: 80%, respectiv 320 probe.
- Test: 20%, respectiv 80 probe.
- Impartire stratificata.
- Random state: 42.
- Random Forest: 200 arbori.
- Logistic Regression: maximum 1000 iteratii.
- StandardScaler ajustat doar pe datele de antrenare.
- Aceeasi impartire train/test pentru ambele modele.

## Rezultatele demo-urilor

- Demo 01 Random Forest: accuracy 94.4%.
- Demo 02 Random Forest: accuracy 100%.
- Demo 03 Logistic Regression: accuracy 100%.

Pentru demo-urile 2 si 3 am corectat caile, am eliminat
coloana sample_id si am utilizat o impartire 50/50,
deoarece datasetul demonstrativ contine doar 6 probe.

Rezultatele de 100% provin dintr-un set de testare foarte mic.

## Rezultatele clasificarii

| Model | Accuracy | Macro F1 | Predictii corecte |
|---|---|---|---|
| Random Forest | 90% | 0.90 | 72/80 |
| Logistic Regression | 88.75% | 0.89 | 71/80 |

Random Forest a obtinut o acuratete usor mai mare,
cu o diferenta de 1.25 puncte procentuale.

## Feature importance

Primele cinci gene dupa importanta in Random Forest:

| Gena | Importanta |
|---|---|
| ULBP2 | 0.036796 |
| ZG16B | 0.031261 |
| KLHDC1 | 0.027506 |
| FAM174A | 0.027442 |
| BBS4 | 0.027338 |

Importanta caracteristicilor reflecta contributia lor
in model, nu demonstreaza cauzalitate biologica.

## KMeans

Am utilizat KMeans cu 2 clustere, fara etichete
in etapa de clustering.

| Label | Cluster 0 | Cluster 1 |
|---|---|---|
| ER_negative | 53 | 147 |
| ER_positive | 170 | 30 |

Clusterele prezinta o asociere cu etichetele reale,
dar separarea claselor nu este perfecta.

## Reflectie

Random Forest este un model neliniar, bazat pe arbori
de decizie, care poate surprinde interactiuni complexe
intre caracteristici.

Logistic Regression este un model liniar care estimeaza
probabilitatea apartenentei la o clasa.

In acest experiment, Random Forest a obtinut un rezultat
usor mai bun, dar diferenta este prea mica pentru a
concluziona ca este general superior.

Invatarea supravegheata utilizeaza etichetele cunoscute,
in timp ce clusteringul nesupravegheat identifica
grupuri pe baza similaritatii dintre probe.

Datasetul echilibrat artificial si selectia aleatorie
a celor 500 de gene reprezinta simplificari didactice.
Rezultatele nu constituie o validare clinica.

## Fisiere generate

- classification_report_dtimeea.txt
- confusion_rf_dtimeea.png
- feature_importance_dtimeea.csv
- cluster_crosstab_dtimeea.csv
- rf_vs_logreg_report_dtimeea.txt

Datele originale si matricea de expresie sunt pastrate
local in data/work/dtimeea/ si nu sunt incluse in Git.
