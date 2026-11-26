# Laborator 5 – Clustering în Bioinformatică

## Rezultate și interpretare

În cadrul acestui laborator am aplicat trei metode de clustering pe dataset-ul Breast Cancer Wisconsin: Hierarchical Clustering, K-means și DBSCAN. Înainte de aplicarea algoritmilor, am preprocesat și standardizat datele.

Prin metoda Hierarchical Clustering am obținut o dendrogramă care evidențiază relațiile de similaritate dintre observații.

Algoritmul K-means, cu K=2, a împărțit cele 569 de observații în două clustere, unul cu 381 de elemente și celălalt cu 188 de elemente.

Algoritmul DBSCAN, folosind eps=1.5 și min_samples=5, a identificat un singur cluster cu 19 elemente, iar celelalte 550 de observații au fost considerate zgomot.

Consider că metoda K-means este cea mai potrivită pentru scopul acestui exercițiu, deoarece permite împărțirea tuturor observațiilor în două grupuri. În cazul DBSCAN, rezultatele sunt influențate de alegerea parametrilor și de dimensionalitatea datelor.

## Reflecție

Atât clustering-ul, cât și arborii filogenetici permit identificarea relațiilor dintre datele biologice pe baza similarităților. Clustering-ul grupează observațiile în funcție de caracteristicile lor, în timp ce arborii filogenetici urmăresc evidențierea relațiilor evolutive dintre secvențe sau organisme.

Cele două metode pot fi utilizate împreună pentru a înțelege mai bine relațiile biologice, pentru identificarea familiilor de gene sau pentru descoperirea unor subtipuri de boli. Totuși, similaritatea statistică nu implică întotdeauna o relație evolutivă directă.