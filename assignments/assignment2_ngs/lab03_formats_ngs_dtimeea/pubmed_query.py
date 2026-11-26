import os
from Bio import Entrez, Medline

Entrez.email = "dobrean.timeea@student.upt.ro"
Entrez.tool = "bioinf_y4_lab3"

HANDLE = "dtimeea"
OUT_FILE = f"pubmed_{HANDLE}.txt"
QUERY = "TP53 AND cancer"
MAX_RESULTS = 5

print(f"Cautare pe PubMed: {QUERY}")

# Cautare ID-uri articole
search_handle = Entrez.esearch(db="pubmed", term=QUERY, retmax=MAX_RESULTS, sort="pub_date")
search_record = Entrez.read(search_handle)
search_handle.close()

id_list = search_record.get("IdList", [])
print(f"Articole identificate: {id_list}")

# Descarcare detalii
fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
records = list(Medline.parse(fetch_handle))
fetch_handle.close()

# Salvare in fisier text
with open(OUT_FILE, "w", encoding="utf-8") as out:
    for idx, rec in enumerate(records, start=1):
        pmid = rec.get("PMID", "N/A")
        title = rec.get("TI", "Fara titlu")
        authors = ", ".join(rec.get("AU", ["Autori nespecificati"]))
        abstract = rec.get("AB", "Rezumat indisponibil.")
        
        out.write(f"--- Articol {idx} (PMID: {pmid}) ---\n")
        out.write(f"Titlu: {title}\n")
        out.write(f"Autori: {authors}\n")
        out.write(f"Rezumat:\n{abstract}\n\n")

print(f"Fisier creat cu succes: {OUT_FILE}")