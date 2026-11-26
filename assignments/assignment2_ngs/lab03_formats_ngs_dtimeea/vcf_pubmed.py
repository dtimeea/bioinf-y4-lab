import os
from pathlib import Path
from Bio import Entrez

Entrez.email = "dobrean.timeea@student.upt.ro"
Entrez.tool = "bioinf_y4_lab3"

HANDLE = "dtimeea"
OUT_FILE = f"variants_{HANDLE}.txt"

# 1. Cautam fisierul VCF (prioritate folderul curent)
vcf_path = Path("variants.vcf")

if not vcf_path.exists() or vcf_path.stat().st_size == 0:
    for p in [Path(f"data/work/{HANDLE}/lab03/variants.vcf"), Path("data/sample/variants.vcf")]:
        if p.exists():
            vcf_path = p
            break

# 2. Daca nu exista sau nu are variante, cream variants.vcf cu date reale TP53
need_create = True
if vcf_path.exists():
    with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                need_create = False
                break

if need_create:
    vcf_path = Path("variants.vcf")
    with open(vcf_path, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr17\t7674220\trs1042522\tC\tG\t100\tPASS\tGENE=TP53\n")
        f.write("chr17\t7675088\t.\tC\tT\t100\tPASS\tGENE=TP53\n")
    print("S-a creat automat variants.vcf cu 2 variante TP53.")

print(f"Fisier VCF utilizat: {vcf_path.resolve()}")

# 3. Parsare variante din VCF
variants = []
with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 5:
            chrom = parts[0]
            pos = parts[1]
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            variants.append((chrom, pos, var_id, ref, alt))
            if len(variants) >= 3:
                break

print(f"Numar variante extrase: {len(variants)}")

# 4. Cautare pe PubMed si salvare in variants_dtimeea.txt
with open(OUT_FILE, "w", encoding="utf-8") as out:
    for chrom, pos, var_id, ref, alt in variants:
        if var_id != "." and var_id.lower().startswith("rs"):
            query = var_id
        else:
            chrom_clean = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
            query = f"{chrom_clean}:{pos} AND TP53"
            
        print(f"Cautare PubMed -> '{query}'")
        
        try:
            handle_search = Entrez.esearch(db="pubmed", term=query, retmax=3)
            res = Entrez.read(handle_search)
            handle_search.close()
            pmids = res.get("IdList", [])
        except Exception as e:
            print(f"Eroare PubMed: {e}")
            pmids = []
        
        out.write(f"Varianta: Chrom={chrom}, Pos={pos}, ID={var_id}, Ref={ref}, Alt={alt}\n")
        out.write(f"Interogare utilizata: {query}\n")
        out.write(f"Articole identificate (PMID): {', '.join(pmids) if pmids else 'Niciun articol gasit'}\n\n")

print(f"Gata! Verifica acum continutul fisierului {OUT_FILE}.")