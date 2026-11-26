from pathlib import Path
from io import StringIO
import requests
from Bio import SeqIO

folder = Path(__file__).resolve().parent
source = folder.parents[1] / "data/work/dtimeea/lab04/p53_raw.fasta"

records = list(SeqIO.parse(source, "fasta"))

species = {
    "Cat": 9685,
    "Pig": 9823,
    "Horse": 9796,
    "Rabbit": 9986,
    "Macaque": 9544
}

url = "https://rest.uniprot.org/uniprotkb/search"

for name, taxid in species.items():
    query = (
        f"gene_exact:TP53 AND organism_id:{taxid}"
        " AND reviewed:true"
    )

    response = requests.get(
        url,
        params={
            "query": query,
            "format": "fasta",
            "size": 10
        },
        timeout=60
    )
    response.raise_for_status()

    candidates = list(
        SeqIO.parse(StringIO(response.text), "fasta")
    )

    candidates = [
        r for r in candidates
        if 350 <= len(r.seq) <= 420
        and "Fragment" not in r.description
    ]

    if not candidates:
        print(f"{name}: nicio secventa adecvata")
        continue

    record = candidates[0]
    accession = record.id.split("|")[1]
    record.id = f"{name}_{accession}"
    record.description = ""

    records.append(record)
    print(f"{record.id}: {len(record.seq)} aminoacizi")

output = folder / "p53_multi_dtimeea.fasta"
SeqIO.write(records, output, "fasta")

print(f"\nTotal secvente: {len(records)}")
print(f"Fisier salvat: {output}")
