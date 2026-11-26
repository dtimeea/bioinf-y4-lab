
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def main():
    # TODO 1: Încărcăm fișierul multi-FASTA propriu
    fasta = Path("data/work/dtimeea/lab04/your_sequences.fasta")

    if not fasta.exists():
        raise FileNotFoundError(
            f"Nu există fișierul {fasta}. "
            "Adaugă mai întâi un fișier multi-FASTA cu cel puțin 3 secvențe aliniate."
        )

    # Citim secvențele aliniate
    alignment = AlignIO.read(str(fasta), "fasta")
    print(f"Secvențe încărcate: {len(alignment)}")

    if len(alignment) < 3:
        raise ValueError("Fișierul trebuie să conțină cel puțin 3 secvențe.")

    # TODO 2: Calculăm matricea de distanțe
    calculator = DistanceCalculator("identity")
    distances = calculator.get_distance(alignment)

    print("\nMatricea de distanțe:")
    print(distances)

    # TODO 3: Construim arborele Neighbor-Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distances)

    # TODO 4: Salvăm arborele în format Newick
    output_file = Path(
        "labs/04_phylogenetics/submissions/dtimeea/tree_dtimeea.nwk"
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output_file), "newick")
    print(f"\nArborele a fost salvat în: {output_file}")

    # TODO 5: Vizualizăm arborele în terminal
    print("\nArborele Neighbor-Joining:")
    Phylo.draw_ascii(tree)


if __name__ == "__main__":
    main()
