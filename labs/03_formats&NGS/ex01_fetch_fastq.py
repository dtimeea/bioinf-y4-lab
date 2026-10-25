
import argparse
import json
import urllib.parse
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def main():
    parser = argparse.ArgumentParser(
        description="Download a FASTQ file from ENA."
    )
    parser.add_argument(
        "accession",
        help="ENA/SRA run accession, e.g. ERR164409 or SRR..."
    )
    parser.add_argument(
        "--handle",
        required=True,
        help="Your GitHub username"
    )
    args = parser.parse_args()

    # Create the output directory
    output_dir = Path("data") / "work" / args.handle / "lab03"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Query ENA for the FASTQ URL
    params = urllib.parse.urlencode({
        "accession": args.accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp",
        "format": "json",
    })

    url = f"{ENA_API}?{params}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "bioinf-y4-lab-student"}
    )

    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            rows = json.load(response)
    except Exception as exc:
        raise SystemExit(f"Could not query ENA: {exc}")

    if not rows or not rows[0].get("fastq_ftp"):
        raise SystemExit(
            f"No FASTQ URL found for accession {args.accession}"
        )

    # Use the first FASTQ URL returned by ENA
    fastq_url = rows[0]["fastq_ftp"].split(";")[0].strip()

    # Ensure the URL has a valid protocol
    if fastq_url.startswith("ftp://"):
        fastq_url = "https://" + fastq_url[len("ftp://"):]
    elif not fastq_url.startswith(("http://", "https://")):
        fastq_url = "https://" + fastq_url

    # Download to a temporary file first
    output_path = output_dir / "your_reads.fastq.gz"
    temp_path = output_dir / "your_reads.fastq.gz.part"

    try:
        request = urllib.request.Request(
            fastq_url,
            headers={"User-Agent": "bioinf-y4-lab-student"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            with open(temp_path, "wb") as output_file:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output_file.write(chunk)

        temp_path.replace(output_path)

    except Exception as exc:
        temp_path.unlink(missing_ok=True)
        raise SystemExit(f"Download failed: {exc}")

    print(f"Downloaded: {output_path.resolve()}")


if __name__ == "__main__":
    main()
