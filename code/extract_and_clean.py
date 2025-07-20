#!/usr/bin/env python
"""
Reads seqdump.txt, creates:
  • cleaned.fasta   (MAFFT‑friendly IDs)
  • metadata.csv    (accession, isolate, country, year, label)
"""

import re, csv
from pathlib import Path
from Bio import SeqIO

IN  = Path("data/seqdump.txt")
OUT_FASTA = Path("results/cleaned.fasta")
OUT_META  = Path("results/metadata.csv")

# --- simple country dictionary (extend as needed) -------------
countries = [
    "Nigeria", "Kenya", "Cameroon", "Senegal", "USA",
    "China", "Brazil", "South_Africa", "India"
]
# --------------------------------------------------------------

records, meta_rows = [], []

for rec in SeqIO.parse(IN, "fasta"):

    header = rec.description
    accession = header.split()[0].split(":")[0]          # e.g. AJ583713.1
    isolate   = re.search(r"[Ii]solate[ _-]?([A-Za-z0-9.-]+)",
                          header)
    isolate   = isolate.group(1) if isolate else "NA"

    # crude year extraction: first 4‑digit number ≥ 1950
    year = re.search(r"(19\d{2}|20\d{2})", header)
    year = year.group(1) if year else "NA"

    country = "NA"
    for c in countries:
        if c.lower().replace("_", " ") in header.lower():
            country = c.replace("_", " ")
            break

    # build label: Accession|Country|Isolate|Year
    label = f"{accession}|{country}|{isolate}|{year}"

    # keep FASTA id simple for MAFFT/IQ‑TREE
    rec.id = accession
    rec.name = accession
    rec.description = ""        # strip long header
    records.append(rec)

    meta_rows.append(
        dict(accession=accession,
             isolate=isolate,
             country=country,
             year=year,
             label=label)
    )

# write cleaned FASTA
OUT_FASTA.parent.mkdir(exist_ok=True, parents=True)
SeqIO.write(records, OUT_FASTA, "fasta")
print(f"✅ cleaned FASTA → {OUT_FASTA}")

# write metadata CSV
with OUT_META.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh, fieldnames=["accession", "country", "isolate", "year", "label"]
    )
    writer.writeheader(); writer.writerows(meta_rows)
print(f"✅ metadata CSV  → {OUT_META}")
