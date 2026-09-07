"""Experimental CSQ renderer for numeric dbNSFP values (no reference lookup).

Models the per-field formatting hook currently missing from vepyr. Uses only
the native numeric VCF output; never reads original plugin cells or VEP values.
The matching Rust policies are exhaustively checked by probe_numeric_sql.py.
"""

from functools import lru_cache
import gzip
import hashlib

FIXED = {
    "MetaSVM_score": 4,
    "MetaLR_score": 4,
    "phyloP100way_vertebrate": 6,
    "phastCons100way_vertebrate": 6,
}
FIELDS = {*FIXED, "GERP++_RS", "CADD_phred"}


@lru_cache(maxsize=100000)
def render(field, value):
    if not value:
        return value
    number = float(value)
    if field in FIXED:
        return f"{number:.{FIXED[field]}f}"
    if field == "CADD_phred":
        scale = 3 if number < 10 else 2 if number < 20 else 1 if number < 30 else 0
        return f"{number:.{scale}f}"
    if number and abs(number) < 0.001:
        mantissa, exponent = f"{number:.2e}".split("e")
        mantissa = mantissa.rstrip("0").rstrip(".")
        if "." not in mantissa:
            mantissa += ".0"
        return f"{mantissa}E{int(exponent)}"
    return value if "." in value else value + ".0"


def render_vcf(source, dest):
    indexes = []
    digest = hashlib.md5()
    records = changed = 0
    with gzip.open(source, "rb") as src, gzip.open(dest, "wb", compresslevel=1) as out:
        for raw in src:
            if raw.startswith(b"#"):
                if raw.startswith(b"##INFO=<ID=CSQ,"):
                    fields = (
                        raw.decode().split("Format: ", 1)[1].split('"', 1)[0].split("|")
                    )
                    indexes = [
                        (i, name) for i, name in enumerate(fields) if name in FIELDS
                    ]
                    assert len(indexes) == len(FIELDS), indexes
                out.write(raw)
                continue
            assert indexes
            records += 1
            line = raw.decode().rstrip("\n").split("\t")
            info = line[7].split(";")
            for i, entry in enumerate(info):
                if not entry.startswith("CSQ="):
                    continue
                csqs = entry[4:].split(",")
                for j, csq in enumerate(csqs):
                    values = csq.split("|")
                    for index, field in indexes:
                        values[index] = render(field, values[index])
                    csqs[j] = "|".join(values)
                info[i] = "CSQ=" + ",".join(csqs)
            line[7] = ";".join(info)
            formatted = ("\t".join(line) + "\n").encode()
            changed += formatted != raw
            out.write(formatted)
            digest.update(formatted)
    return {"records": records, "body_md5": digest.hexdigest()}, changed
