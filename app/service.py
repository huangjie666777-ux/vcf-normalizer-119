"""Business orchestration: reference loading, normalization, sorting and
VCF rendering. All-or-nothing: any InputError aborts the whole batch.
"""

from __future__ import annotations

from dataclasses import dataclass

import pysam

from .fasta import parse_fasta
from .normalize import normalize_variant
from .vcf_io import VcfDocument, parse_vcf


@dataclass
class NormResult:
    vcf_text: str
    n_records: int
    n_changed: int


def _update_end(info: str, end: int) -> str:
    """Rewrite INFO/END to the normalized REF end coordinate, if present."""
    if info in (".", ""):
        return info
    parts = info.split(";")
    for i, part in enumerate(parts):
        if part.startswith("END="):
            parts[i] = f"END={end}"
    return ";".join(parts)


def normalize_vcf(fasta_text: str, vcf_text: str, max_records: int) -> NormResult:
    fasta_records = parse_fasta(fasta_text)
    ref_by_name = {r.name: r.sequence for r in fasta_records}
    chrom_order = {r.name: i for i, r in enumerate(fasta_records)}

    doc: VcfDocument = parse_vcf(vcf_text, set(ref_by_name), max_records)

    rendered: list[tuple[tuple, str]] = []
    n_changed = 0
    for rec in doc.records:
        ref_seq = ref_by_name[rec.chrom]
        pos, ref, alt = normalize_variant(
            ref_seq, rec.pos, rec.ref, rec.alt, rec.line_no, rec.chrom
        )
        if (pos, ref, alt) != (rec.pos, rec.ref, rec.alt):
            n_changed += 1
        fields = list(rec.fields)
        fields[1] = str(pos)
        fields[3] = ref
        fields[4] = alt
        fields[7] = _update_end(fields[7], pos + len(ref) - 1)
        key = (chrom_order[rec.chrom], pos, ref, alt)
        rendered.append((key, "\t".join(fields)))

    # Stable sort: equal keys keep the original input order; duplicates are
    # never merged.
    rendered.sort(key=lambda item: item[0])

    out_lines = list(doc.meta_lines)
    out_lines.append(doc.column_header)
    out_lines.extend(line for _, line in rendered)
    return NormResult(
        vcf_text="\n".join(out_lines) + "\n",
        n_records=len(doc.records),
        n_changed=n_changed,
    )


def verify_with_faidx(fasta_path: str, result: NormResult) -> None:
    """Cross-check every normalized REF against the faidx-indexed FASTA via
    pysam random access. Raises AssertionError on any disagreement."""
    pysam.faidx(fasta_path)
    with pysam.FastaFile(fasta_path) as fa:
        for line in result.vcf_text.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            chrom, pos, ref = f[0], int(f[1]), f[3]
            fetched = fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
            if fetched != ref:
                raise AssertionError(
                    f"normalized REF mismatch at {chrom}:{pos}: "
                    f"{ref} != {fetched}"
                )
