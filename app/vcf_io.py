"""VCF 4.3 parsing and per-record validation.

Only bi-allelic SNVs and pure insertions/deletions are accepted. Multi-ALT,
symbolic alleles, star alleles and complex replacements are rejected. Every
error reports the original input line number; nothing is skipped silently.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .errors import InputError

_BASES = re.compile(r"[ACGTacgt]+")
_GT = re.compile(r"(\.|\d+)([|/](\.|\d+))*")
_FIXED_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


@dataclass
class VcfRecord:
    line_no: int
    fields: list[str]
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass
class VcfDocument:
    meta_lines: list[str] = field(default_factory=list)  # ##... lines
    column_header: str = ""  # #CHROM line, kept verbatim
    samples: list[str] = field(default_factory=list)
    records: list[VcfRecord] = field(default_factory=list)


def _validate_gt(value: str, line_no: int, sample: str) -> None:
    if not _GT.fullmatch(value):
        raise InputError(
            line_no, f"invalid GT value '{value}' in sample {sample}", "vcf"
        )
    for allele in re.split(r"[|/]", value):
        if allele != "." and int(allele) > 1:
            raise InputError(
                line_no,
                f"GT allele index {allele} out of range for a single-ALT record "
                f"in sample {sample}",
                "vcf",
            )


def _validate_genotypes(doc: VcfDocument, rec: VcfRecord) -> None:
    n_fields = len(rec.fields)
    if not doc.samples:
        if n_fields != 8:
            raise InputError(
                rec.line_no,
                f"expected 8 columns (header declares no samples), got {n_fields}",
                "vcf",
            )
        return
    expected = 9 + len(doc.samples)
    if n_fields != expected:
        raise InputError(
            rec.line_no,
            f"expected {expected} columns (8 fixed + FORMAT + "
            f"{len(doc.samples)} samples), got {n_fields}",
            "vcf",
        )
    format_keys = rec.fields[8].split(":")
    if not rec.fields[8] or any(k == "" for k in format_keys):
        raise InputError(rec.line_no, "malformed FORMAT column", "vcf")
    if "GT" not in format_keys:
        return
    gt_idx = format_keys.index("GT")
    for sample_name, sample_col in zip(doc.samples, rec.fields[9:]):
        subfields = sample_col.split(":")
        if gt_idx >= len(subfields):
            raise InputError(
                rec.line_no,
                f"sample {sample_name} has no GT subfield declared by FORMAT",
                "vcf",
            )
        _validate_gt(subfields[gt_idx], rec.line_no, sample_name)


def parse_vcf(text: str, chroms: set, max_records: int) -> VcfDocument:
    doc = VcfDocument()
    header_seen = False
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        if not line:
            continue
        if line.startswith("##"):
            if header_seen:
                raise InputError(
                    line_no, "meta line found after #CHROM header", "vcf"
                )
            doc.meta_lines.append(line)
            continue
        if line.startswith("#"):
            if header_seen:
                raise InputError(line_no, "duplicate #CHROM header line", "vcf")
            header_seen = True
            doc.column_header = line
            cols = line.lstrip("#").split("\t")
            if cols[:8] != _FIXED_COLUMNS:
                raise InputError(
                    line_no,
                    "column header must start with "
                    + "\t".join(_FIXED_COLUMNS),
                    "vcf",
                )
            if len(cols) == 9:
                raise InputError(
                    line_no, "FORMAT column present without sample columns", "vcf"
                )
            if len(cols) > 9:
                if cols[8] != "FORMAT":
                    raise InputError(
                        line_no, "9th column must be FORMAT", "vcf"
                    )
                doc.samples = cols[9:]
            continue
        if not header_seen:
            raise InputError(
                line_no, "data line before #CHROM header", "vcf"
            )
        fields = line.split("\t")
        if len(fields) < 8:
            raise InputError(
                line_no, f"expected at least 8 columns, got {len(fields)}", "vcf"
            )
        chrom, pos_s, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        if chrom not in chroms:
            raise InputError(
                line_no, f"CHROM '{chrom}' not present in the reference", "vcf"
            )
        if not pos_s.isdigit():
            raise InputError(line_no, f"POS '{pos_s}' is not an integer", "vcf")
        pos = int(pos_s)
        if pos < 1:
            raise InputError(line_no, "POS must be >= 1 (1-based)", "vcf")
        ref_u = ref.upper()
        alt_u = alt.upper()
        if not _BASES.fullmatch(ref):
            raise InputError(
                line_no, f"REF '{ref}' must contain only A/C/G/T", "vcf"
            )
        if "," in alt:
            raise InputError(
                line_no, "multi-ALT records are not supported", "vcf"
            )
        if alt == "*":
            raise InputError(line_no, "star allele is not supported", "vcf")
        if alt.startswith("<") or "[" in alt or "]" in alt:
            raise InputError(
                line_no, f"symbolic/breakend ALT '{alt}' is not supported", "vcf"
            )
        if not _BASES.fullmatch(alt):
            raise InputError(
                line_no, f"ALT '{alt}' must contain only A/C/G/T", "vcf"
            )
        if ref_u == alt_u:
            raise InputError(line_no, "REF and ALT must differ", "vcf")
        # Pure indels reduce to an empty side after trimming shared
        # prefix/suffix; anything else is a complex replacement.
        trim_r, trim_a = ref_u, alt_u
        while trim_r and trim_a and trim_r[-1] == trim_a[-1]:
            trim_r, trim_a = trim_r[:-1], trim_a[:-1]
        while trim_r and trim_a and trim_r[0] == trim_a[0]:
            trim_r, trim_a = trim_r[1:], trim_a[1:]
        if trim_r and trim_a and (len(trim_r) > 1 or len(trim_a) > 1):
            raise InputError(
                line_no,
                f"complex replacement REF={ref} ALT={alt} is not supported "
                "(only SNVs and pure indels)",
                "vcf",
            )
        rec = VcfRecord(
            line_no=line_no,
            fields=fields,
            chrom=chrom,
            pos=pos,
            ref=ref_u,
            alt=alt_u,
        )
        _validate_genotypes(doc, rec)
        doc.records.append(rec)
        if len(doc.records) > max_records:
            raise InputError(
                line_no,
                f"record count exceeds the limit of {max_records}",
                "vcf",
            )
    if not header_seen:
        raise InputError(0, "missing #CHROM column header line", "vcf")
    return doc
