"""VCF 4.3 parsing, validation, normalization and rendering."""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .normalize import normalize
from .reference import Reference

MAX_RECORDS = 100_000

_GT_PATTERN = re.compile(r"^(\.|\d+)([/|](\.|\d+))*$")
_ACGT = re.compile(r"^[ACGT]+$")


class VCFError(ValueError):
    """Raised when the uploaded VCF is invalid; message cites the input line."""


@dataclass
class Record:
    lineno: int
    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str
    qual: str
    filt: str
    info: str
    fmt: str | None
    samples: list[str]
    sort_key: tuple = field(default=None)  # type: ignore[assignment]


@dataclass
class VCFDocument:
    meta_lines: list[str]
    column_header: str
    sample_names: list[str]
    records: list[Record]


def _fail(lineno: int, reason: str) -> None:
    raise VCFError(f"line {lineno}: {reason}")


def _validate_gt(lineno: int, fmt: str, samples: list[str]) -> None:
    keys = fmt.split(":")
    if not keys or keys[0] != "GT":
        _fail(lineno, "FORMAT must start with GT when sample columns are present")
    for sample in samples:
        values = sample.split(":")
        if len(values) > len(keys):
            _fail(lineno, "sample column has more fields than FORMAT declares")
        gt = values[0] if values else ""
        if not _GT_PATTERN.match(gt):
            _fail(lineno, f"invalid GT value '{gt}'")


def parse_vcf(text: str, reference: Reference) -> VCFDocument:
    meta_lines: list[str] = []
    column_header: str | None = None
    sample_names: list[str] = []
    records: list[Record] = []

    lines = text.splitlines()
    for lineno, raw in enumerate(lines, start=1):
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        if line.startswith("##"):
            if column_header is not None:
                _fail(lineno, "meta line found after the #CHROM header")
            meta_lines.append(line)
            continue
        if line.startswith("#"):
            if column_header is not None:
                _fail(lineno, "duplicate #CHROM header line")
            cols = line.split("\t")
            fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            if len(cols) < 8 or cols[:8] != fixed:
                _fail(lineno, "header must start with the 8 fixed VCF columns")
            if len(cols) > 8:
                if len(cols) < 10 or cols[8] != "FORMAT":
                    _fail(lineno, "sample columns require a FORMAT column")
                sample_names = cols[9:]
                if not sample_names:
                    _fail(lineno, "FORMAT column present without any sample")
            column_header = line
            continue
        if column_header is None:
            _fail(lineno, "record found before the #CHROM header line")

        fields = line.split("\t")
        expected = 8 + (1 + len(sample_names) if sample_names else 0)
        if len(fields) != expected:
            _fail(
                lineno,
                f"expected {expected} tab-separated fields, got {len(fields)}",
            )
        chrom, pos_s, vid, ref, alt, qual, filt, info = fields[:8]
        fmt = fields[8] if sample_names else None
        samples = fields[9:] if sample_names else []

        if chrom not in reference.sequences:
            _fail(lineno, f"chromosome '{chrom}' not present in the reference")
        seq = reference.sequences[chrom]
        if not pos_s.isdigit():
            _fail(lineno, f"POS '{pos_s}' is not a positive integer")
        pos = int(pos_s)
        if pos < 1:
            _fail(lineno, "POS must be 1-based (>= 1)")
        if pos + len(ref) - 1 > len(seq):
            _fail(lineno, "REF extends past the end of the reference sequence")
        if not _ACGT.match(ref):
            _fail(lineno, "REF must contain only A/C/G/T")
        if "," in alt:
            _fail(lineno, "multi-ALT records are not supported")
        if alt == "*":
            _fail(lineno, "star (*) alleles are not supported")
        if alt.startswith("<") or alt.endswith(">"):
            _fail(lineno, "symbolic alleles are not supported")
        if "[" in alt or "]" in alt:
            _fail(lineno, "breakend alleles are not supported")
        if not _ACGT.match(alt):
            _fail(lineno, "ALT must contain only A/C/G/T")
        if ref == alt:
            _fail(lineno, "REF and ALT must differ")
        if len(ref) > 1 and len(alt) > 1 and len(ref) != len(alt):
            _fail(lineno, "complex substitutions (MNP/indel mix) are not supported")
        observed = seq[pos - 1 : pos - 1 + len(ref)]
        if observed != ref:
            _fail(
                lineno,
                f"REF '{ref}' does not match reference '{observed}' at {chrom}:{pos}",
            )
        if sample_names:
            _validate_gt(lineno, fmt or "", samples)

        records.append(
            Record(
                lineno=lineno,
                chrom=chrom,
                pos=pos,
                vid=vid,
                ref=ref,
                alt=alt,
                qual=qual,
                filt=filt,
                info=info,
                fmt=fmt,
                samples=samples,
            )
        )

    if column_header is None:
        raise VCFError("missing #CHROM header line")
    if len(records) > MAX_RECORDS:
        raise VCFError(f"too many records: {len(records)} > {MAX_RECORDS}")
    return VCFDocument(meta_lines, column_header, sample_names, records)


def _update_end(info: str, end: int) -> str:
    if info == "." or "END=" not in info:
        return info
    parts = [
        f"END={end}" if p.startswith("END=") else p
        for p in info.split(";")
    ]
    return ";".join(parts)


def normalize_document(doc: VCFDocument, reference: Reference) -> list[str]:
    """Normalize every record and return sorted output lines (records only)."""
    for rec in doc.records:
        seq = reference.sequences[rec.chrom]
        nv = normalize(rec.pos, rec.ref, rec.alt, seq)
        rec.pos, rec.ref, rec.alt = nv.pos, nv.ref, nv.alt
        rec.info = _update_end(rec.info, rec.pos + len(rec.ref) - 1)
        rec.sort_key = (reference.rank[rec.chrom], rec.pos, rec.ref, rec.alt)

    ordered = sorted(
        enumerate(doc.records), key=lambda item: (item[1].sort_key, item[0])
    )
    out: list[str] = []
    for _, rec in ordered:
        fields = [
            rec.chrom,
            str(rec.pos),
            rec.vid,
            rec.ref,
            rec.alt,
            rec.qual,
            rec.filt,
            rec.info,
        ]
        if rec.fmt is not None:
            fields.append(rec.fmt)
            fields.extend(rec.samples)
        out.append("\t".join(fields))
    return out


def render(doc: VCFDocument, record_lines: list[str]) -> str:
    return "\n".join([*doc.meta_lines, doc.column_header, *record_lines]) + "\n"
