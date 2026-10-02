"""Variant normalization: trimming and reference-guided left alignment.

The transform is idempotent: normalizing an already normalized record
returns it unchanged. The mutated reference sequence is always preserved.
"""

from __future__ import annotations

from .errors import InputError


def check_ref_matches(
    ref_seq: str, pos: int, ref: str, line_no: int, chrom: str
) -> None:
    """Verify 1-based POS/REF against the reference sequence."""
    end = pos + len(ref) - 1
    if end > len(ref_seq):
        raise InputError(
            line_no,
            f"REF end coordinate {end} exceeds length of '{chrom}' "
            f"({len(ref_seq)})",
            "vcf",
        )
    actual = ref_seq[pos - 1 : end]
    if actual != ref:
        raise InputError(
            line_no,
            f"REF '{ref}' does not match reference '{actual}' at "
            f"{chrom}:{pos}",
            "vcf",
        )


def normalize_variant(
    ref_seq: str, pos: int, ref: str, alt: str, line_no: int, chrom: str
) -> tuple[int, str, str]:
    """Return (pos, ref, alt) trimmed and left-aligned against ref_seq.

    pos is 1-based. ref/alt are upper-case ACGT strings, already validated
    as SNV or pure indel.
    """
    check_ref_matches(ref_seq, pos, ref, line_no, chrom)

    # Trim common suffix, then common prefix.
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1

    if ref and alt:
        # SNV: fully normalized once trimmed.
        return pos, ref, alt

    # Pure indel: exactly one side is empty now.
    indel = ref if ref else alt
    is_deletion = bool(ref)

    # Roll the indel to its leftmost equivalent position, never crossing
    # the sequence start.
    while pos > 1 and ref_seq[pos - 2] == indel[-1]:
        indel = ref_seq[pos - 2] + indel[:-1]
        pos -= 1

    if pos > 1:
        # Anchor on the preceding reference base.
        anchor = ref_seq[pos - 2]
        if is_deletion:
            return pos - 1, anchor + indel, anchor
        return pos - 1, anchor, anchor + indel

    # At the sequence start: anchor on the following reference base.
    if is_deletion:
        if len(indel) >= len(ref_seq):
            raise InputError(
                line_no,
                "deletion spans the entire reference sequence and cannot "
                "be anchored",
                "vcf",
            )
        anchor = ref_seq[len(indel)]
        return 1, indel + anchor, anchor
    anchor = ref_seq[0]
    return 1, anchor, indel + anchor
