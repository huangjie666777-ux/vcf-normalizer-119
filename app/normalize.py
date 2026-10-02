"""Variant normalization: trimming and left-alignment against a reference."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NormalizedVariant:
    pos: int  # 1-based
    ref: str
    alt: str


def normalize(pos: int, ref: str, alt: str, seq: str) -> NormalizedVariant:
    """Normalize a single-ALT SNV or pure indel.

    pos is 1-based; seq is the full reference sequence of the contig.
    Strips removable common prefix/suffix, left-aligns pure indels to the
    leftmost equivalent position (never past the sequence start) and keeps
    the anchor base required by VCF. The mutated sequence is preserved, and
    normalizing the result again is a fixed point.
    """
    start = pos - 1  # 0-based offset of ref[0] in seq
    r = ref
    a = alt

    # Trim common suffix, then common prefix.
    while r and a and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    while r and a and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        start += 1

    if not r or not a:
        # Pure indel: rotate the indel sequence left while the preceding
        # reference base matches the indel's last base. Each rotation keeps
        # the mutated sequence identical.
        indel = r if r else a
        insertion = not r
        while start > 0 and seq[start - 1] == indel[-1]:
            indel = indel[-1] + indel[:-1]
            start -= 1
        if insertion:
            r, a = "", indel
        else:
            r, a = indel, ""
        # Restore the mandatory anchor base.
        if start > 0:
            base = seq[start - 1]
            r = base + r
            a = base + a
            start -= 1
        else:
            # At the very sequence start: anchor on the right instead.
            base = seq[start + len(r)]
            r = r + base
            a = a + base

    return NormalizedVariant(pos=start + 1, ref=r, alt=a)
