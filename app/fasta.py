"""FASTA reference parsing and validation.

Rules: multi-line and mixed-case sequence allowed; the sequence name is the
first whitespace-delimited word of the header; empty sequences, duplicate
names and non-ACGT characters are rejected.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import InputError

_BASES = re.compile(r"[ACGTacgt]+")


@dataclass
class FastaRecord:
    name: str
    sequence: str  # upper-cased


def parse_fasta(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    seen: set[str] = set()
    name: str | None = None
    header_line = 0
    chunks: list[str] = []

    def flush() -> None:
        nonlocal name, chunks, header_line
        if name is None:
            return
        seq = "".join(chunks).upper()
        if not seq:
            raise InputError(header_line, f"sequence '{name}' is empty", "fasta")
        records.append(FastaRecord(name=name, sequence=seq))
        name = None
        chunks = []

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            title = line[1:].strip()
            if not title:
                raise InputError(line_no, "empty FASTA header", "fasta")
            candidate = title.split()[0]
            if candidate in seen:
                raise InputError(
                    line_no, f"duplicate sequence name '{candidate}'", "fasta"
                )
            seen.add(candidate)
            name = candidate
            header_line = line_no
        else:
            if name is None:
                raise InputError(
                    line_no, "sequence data before first header", "fasta"
                )
            if not _BASES.fullmatch(line):
                raise InputError(
                    line_no,
                    "invalid characters in sequence (only A/C/G/T allowed)",
                    "fasta",
                )
            chunks.append(line)
    flush()
    if not records:
        raise InputError(0, "no sequences found", "fasta")
    return records
