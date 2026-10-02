"""FASTA reference parsing and validation."""

from __future__ import annotations

from dataclasses import dataclass


class ReferenceError(ValueError):
    """Raised when the uploaded FASTA reference is invalid."""


@dataclass
class Reference:
    """Ordered collection of reference sequences."""

    sequences: dict[str, str]
    order: list[str]
    rank: dict[str, int]

    @classmethod
    def parse(cls, text: str) -> "Reference":
        sequences: dict[str, str] = {}
        order: list[str] = []
        name: str | None = None
        chunks: list[str] = []

        def flush() -> None:
            if name is None:
                return
            seq = "".join(chunks).upper()
            if not seq:
                raise ReferenceError(f"sequence '{name}' is empty")
            bad = sorted(set(seq) - set("ACGT"))
            if bad:
                raise ReferenceError(
                    f"sequence '{name}' contains non-ACGT characters: {''.join(bad)}"
                )
            sequences[name] = seq
            order.append(name)

        for lineno, raw in enumerate(text.splitlines(), start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                title = line[1:].strip()
                if not title:
                    raise ReferenceError(f"line {lineno}: empty FASTA header")
                name = title.split()[0]
                if name in sequences:
                    raise ReferenceError(f"duplicate sequence name '{name}'")
                chunks = []
            else:
                if name is None:
                    raise ReferenceError(
                        f"line {lineno}: sequence data before any FASTA header"
                    )
                chunks.append(line)
        flush()
        if not order:
            raise ReferenceError("no sequences found in FASTA reference")
        rank = {n: i for i, n in enumerate(order)}
        return cls(sequences=sequences, order=order, rank=rank)
