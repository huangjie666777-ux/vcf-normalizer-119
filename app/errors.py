"""Shared error types carrying the original input line number."""

from __future__ import annotations


class InputError(ValueError):
    """Validation error tied to a line of an uploaded input file."""

    def __init__(self, line_no: int, reason: str, source: str = "vcf") -> None:
        self.line_no = line_no
        self.reason = reason
        self.source = source
        where = f"{source} line {line_no}" if line_no > 0 else source
        super().__init__(f"{where}: {reason}")

    def as_dict(self) -> dict:
        return {"source": self.source, "line": self.line_no, "reason": self.reason}
