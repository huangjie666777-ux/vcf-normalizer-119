"""HTTP layer: upload FASTA+VCF, validate, normalize, download."""

from __future__ import annotations

import os
import tempfile
import uuid

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse
from starlette.background import BackgroundTask

from .reference import Reference, ReferenceError
from .vcfio import MAX_RECORDS, VCFError, normalize_document, parse_vcf, render

MAX_UPLOAD_BYTES = 10 * 1024 * 1024  # 10 MiB per file

app = FastAPI(title="VCF Normalizer", version="1.0.0")


async def _read_limited(upload: UploadFile, label: str) -> str:
    data = await upload.read(MAX_UPLOAD_BYTES + 1)
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"{label} exceeds the {MAX_UPLOAD_BYTES}-byte upload limit",
        )
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(
            status_code=422, detail=f"{label} is not valid UTF-8 text"
        ) from exc


@app.post("/normalize")
async def normalize_endpoint(
    fasta: UploadFile = File(...), vcf: UploadFile = File(...)
) -> FileResponse:
    fasta_text = await _read_limited(fasta, "FASTA reference")
    vcf_text = await _read_limited(vcf, "VCF file")

    try:
        reference = Reference.parse(fasta_text)
    except ReferenceError as exc:
        return JSONResponse(status_code=422, content={"error": f"FASTA: {exc}"})

    try:
        doc = parse_vcf(vcf_text, reference)
        record_lines = normalize_document(doc, reference)
    except VCFError as exc:
        return JSONResponse(status_code=422, content={"error": f"VCF: {exc}"})

    output = render(doc, record_lines)
    fd, tmp_path = tempfile.mkstemp(
        prefix=f"vcfnorm-{uuid.uuid4().hex}-", suffix=".vcf"
    )
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(output)

    def _cleanup() -> None:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    # The background task runs after the response body has been sent, so each
    # concurrent request gets its own temp file that is always removed.
    return FileResponse(
        tmp_path,
        media_type="text/plain",
        filename="normalized.vcf",
        background=BackgroundTask(_cleanup),
    )


@app.get("/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/limits")
async def limits() -> dict[str, int]:
    return {"max_upload_bytes": MAX_UPLOAD_BYTES, "max_records": MAX_RECORDS}
