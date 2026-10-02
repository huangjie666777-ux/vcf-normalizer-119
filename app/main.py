"""HTTP delivery layer: upload limits, per-request temp files, all-or-nothing
VCF download.
"""

from __future__ import annotations

import os
import tempfile

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import JSONResponse, Response

from .errors import InputError
from .service import normalize_vcf, verify_with_faidx

MAX_UPLOAD_BYTES = int(os.environ.get("VCFNORM_MAX_UPLOAD_BYTES", 16 * 1024 * 1024))
MAX_RECORDS = int(os.environ.get("VCFNORM_MAX_RECORDS", 100_000))

app = FastAPI(title="vcf-normalizer", version="1.0.0")


async def _read_limited(upload: UploadFile, label: str) -> bytes:
    data = bytearray()
    while chunk := await upload.read(1024 * 1024):
        data.extend(chunk)
        if len(data) > MAX_UPLOAD_BYTES:
            raise HTTPException(
                status_code=413,
                detail={
                    "errors": [
                        {
                            "source": label,
                            "line": 0,
                            "reason": f"upload exceeds the "
                            f"{MAX_UPLOAD_BYTES}-byte limit",
                        }
                    ]
                },
            )
    return bytes(data)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/normalize")
async def normalize(fasta: UploadFile = File(...), vcf: UploadFile = File(...)):
    fasta_bytes = await _read_limited(fasta, "fasta")
    vcf_bytes = await _read_limited(vcf, "vcf")
    try:
        fasta_text = fasta_bytes.decode("utf-8")
        vcf_text = vcf_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(
            status_code=422,
            detail={"errors": [{"source": "upload", "line": 0,
                                "reason": f"files must be UTF-8 text: {exc}"}]},
        )

    # Per-request temp directory: concurrent requests never share state and
    # everything is removed on success or failure.
    with tempfile.TemporaryDirectory(prefix="vcfnorm_") as tmpdir:
        try:
            result = normalize_vcf(fasta_text, vcf_text, MAX_RECORDS)
            fasta_path = os.path.join(tmpdir, "reference.fa")
            with open(fasta_path, "w", encoding="utf-8") as fh:
                fh.write(fasta_text)
            verify_with_faidx(fasta_path, result)
        except InputError as exc:
            return JSONResponse(
                status_code=422, content={"errors": [exc.as_dict()]}
            )
    return Response(
        content=result.vcf_text,
        media_type="text/plain",
        headers={
            "Content-Disposition": 'attachment; filename="normalized.vcf"',
            "X-Records": str(result.n_records),
            "X-Records-Changed": str(result.n_changed),
        },
    )
