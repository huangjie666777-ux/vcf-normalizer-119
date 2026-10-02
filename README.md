# VCF Normalizer Backend

FastAPI backend that normalizes VCF 4.3 variants against a FASTA reference so
that equivalent variants share one canonical representation (trimming +
left alignment), for sequencing summary pipelines. No frontend.

## Scope

- Accepts **single-ALT SNVs and pure insertions/deletions** only.
  Multi-ALT, symbolic alleles (`<DEL>`, breakends), star alleles (`*`) and
  complex replacements (e.g. `REF=AT ALT=GC`, `REF=C ALT=TG`) are rejected.
- FASTA: multi-line and mixed-case allowed; the sequence name is the first
  word of the header. Empty sequences, duplicate names and non-ACGT
  characters are rejected.
- POS is 1-based. CHROM must exist in the reference, REF must lie within the
  sequence and match it exactly. Every validation error reports the **original
  input line number** and reason; nothing is skipped.
- Normalization: shared prefix/suffix trimmed, indels left-aligned against
  the reference with the required anchor base kept, never crossing the
  sequence start. The transform is idempotent and preserves the mutated
  sequence.
- Output records are stably sorted by FASTA sequence order, normalized POS,
  REF, ALT. Duplicate records are kept, never merged. Header lines, sample
  order, GT phasing/missing values and other annotations are preserved;
  `INFO/END` is updated to the normalized REF end coordinate.
- **All-or-nothing**: the normalized VCF is downloadable only when every
  record validates; any failure returns `422` with a JSON error body.
- Limits (env-configurable): `VCFNORM_MAX_UPLOAD_BYTES` (default 16 MiB per
  file), `VCFNORM_MAX_RECORDS` (default 100000). Each request works in its
  own temporary directory, cleaned up on success or failure.

## Layout

- `app/fasta.py` — reference reading/validation
- `app/vcf_io.py` — VCF parsing and per-record validation
- `app/normalize.py` — trimming + left alignment
- `app/service.py` — orchestration, sorting, rendering, pysam cross-check
- `app/main.py` — HTTP delivery (upload limits, temp files, download)

## Commands

```bash
.venv/bin/python -m pytest tests -q          # tests
.venv/bin/python -m compileall app           # build check
.venv/bin/uvicorn app.main:app --port 8000   # run
```

## Demo (tandem-repeat example)

`examples/repeat_region.vcf` contains a deletion and an insertion inside a
CA tandem repeat that are not left-aligned, plus an SNV:

```bash
curl -s -F fasta=@examples/repeat_region.fa -F vcf=@examples/repeat_region.vcf \
  http://127.0.0.1:8000/normalize -o normalized.vcf -D -
```

The deletion `chr1:7 ACA>A` normalizes to `chr1:3 GCA>G` (leftmost copy of
the repeat) and `INFO/END` is updated accordingly.

Error responses are JSON, e.g.:

```json
{"errors": [{"source": "vcf", "line": 6, "reason": "multi-ALT records are not supported"}]}
```
