# VCF Normalizer Backend

FastAPI backend that normalizes VCF 4.3 variants against an uploaded FASTA
reference so that equivalent variants share one canonical representation
(common-affix trimming + reference-guided left alignment of pure indels).

## Scope

- Accepts single-ALT SNVs and pure insertions/deletions (REF/ALT only A/C/G/T).
- Rejects multi-ALT, symbolic alleles (`<DEL>`), star alleles (`*`),
  breakends and complex substitutions, with the original input line number.
- FASTA: multi-line and mixed-case allowed; sequence name is the first word
  of the header; empty sequences, duplicate names and non-ACGT bases are
  rejected.
- POS is 1-based; chromosome existence, coordinate range and REF/reference
  concordance are verified for every record.
- Indels are shifted to the leftmost equivalent position without crossing
  the sequence start; the anchor base is kept; normalization is idempotent
  and preserves the mutated sequence.
- Records are stably sorted by FASTA sequence order, normalized POS, REF,
  ALT; duplicates are kept, never merged. Headers, sample order, GT phasing
  and missing values are preserved; `INFO/END` is updated to the
  normalized REF end coordinate.
- All-or-nothing: any invalid record returns a JSON error (HTTP 422) and no
  file is delivered. Uploads are limited to 10 MiB per file and 100 000
  records; each request writes its own temp file which is removed after the
  response is sent.

## Layout

- `app/reference.py` — FASTA parsing/validation
- `app/vcfio.py` — VCF parsing, per-record validation, sorting, rendering
- `app/normalize.py` — trimming + left alignment
- `app/main.py` — HTTP layer (upload limits, temp files, JSON errors)
- `examples/repeat.fa` / `examples/repeat.vcf` — tandem-repeat demo input

## Commands

```bash
# tests
.venv/bin/python -m pytest tests -q

# run
.venv/bin/python -m uvicorn app.main:app --port 8000

# demo (tandem-repeat example)
curl -s -F fasta=@examples/repeat.fa -F vcf=@examples/repeat.vcf \
  http://127.0.0.1:8000/normalize -o normalized.vcf

# error case returns JSON, e.g. a multi-ALT record:
# {"error":"VCF: line 8: multi-ALT records are not supported"}
```
