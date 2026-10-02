import pytest
from fastapi.testclient import TestClient

from app.main import app

FASTA = ">chr1\nGGGTTAAACCC\nTATATATATAT\nACCGGTTAACG\n>chr2\nacgtacgtac\ngtacgtacgt\n"

HEADER = (
    "##fileformat=VCFv4.3\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
)

client = TestClient(app)


def post(vcf_text, fasta_text=FASTA):
    return client.post(
        "/normalize",
        files={
            "fasta": ("ref.fa", fasta_text, "text/plain"),
            "vcf": ("in.vcf", vcf_text, "text/plain"),
        },
    )


def test_success_normalizes_sorts_and_updates_end():
    vcf = HEADER + (
        "chr1\t18\td\tTAT\tT\t50\tPASS\tEND=20\tGT\t0/1\t1|1\n"
        "chr1\t16\ti\tT\tTATAT\t60\tPASS\t.\tGT\t0|1\t./.\n"
        "chr1\t4\ts\tT\tC\t99\tPASS\tEND=4\tGT\t1/1\t0/0\n"
        "chr2\t11\tj\tG\tGTAC\t40\tPASS\t.\tGT\t0/1\t0/1\n"
    )
    resp = post(vcf)
    assert resp.status_code == 200
    lines = resp.text.strip().split("\n")
    assert lines[0].startswith("##fileformat")
    assert lines[3].startswith("#CHROM")
    body = lines[4:]
    assert body[0].split("\t")[:5] == ["chr1", "4", "s", "T", "C"]
    assert body[1].split("\t")[:5] == ["chr1", "11", "i", "C", "CTATA"]
    assert body[2].split("\t")[:5] == ["chr1", "11", "d", "CTA", "C"]
    assert body[3].split("\t")[:5] == ["chr2", "11", "j", "G", "GTAC"]
    # END updated to normalized REF end
    assert "END=4" in body[0]
    # phasing and missing GT preserved
    assert body[1].split("\t")[9:] == ["0|1", "./."]
    assert body[2].split("\t")[9:] == ["0/1", "1|1"]


def test_duplicate_records_kept_stable():
    vcf = HEADER + (
        "chr1\t4\ta\tT\tC\t.\t.\t.\tGT\t0/1\t0/0\n"
        "chr1\t4\tb\tT\tC\t.\t.\t.\tGT\t0/1\t0/0\n"
    )
    resp = post(vcf)
    assert resp.status_code == 200
    ids = [ln.split("\t")[2] for ln in resp.text.strip().split("\n") if not ln.startswith("#")]
    assert ids == ["a", "b"]


@pytest.mark.parametrize(
    "record, fragment",
    [
        ("chr1\t4\tx\tT\tC,G\t.\t.\t.\tGT\t0/1\t0/0", "multi-ALT"),
        ("chr1\t4\tx\tT\t<DEL>\t.\t.\t.\tGT\t0/1\t0/0", "symbolic"),
        ("chr1\t4\tx\tT\t*\t.\t.\t.\tGT\t0/1\t0/0", "star"),
        ("chr1\t4\tx\tTAA\tTC\t.\t.\t.\tGT\t0/1\t0/0", "complex"),
        ("chr1\t4\tx\tA\tC\t.\t.\t.\tGT\t0/1\t0/0", "does not match"),
        ("chr9\t4\tx\tT\tC\t.\t.\t.\tGT\t0/1\t0/0", "not present"),
        ("chr1\t0\tx\tG\tC\t.\t.\t.\tGT\t0/1\t0/0", "1-based"),
        ("chr1\t34\tx\tA\tC\t.\t.\t.\tGT\t0/1\t0/0", "past the end"),
        ("chr1\t4\tx\tT\tT\t.\t.\t.\tGT\t0/1\t0/0", "must differ"),
        ("chr1\t4\tx\tT\tN\t.\t.\t.\tGT\t0/1\t0/0", "A/C/G/T"),
        ("chr1\t4\tx\tT\tC\t.\t.\t.\tGT\t0/1\tbad", "invalid GT"),
    ],
)
def test_invalid_records_rejected_with_line_number(record, fragment):
    resp = post(HEADER + record + "\n")
    assert resp.status_code == 422
    body = resp.json()
    assert "line 5" in body["error"]
    assert fragment in body["error"]


def test_bad_fasta_rejected():
    resp = post(HEADER, fasta_text=">chr1\nACGN\n")
    assert resp.status_code == 422
    assert "FASTA" in resp.json()["error"]


def test_missing_header_rejected():
    resp = post("chr1\t4\tx\tT\tC\t.\t.\t.\n")
    assert resp.status_code == 422
