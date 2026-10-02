import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

FASTA = b">chr1 demo\nTTGCACACACACAGG\n"
HEADER = (
    "##fileformat=VCFv4.3\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
)


def post(fasta=FASTA, vcf=b""):
    return client.post(
        "/normalize",
        files={"fasta": ("ref.fa", fasta), "vcf": ("in.vcf", vcf)},
    )


def vcf_line(body: str) -> bytes:
    return (HEADER + body).encode()


def test_health():
    assert client.get("/health").json() == {"status": "ok"}


def test_success_downloads_vcf():
    r = post(vcf=vcf_line("chr1\t7\t.\tACA\tA\t.\t.\t.\tGT\t0|1\n"))
    assert r.status_code == 200
    assert "attachment" in r.headers["content-disposition"]
    assert "chr1\t3\t.\tGCA\tG" in r.text
    assert "0|1" in r.text  # phasing preserved


@pytest.mark.parametrize(
    "alt,fragment",
    [
        ("C,A", "multi-ALT"),
        ("*", "star allele"),
        ("<DEL>", "symbolic"),
        ("TG", "complex replacement"),
    ],
)
def test_rejected_alts(alt, fragment):
    r = post(vcf=vcf_line(f"chr1\t4\t.\tC\t{alt}\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    err = r.json()["errors"][0]
    assert err["line"] == 3
    assert fragment in err["reason"]


def test_ref_mismatch_reports_line():
    r = post(vcf=vcf_line("chr1\t4\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    err = r.json()["errors"][0]
    assert err["line"] == 3
    assert "does not match reference" in err["reason"]


def test_bad_gt_reports_line_and_sample():
    r = post(vcf=vcf_line("chr1\t4\t.\tC\tT\t.\t.\t.\tGT\t0/2\n"))
    assert r.status_code == 422
    err = r.json()["errors"][0]
    assert err["line"] == 3
    assert "s1" in err["reason"]


def test_unknown_chrom():
    r = post(vcf=vcf_line("chrX\t4\t.\tC\tT\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    assert "not present in the reference" in r.json()["errors"][0]["reason"]


def test_all_or_nothing():
    vcf = vcf_line(
        "chr1\t4\t.\tC\tT\t.\t.\t.\tGT\t0/1\n"
        "chr1\t5\t.\tA\t*\t.\t.\t.\tGT\t0/1\n"
    )
    r = post(vcf=vcf)
    assert r.status_code == 422
    assert "attachment" not in r.headers.get("content-disposition", "")


def test_fasta_duplicate_name():
    r = post(fasta=b">chr1\nACGT\n>chr1 again\nACGT\n",
             vcf=vcf_line("chr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    err = r.json()["errors"][0]
    assert err["source"] == "fasta"
    assert "duplicate" in err["reason"]


def test_fasta_invalid_base():
    r = post(fasta=b">chr1\nACGN\n",
             vcf=vcf_line("chr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    assert r.json()["errors"][0]["line"] == 2


def test_fasta_empty_sequence():
    r = post(fasta=b">chr1\n>chr2\nACGT\n",
             vcf=vcf_line("chr2\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"))
    assert r.status_code == 422
    assert "empty" in r.json()["errors"][0]["reason"]
