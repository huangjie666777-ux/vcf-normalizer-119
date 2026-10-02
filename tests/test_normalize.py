from app.normalize import normalize_variant
from app.service import normalize_vcf

REF = "TTGCACACACACAGG"  # chr1 in examples/repeat_region.fa


def test_snv_unchanged():
    assert normalize_variant(REF, 4, "C", "T", 1, "chr1") == (4, "C", "T")


def test_deletion_left_aligned_in_tandem_repeat():
    # Delete one CA unit at a non-leftmost copy.
    pos, ref, alt = normalize_variant(REF, 7, "ACA", "A", 1, "chr1")
    assert (pos, ref, alt) == (3, "GCA", "G")


def test_insertion_left_aligned_in_tandem_repeat():
    pos, ref, alt = normalize_variant(REF, 9, "A", "ACA", 1, "chr1")
    assert (pos, ref, alt) == (3, "G", "GCA")


def test_common_prefix_and_suffix_trimmed():
    # REF/ALT share prefix A and suffix G; reduces to SNV C>T.
    assert normalize_variant("AACGT", 2, "ACG", "ATG", 1, "x") == (3, "C", "T")


def test_deletion_at_sequence_start_uses_following_anchor():
    assert normalize_variant("ACGT", 1, "AC", "C", 1, "x") == (1, "AC", "C")


def test_does_not_cross_sequence_start():
    # Homopolymer at the very start: left-align stops at position 1.
    pos, ref, alt = normalize_variant("AAAG", 2, "AA", "A", 1, "x")
    assert (pos, ref, alt) == (1, "AA", "A")


def test_idempotent():
    for args in [(REF, 7, "ACA", "A"), (REF, 9, "A", "ACA"), ("AAAG", 2, "AA", "A")]:
        once = normalize_variant(args[0], *args[1:], 1, "chr1")
        twice = normalize_variant(args[0], *once, 1, "chr1")
        assert once == twice


def test_mutated_sequence_preserved():
    ref = REF
    pos, r, a = normalize_variant(ref, 7, "ACA", "A", 1, "chr1")
    alt_seq = ref[: pos - 1] + a + ref[pos - 1 + len(r) :]
    orig_alt = ref[:6] + "A" + ref[6 + 3 :]
    assert alt_seq == orig_alt


FASTA = ">chr1\n" + REF + "\n>chr2\nGGGAAAAAGCC\n"


def test_service_sorts_stably_and_keeps_duplicates():
    vcf = (
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr2\t4\tb\tA\tG\t.\t.\t.\n"
        "chr1\t7\td1\tACA\tA\t.\t.\t.\n"
        "chr1\t3\td2\tGCA\tG\t.\t.\t.\n"
        "chr1\t4\ts\tC\tT\t.\t.\t.\n"
    )
    out = normalize_vcf(FASTA, vcf, 1000).vcf_text.splitlines()
    body = [l for l in out if not l.startswith("#")]
    keys = [(l.split("\t")[0], int(l.split("\t")[1])) for l in body]
    assert keys == [("chr1", 3), ("chr1", 3), ("chr1", 4), ("chr2", 4)]
    # d1 normalized to the same key as d2; d2 came first in input? No:
    # d1 (line 3) sorts before d2 (line 4)? Stable: d1 appears before d2.
    ids = [l.split("\t")[2] for l in body]
    assert ids == ["d1", "d2", "s", "b"]


def test_service_updates_end():
    vcf = (
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t7\t.\tACA\tA\t.\t.\tEND=9;DP=30\n"
    )
    out = normalize_vcf(FASTA, vcf, 1000).vcf_text
    assert "chr1\t3\t.\tGCA\tG\t.\t.\tEND=5;DP=30" in out
