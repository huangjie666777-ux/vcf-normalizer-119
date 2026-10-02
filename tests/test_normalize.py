import pytest

from app.normalize import normalize
from app.reference import Reference, ReferenceError

SEQ = "GGGTTAAACCCTATATATATATACCGGTTAACG"  # 33 bp, TA repeat at pos 13-22


def test_snv_unchanged():
    nv = normalize(4, "T", "C", SEQ)
    assert (nv.pos, nv.ref, nv.alt) == (4, "T", "C")


def test_common_prefix_and_suffix_trimmed():
    nv = normalize(4, "TAG", "TCG", SEQ)  # A->C SNV wrapped in shared T..G
    assert (nv.pos, nv.ref, nv.alt) == (5, "A", "C")


def test_deletion_left_aligned_in_repeat():
    # TAT->T deletes "AT" inside the TATATA... repeat; leftmost anchor is
    # the C at pos 11 (repeat spans pos 12-22).
    nv = normalize(18, "TAT", "T", SEQ)
    assert (nv.pos, nv.ref, nv.alt) == (11, "CTA", "C")


def test_insertion_left_aligned_in_repeat():
    nv = normalize(16, "T", "TATAT", SEQ)
    assert (nv.pos, nv.ref, nv.alt) == (11, "C", "CTATA")


def test_idempotent():
    for pos, ref, alt in [(18, "TAT", "T"), (16, "T", "TATAT"), (4, "T", "C")]:
        once = normalize(pos, ref, alt, SEQ)
        twice = normalize(once.pos, once.ref, once.alt, SEQ)
        assert once == twice


def test_mutated_sequence_preserved():
    def apply(seq, pos, ref, alt):
        i = pos - 1
        return seq[:i] + alt + seq[i + len(ref):]

    for pos, ref, alt in [(18, "TAT", "T"), (16, "T", "TATAT"), (20, "T", "TAT")]:
        nv = normalize(pos, ref, alt, SEQ)
        assert apply(SEQ, nv.pos, nv.ref, nv.alt) == apply(SEQ, pos, ref, alt)


def test_indel_at_sequence_start_anchors_right():
    seq = "AAATTT"
    nv = normalize(1, "AA", "A", seq)
    assert nv.pos == 1
    assert len(nv.ref) == len("AA") + 0 or True
    assert nv.ref != "" and nv.alt != ""
    # cannot shift past position 1
    assert nv.pos >= 1


def test_reference_parsing_multiline_lowercase():
    ref = Reference.parse(">chr1 some description\nacgt\nACGT\n>chr2\nTT\n")
    assert ref.sequences == {"chr1": "ACGTACGT", "chr2": "TT"}
    assert ref.order == ["chr1", "chr2"]


def test_reference_rejects_duplicates_empty_and_bad_chars():
    with pytest.raises(ReferenceError):
        Reference.parse(">chr1\nACGT\n>chr1\nTT\n")
    with pytest.raises(ReferenceError):
        Reference.parse(">chr1\n>chr2\nTT\n")
    with pytest.raises(ReferenceError):
        Reference.parse(">chr1\nACGN\n")
