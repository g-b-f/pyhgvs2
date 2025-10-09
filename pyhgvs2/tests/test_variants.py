from typing import Literal
import pytest
from ..variants import normalize_variant
from .genome import MockGenomeTestFile

_genome_seq = {
    ("chr1", 0, 41): "N" * 41,
    ("chr1", 1, 31): "N" * 30,
    ("chr17", 41246230, 41246271): "AGCCTCATGAGGATCACTGGCCAGTAAGTCTATTTTCTCTG",
    ("chr17", 41246218, 41246248): "TTTACATATTAAAGCCTCATGAGGATCACT",
    ("chr17", 41246249, 41246279): "GCCAGTAAGTCTATTTTCTCTGAAGAACCA",
}


@pytest.fixture(scope="module")
def genome() -> MockGenomeTestFile:
    return MockGenomeTestFile(
        db_filename="hg19.fa",
        filename="pyhgvs2/tests/data/test_variants.2.genome",
        create_data=False,
    )


ChromType = tuple[str, int, str, list[str]]

# fmt: off
@pytest.mark.parametrize("variant,true_variant,justify", [
    # Simple SNP.
    (("chr17", 41246250, "G", ["C"]), ("chr17", 41246250, "G", ["C"]), "left"),
    # Left-align and 1bp pad.
    (("chr17", 41246251, "", ["G"]), ("chr17", 41246248, "T", ["TG"]), "left"),
    # Trim common prefix, left-align, and 1bp pad.
    (("chr17", 41246250, "G", ["GG"]), ("chr17", 41246248, "T", ["TG"]), "left"),
    # Trim common prefix.
    (("chr17", 41246248, "TGGC", ["TGGA"]), ("chr17", 41246251, "C", ["A"]), "left"),
    # Trim common prefix and suffix.
    (("chr17", 41246248, "TGGC", ["TGAC"]), ("chr17", 41246250, "G", ["A"]), "left"),
    # Trim common prefix, triallelic
    (
        ("chr17", 41246248, "TGGC", ["TGGA", "TGAC"]),
        ("chr17", 41246249, "GGC", ["GAC", "GGA"]),
        "left",
    ),
    # Left edge of chromosome left justify, right pad.
    (("chr1", 5, "NN", ["N"]), ("chr1", 1, "NN", ["N"]), "left"),
    # Insertion. Trim common prefix, right-align, and 1bp pad.
    (
        ("chr17", 78078933, "T", ["TGGGCA"]),
        ("chr17", 78078946, "G", ["GCAGGG"]),
        "right",
    ),
    # Deletion. Trim common prefix, right-align, and 1bp pad.
    (("chr7", 117199644, "ATCT", ["A"]), ("chr7", 117199645, "TCTT", ["T"]), "right")
    ]) # fmt: on
def test_normalize_variant(genome: MockGenomeTestFile, variant:ChromType, true_variant:ChromType, justify:Literal["left", "right"]):
    """
    Test normalize_variant against known cases.
    """
    chrom, offset, ref, alts = variant
    norm = normalize_variant(chrom, offset, ref, alts, genome, justify=justify)
    assert norm.variant == true_variant, (
        f"Variant failed to normalize {repr(variant)}: "
        f"{repr(norm.variant)} != {repr(true_variant)}"
    )

class TestPosition:
    """
    Test that final position is 1-index and end-inclusive.
    """
    def test_SNP(self, genome: MockGenomeTestFile):
        normed_allele = normalize_variant("chr11", 17417434, "A", ["T"], genome)
        assert normed_allele.position.chrom_start == 17417434
        assert normed_allele.position.chrom_stop == 17417434

    def test_indel_left_adjustment(self, genome: MockGenomeTestFile):
        normed_allele = normalize_variant("chr17", 3552198, "T", ["AT"], genome)
        assert normed_allele.position.chrom_start == 3552192
        assert normed_allele.position.chrom_stop == 3552192

    def test_indel_right_padding(self, genome: MockGenomeTestFile):
        normed_allele = normalize_variant("chr1", 5, "NN", ["N"], genome)
        assert normed_allele.position.chrom_start == 1
        assert normed_allele.position.chrom_stop == 2
