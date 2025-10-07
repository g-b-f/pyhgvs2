import re

class HGVSRegex:
    """
    All regular expression for HGVS names.
    """

    @staticmethod
    def _compile_regex(regex) -> re.Pattern[str]:
        return re.compile("^" + regex + "$")

    @staticmethod
    def _compile_regexes(regexes) -> list[re.Pattern[str]]:
        return [re.compile("^" + regex + "$") for regex in regexes]

    # DNA syntax
    # http://www.hgvs.org/mutnomen/standards.html#nucleotide
    BASE = r"[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]|\d+"
    BASES = r"[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+"
    DNA_REF = "(?P<ref>" + BASES + ")"
    DNA_ALT = "(?P<alt>" + BASES + ")"

    # Mutation types
    EQUAL = "(?P<mutation_type>=)"
    SUB = "(?P<mutation_type>>)"
    INS = "(?P<mutation_type>ins)"
    DEL = "(?P<mutation_type>del)"
    DUP = "(?P<mutation_type>dup)"
    INV = "(?P<mutation_type>inv)"

    # Simple coordinate syntax
    COORD_START = r"(?P<start>\d+)"
    COORD_END = r"(?P<end>\d+)"
    COORD_RANGE = COORD_START + "_" + COORD_END

    # cDNA coordinate syntax
    CDNA_COORD = (
        r"(?P<coord_prefix>|-|\*)(?P<coord>\d+)"
        r"((?P<offset_prefix>-|\+)(?P<offset>\d+))?"
    )
    CDNA_START = (
        r"(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)"
        r"((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)"
    )
    CDNA_END = (
        r"(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)"
        r"((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)"
    )
    CDNA_RANGE = CDNA_START + "_" + CDNA_END

    # cDNA allele syntax
    CDNA_ALLELE = [
        # No change
        CDNA_START + DNA_REF + EQUAL,
        # Substitution
        CDNA_START + DNA_REF + SUB + DNA_ALT,
        # 1bp insertion, deletion, duplication
        CDNA_START + INS + DNA_ALT,
        CDNA_START + DEL + DNA_REF,
        CDNA_START + DUP + DNA_REF,
        CDNA_START + DEL,
        CDNA_START + DUP,
        # Insertion, deletion, duplication, inversion
        CDNA_RANGE + INS + DNA_ALT,
        CDNA_RANGE + DEL + DNA_REF,
        CDNA_RANGE + DUP + DNA_REF,
        CDNA_RANGE + DEL,
        CDNA_RANGE + DUP,
        CDNA_RANGE + INV,
        # Indels
        "(?P<delins>" + CDNA_START + "del" + DNA_REF + "ins" + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + "del" + DNA_REF + "ins" + DNA_ALT + ")",
        "(?P<delins>" + CDNA_START + "delins" + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + "delins" + DNA_ALT + ")",
    ]

    CDNA_ALLELE_REGEXES = _compile_regexes(CDNA_ALLELE)
    # Peptide syntax
    PEP = "([A-Z]([a-z]{2}))+"
    PEP_REF = "(?P<ref>" + PEP + ")"
    PEP_REF2 = "(?P<ref2>" + PEP + ")"
    PEP_ALT = "(?P<alt>" + PEP + ")"

    PEP_EXTRA = r"(?P<extra>(|=|\?)(|fs))"

    # Peptide allele syntax
    # fmt: off
    PEP_ALLELE = [
        # No peptide change
        # Example: Glu1161=
        PEP_REF + COORD_START + PEP_EXTRA,
        # Peptide change
        # Example: Glu1161Ser
        PEP_REF + COORD_START + PEP_ALT + PEP_EXTRA,
        # Peptide indel
        # Example: Glu1161_Ser1164?fs
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END + PEP_EXTRA + ")",
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END + PEP_ALT + PEP_EXTRA + ")",
    ]
    # fmt: on
    PEP_ALLELE_REGEXES = _compile_regexes(PEP_ALLELE)

    # Genomic allele syntax
    GENOMIC_ALLELE = [
        # No change
        COORD_START + DNA_REF + EQUAL,
        # Substitution
        COORD_START + DNA_REF + SUB + DNA_ALT,
        # 1bp insertion, deletion, duplication
        COORD_START + INS + DNA_ALT,
        COORD_START + DEL + DNA_REF,
        COORD_START + DUP + DNA_REF,
        COORD_START + DEL,
        COORD_START + DUP,
        # Insertion, deletion, duplication. inversion
        COORD_RANGE + INS + DNA_ALT,
        COORD_RANGE + DEL + DNA_REF,
        COORD_RANGE + DUP + DNA_REF,
        COORD_RANGE + DEL,
        COORD_RANGE + DUP,
        COORD_RANGE + INV,
        # Indels
        "(?P<delins>" + COORD_START + "del" + DNA_REF + "ins" + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + "del" + DNA_REF + "ins" + DNA_ALT + ")",
        "(?P<delins>" + COORD_START + "delins" + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + "delins" + DNA_ALT + ")",
    ]

    GENOMIC_ALLELE_REGEXES = _compile_regexes(GENOMIC_ALLELE)

    # Repeated sequence syntax
    REPEAT_UNIT = r"(?P<repeat_unit>" + BASES + ")"
    REPEAT_COUNT = r"\[(?P<repeat_count>\d+)\]"

    # Pattern for mixed repeats: sequence1[count1]sequence2[count2]...
    # Need to use non-capturing group to avoid the alternation issue with BASES
    MIXED_REPEAT = (
        r"(?P<mixed_repeat>(?:(?:[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)\[\d+\])+)"
    )

    # Repeated sequence patterns for cDNA
    CDNA_REPEAT_ALLELE = [
        # Simple repeat at single position: c.123CAG[26]
        CDNA_START + REPEAT_UNIT + REPEAT_COUNT,
        # Range repeat: c.123_456CAG[26]
        CDNA_RANGE + REPEAT_UNIT + REPEAT_COUNT,
        # Mixed repeat at single position: c.123CAG[19]CAA[4]
        CDNA_START + MIXED_REPEAT,
        # Mixed repeat at range: c.123_456CAG[19]CAA[4]
        CDNA_RANGE + MIXED_REPEAT,
    ]

    CDNA_REPEAT_ALLELE_REGEXES = _compile_regexes(CDNA_REPEAT_ALLELE)

    # Repeated sequence patterns for genomic
    GENOMIC_REPEAT_ALLELE = [
        # Simple repeat at single position: g.123CAG[26]
        COORD_START + REPEAT_UNIT + REPEAT_COUNT,
        # Range repeat: g.123_456CAG[26]
        COORD_RANGE + REPEAT_UNIT + REPEAT_COUNT,
        # Mixed repeat at single position: g.123CAG[19]CAA[4]
        COORD_START + MIXED_REPEAT,
        # Mixed repeat at range: g.123_456CAG[19]CAA[4]
        COORD_RANGE + MIXED_REPEAT,
    ]

    GENOMIC_REPEAT_ALLELE_REGEXES = _compile_regexes(GENOMIC_REPEAT_ALLELE)

    # Prefix parsing patterns
    TRANSCRIPT_GENE_PARENS = r"^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$"
    GENE_TRANSCRIPT_BRACES = r"^(?P<gene>[^{]+)\{(?P<transcript>[^}]+)\}$"
    MIXED_REPEAT_PARSER = r"([ACGTBDHKMNRSVWY]+|\d+)\[(\d+)\]"

    TRANSCRIPT_GENE_PARENS_REGEX = _compile_regex(TRANSCRIPT_GENE_PARENS)
    GENE_TRANSCRIPT_BRACES_REGEX = _compile_regex(GENE_TRANSCRIPT_BRACES)
