import re

# from .constants import CHROM_PREFIX, CDNA_STOP_CODON, CDNA_START_CODON


class HGVSRegex:
    """
    All regular expression for HGVS names.
    """

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

    CDNA_ALLELE_REGEXES = [re.compile("^" + regex + "$") for regex in CDNA_ALLELE]

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
    PEP_ALLELE_REGEXES = [re.compile("^" + regex + "$") for regex in PEP_ALLELE]

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

    GENOMIC_ALLELE_REGEXES = [re.compile("^" + regex + "$") for regex in GENOMIC_ALLELE]
