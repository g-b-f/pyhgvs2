r"""
Methods for manipulating HGVS names

Recommendations for the HGVS naming standard:
http://www.hgvs.org/mutnomen/standards.html

Definition of which transcript to use coding variants:
ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene


HGVS language currently implemented.

HGVS = ALLELE
     | PREFIX_NAME : ALLELE

PREFIX_NAME = TRANSCRIPT
            | TRANSCRIPT '(' GENE ')'

TRANSCRIPT = TRANSCRIPT_NAME
           | TRANSCRIPT_NAME '.' TRANSCRIPT_VERSION

TRANSCRIPT_VERSION = NUMBER

ALLELE = 'c.' CDNA_ALLELE    # cDNA
       | 'g.' GENOMIC_ALLELE # genomic
       | 'm.' MIT_ALLELE     # mitochondrial sequence
       | 'n.' NC_ALLELE      # non-coding RNA reference sequence
       | 'r.' RNA_ALLELE     # RNA sequence (like r.76a>u)
       | 'p.' PROTEIN_ALLELE # protein sequence (like  p.Lys76Asn)

NC_ALLELE =
RNA_ALLELE =
CDNA_ALLELE = CDNA_COORD SINGLE_BASE_CHANGE
            | CDNA_COORD_RANGE MULTI_BASE_CHANGE

GENOMIC_ALLELE =
MIT_ALLELE = COORD SINGLE_BASE_CHANGE
           | COORD_RANGE MULTI_BASE_CHANGE

SINGLE_BASE_CHANGE = CDNA_ALLELE = CDNA_COORD BASE '='        # no change
                   | CDNA_COORD BASE '>' BASE                 # substitution
                   | CDNA_COORD 'ins' BASE                    # 1bp insertion
                   | CDNA_COORD 'del' BASE                    # 1bp deletion
                   | CDNA_COORD 'dup' BASE                    # 1bp duplication
                   | CDNA_COORD 'ins'                         # 1bp insertion
                   | CDNA_COORD 'del'                         # 1bp deletion
                   | CDNA_COORD 'dup'                         # 1bp duplication
                   | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
                   | CDNA_COORD 'delins' BASE                 # 1bp indel

MULTI_BASE_CHANGE = COORD_RANGE 'del' BASES             # deletion
                  | COORD_RANGE 'ins' BASES             # insertion
                  | COORD_RANGE 'dup' BASES             # duplication
                  | COORD_RANGE 'del'                   # deletion
                  | COORD_RANGE 'dup'                   # duplication
                  | COORD_RANGE 'del' BASES 'ins' BASES # indel
                  | COORD_RANGE 'delins' BASES          # indel


AMINO1 = [GAVLIMFWPSTCYNQDEKRH]

AMINO3 = 'Gly' | 'Ala' | 'Val' | 'Leu' | 'Ile' | 'Met' | 'Phe' | 'Trp' | 'Pro'
       | 'Ser' | 'Thr' | 'Cys' | 'Tyr' | 'Asn' | 'Gln' | 'Asp' | 'Glu' | 'Lys'
       | 'Arg' | 'His'

PROTEIN_ALLELE = AMINO3 COORD '='               # no peptide change
               | AMINO1 COORD '='               # no peptide change
               | AMINO3 COORD AMINO3 PEP_EXTRA  # peptide change
               | AMINO1 COORD AMINO1 PEP_EXTRA  # peptide change
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA        # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA        # indel
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA AMINO3 # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA AMINO1 # indel

# A genomic range:
COORD_RANGE = COORD '_' COORD

# A cDNA range:
CDNA_COORD_RANGE = CDNA_COORD '_' CDNA_COORD

# A cDNA coordinate:
CDNA_COORD = COORD_PREFIX COORD
           | COORD_PREFIX COORD OFFSET_PREFIX OFFSET
COORD_PREFIX = '' | '-' | '*'
COORD = NUMBER
OFFSET_PREFIX = '-' | '+'
OFFSET = NUMBER

# Primitives:
NUMBER = \d+
BASE = [ACGT]
BASES = BASE+

"""

from .functions import (
    format_hgvs_name,
    hgvs_justify_dup,
    hgvs_justify_indel,
    hgvs_normalize_variant,
    parse_hgvs_name,
    variant_to_hgvs_name,
)
from .hgvsg_name import HGVSName
from .models import CDNACoord

__all__ = [
    "HGVSName",
    "hgvs_justify_dup",
    "hgvs_justify_indel",
    "hgvs_normalize_variant",
    "parse_hgvs_name",
    "variant_to_hgvs_name",
    "format_hgvs_name",
    "CDNACoord",
]
