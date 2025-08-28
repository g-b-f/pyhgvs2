import re
from typing import Optional

from .constants import CHROM_PREFIX
from .exceptions import (
    InvalidKindError,
    NonMatchingAlleleError,
    StartGreaterThanEndError,
)
from .lookups import cdna_to_genomic_coord, get_refseq_type
from .models import CDNACoord, Transcript
from .regex import HGVSRegex
from .variants import revcomp


class HGVSName:
    """
    Represents a HGVS variant name.
    """

    _regexes = HGVSRegex

    def __init__(
        self,
        name="",
        prefix="",
        chrom="",
        transcript="",
        gene="",
        kind="",
        mutation_type=None,
        start=0,
        end=0,
        ref_allele="",
        ref2_allele="",
        alt_allele="",
        cdna_start=None,
        cdna_end=None,
        pep_extra="",
    ):

        # Full HGVS name.
        self.name = name

        # Name parts.
        self.prefix = prefix
        self.chrom = chrom
        self.transcript = transcript
        self.gene = gene
        self.kind = kind
        self.mutation_type = mutation_type
        self.start = start
        self.end = end
        self.ref_allele = ref_allele  # reference allele
        self.ref2_allele = ref2_allele  # reference allele at end of pep indel
        self.alt_allele = alt_allele  # alternate allele

        # cDNA-specific fields
        self.cdna_start = cdna_start if cdna_start else CDNACoord()
        self.cdna_end = cdna_end if cdna_end else CDNACoord()

        # Protein-specific fields
        self.pep_extra = pep_extra

        if name:
            self.parse(name)

    def parse(self, name: str):
        """Parse a HGVS name."""
        # Does HGVS name have transcript/gene prefix?
        if ":" in name:
            prefix, allele = name.split(":", 1)
        else:
            prefix = ""
            allele = name

        self.name = name

        # Parse prefix and allele.
        self.parse_allele(allele)
        self.parse_prefix(prefix, self.kind)
        self._validate()

    def parse_prefix(self, prefix: str, kind):
        """
        Parse a HGVS prefix (gene/transcript/chromosome).

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
          BRCA1{NM_007294.3}:c.2207A>C
        """

        self.prefix = prefix

        # No prefix.
        if prefix == "":
            self.chrom = ""
            self.transcript = ""
            self.gene = ""
            return

        # Transcript and gene given with parens.
        # example: NM_007294.3(BRCA1):c.2207A>C
        match = re.match(r"^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$", prefix)
        if match:
            self.transcript = match.group("transcript")
            self.gene = match.group("gene")
            return

        # Transcript and gene given with braces.
        # example: BRCA1{NM_007294.3}:c.2207A>C
        match = re.match(r"^(?P<gene>[^{]+)\{(?P<transcript>[^}]+)\}$", prefix)
        if match:
            self.transcript = match.group("transcript")
            self.gene = match.group("gene")
            return

        # Determine using Ensembl type.
        if prefix.startswith("ENST"):
            self.transcript = prefix
            return

        # Determine using refseq type.
        refseq_type = get_refseq_type(prefix)
        if refseq_type in ("mRNA", "RNA"):
            self.transcript = prefix
            return

        # Determine using refseq type.
        if prefix.startswith(CHROM_PREFIX) or refseq_type == "genomic":
            self.chrom = prefix
            return

        # Assume gene name.
        self.gene = prefix

    def parse_allele(self, allele: str):
        """
        Parse a HGVS allele description.

        Some examples include:
          cDNA substitution: c.101A>C,
          cDNA indel: c.3428delCinsTA, c.1000_1003delATG, c.1000_1001insATG
          No protein change: p.Glu1161=
          Protein change: p.Glu1161Ser
          Protein frameshift: p.Glu1161_Ser1164?fs
          Genomic substitution: g.1000100A>T
          Genomic indel: g.1000100_1000102delATG
        """
        if "." not in allele:
            raise InvalidKindError(allele)

        # Determine HGVS name kind.
        kind, details = allele.split(".", 1)
        self.kind = kind
        self.mutation_type = None

        if kind == "c":
            self.parse_cdna(details)
        elif kind == "p":
            self.parse_protein(details)
        elif kind == "g":
            self.parse_genome(details)
        else:
            raise NotImplementedError(f"unknown kind: {allele}")

    def parse_cdna(self, details: str):
        """
        Parse a HGVS cDNA name.

        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
        """
        for regex in self._regexes.CDNA_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get("delins"):
                    self.mutation_type = "delins"
                else:
                    self.mutation_type = groups["mutation_type"]

                # Parse coordinates.
                self.cdna_start = CDNACoord(string=groups.get("start"))
                if groups.get("end"):
                    self.cdna_end = CDNACoord(string=groups.get("end"))
                else:
                    self.cdna_end = CDNACoord(string=groups.get("start"))

                # Parse alleles.
                self.ref_allele = groups.get("ref", "")
                self.alt_allele = groups.get("alt", "")

                # Convert numerical allelles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
                return

        raise NonMatchingAlleleError(details, "cDNA")

    def parse_protein(self, details: str):
        """
        Parse a HGVS protein name.

        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs
        """
        for regex in self._regexes.PEP_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups: dict[str, str] = match.groupdict()

                # Parse mutation type.
                if groups.get("delins"):
                    self.mutation_type = "delins"
                else:
                    self.mutation_type = ">"

                # Parse coordinates.
                self.start = int(groups["start"])
                self.end = int(groups.get("end", self.start))

                # Parse alleles.
                self.ref_allele = groups.get("ref", "")
                if groups.get("ref2"):
                    self.ref2_allele = groups.get("ref2")
                    self.alt_allele = groups.get("alt", "")
                else:
                    # If alt is not given, assume matching with ref
                    self.ref2_allele = self.ref_allele
                    self.alt_allele = groups.get("alt", self.ref_allele)

                self.pep_extra = groups.get("extra")
                return

        raise NonMatchingAlleleError(details, "protein")

    def parse_genome(self, details: str):
        """
        Parse a HGVS genomic name.

        Som examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG
        """

        for regex in self._regexes.GENOMIC_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get("delins"):
                    self.mutation_type = "delins"
                else:
                    self.mutation_type = groups["mutation_type"]

                # Parse coordinates.
                self.start = int(groups["start"])
                if groups.get("end"):
                    self.end = int(groups.get("end", self.start))
                else:
                    self.end = self.start

                # Parse alleles.
                self.ref_allele = groups.get("ref", "")
                self.alt_allele = groups.get("alt", "")

                # Convert numerical allelles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
                return

        raise NonMatchingAlleleError(details, "genomic")

    def _validate(self):
        """
        Check for internal inconsistencies in representation
        """
        if self.start > self.end:
            raise StartGreaterThanEndError(self.name)

    def __repr__(self):
        try:
            return f"HGVSName('{self.format()}')"
        except NotImplementedError:
            return f"HGVSName('{self.name}')"

    def __unicode__(self):
        return self.format()

    def format(self, use_prefix=True, use_gene=True) -> str:
        """Generate a HGVS name as a string."""

        if self.kind == "c":
            allele = "c." + self.format_cdna()
        elif self.kind == "p":
            allele = "p." + self.format_protein()
        elif self.kind == "g":
            allele = "g." + self.format_genome()
        else:
            raise NotImplementedError(f"not implemented: '{self.kind}'")

        prefix = self.format_prefix(use_gene=use_gene) if use_prefix else ""

        if prefix:
            return prefix + ":" + allele
        else:
            return allele

    def format_prefix(self, use_gene=True) -> str:
        """
        Generate HGVS trancript/gene prefix.

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
        """

        if self.kind == "g":
            if self.chrom:
                return self.chrom

        if self.transcript:
            if use_gene and self.gene:
                return f"{self.transcript}({self.gene})"
            else:
                return self.transcript
        else:
            if use_gene:
                return self.gene
            else:
                return ""

    def format_cdna_coords(self) -> str:
        """
        Generate HGVS cDNA coordinates string.
        """
        # Format coordinates.
        if self.cdna_start == self.cdna_end:
            return str(self.cdna_start)
        else:
            return f"{self.cdna_start}_{self.cdna_end}"

    def format_dna_allele(self) -> str:
        """
        Generate HGVS DNA allele.
        """
        if self.mutation_type == "=":
            # No change.
            # example: 101A=
            return self.ref_allele + "="

        if self.mutation_type == ">":
            # SNP.
            # example: 101A>C
            return self.ref_allele + ">" + self.alt_allele

        elif self.mutation_type == "delins":
            # Indel.
            # example: 112_117delAGGTCAinsTG, 112_117delinsTG
            return "del" + self.ref_allele + "ins" + self.alt_allele

        elif self.mutation_type in ("del", "dup"):
            # Delete, duplication.
            # example: 1000_1003delATG, 1000_1003dupATG
            return self.mutation_type + self.ref_allele

        elif self.mutation_type == "ins":
            # Insert.
            # example: 1000_1001insATG
            return self.mutation_type + self.alt_allele

        elif self.mutation_type == "inv":
            return self.mutation_type

        else:
            raise ValueError(f"unknown mutation type: '{self.mutation_type}'")

    def format_cdna(self) -> str:
        """
        Generate HGVS cDNA allele.

        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
        """
        return self.format_cdna_coords() + self.format_dna_allele()

    def format_protein(self) -> str:
        """
        Generate HGVS protein name.

        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs
        """
        if (
            self.start == self.end
            and self.ref_allele == self.ref2_allele == self.alt_allele
        ):
            # Match.
            # Example: Glu1161=
            pep_extra = self.pep_extra if self.pep_extra else "="
            return self.ref_allele + str(self.start) + pep_extra

        elif (
            self.start == self.end
            and self.ref_allele == self.ref2_allele
            and self.ref_allele != self.alt_allele
        ):
            # Change.
            # Example: Glu1161Ser
            return self.ref_allele + str(self.start) + self.alt_allele + self.pep_extra

        elif self.start != self.end:
            # Range change.
            # Example: Glu1161_Ser1164?fs
            return (
                self.ref_allele
                + str(self.start)
                + "_"
                + self.ref2_allele
                + str(self.end)
                + self.pep_extra
            )

        else:
            raise NotImplementedError("protein name formatting.")

    def format_coords(self) -> str:
        """
        Generate HGVS cDNA coordinates string.
        """
        if self.start == self.end:
            return str(self.start)
        else:
            return f"{self.start}_{self.end}"

    def format_genome(self) -> str:
        """
        Generate HGVS genomic allele.

        Som examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG
        """
        return self.format_coords() + self.format_dna_allele()

    def get_coords(
        self, transcript: Optional[Transcript] = None
    ) -> tuple[str, int, int]:
        """Return genomic coordinates of reference allele."""
        if self.kind == "c":
            if transcript is None:
                raise ValueError("must pass a transcript")
            chrom = transcript.tx_position.chrom
            start = cdna_to_genomic_coord(transcript, self.cdna_start)
            end = cdna_to_genomic_coord(transcript, self.cdna_end)

            if not transcript.tx_position.is_forward_strand:
                if end > start:
                    raise AssertionError("cdna_start cannot be greater than cdna_end")
                start, end = end, start
            else:
                if start > end:
                    raise AssertionError("cdna_start cannot be greater than cdna_end")

            if self.mutation_type == "ins":
                # Inserts have empty interval.
                if start < end:
                    start += 1
                    end -= 1
                else:
                    end = start - 1
            elif self.mutation_type == "dup":
                end = start - 1

        elif self.kind == "g":
            chrom = self.chrom
            start = self.start
            end = self.end

        else:
            raise NotImplementedError(
                f'Coordinates are not available for this kind of HGVS name "{self.kind}"'
            )

        return chrom, start, end

    def get_vcf_coords(self, transcript=None) -> tuple[str, int, int]:
        """Return genomic coordinates of reference allele in VCF-style."""
        chrom, start, end = self.get_coords(transcript)

        # Inserts and deletes require left-padding by 1 base
        if self.mutation_type in ("=", ">"):
            pass
        elif self.mutation_type in ("del", "ins", "dup", "delins"):
            # Indels have left-padding.
            start -= 1
        else:
            raise NotImplementedError(f"Unknown mutation_type '{self.mutation_type}'")
        return chrom, start, end

    def get_ref_alt(self, is_forward_strand=True) -> tuple[str, str]:
        """Return reference and alternate alleles."""
        if self.kind == "p":
            raise NotImplementedError(
                "get_ref_alt is not implemented for protein HGVS names"
            )

        ref, alt = self.ref_allele, self.alt_allele

        # Represent duplications are inserts.
        if self.mutation_type == "dup":
            alleles = ("", alt[: len(alt) // 2])
        else:
            alleles = (ref, alt)

        if is_forward_strand:
            return alleles
        else:
            return tuple(map(revcomp, alleles))  # type: ignore[return-value]
