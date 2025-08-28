"""
Models for representing genomic elements.
"""

import re
from dataclasses import dataclass
from typing import Mapping, Optional

from .constants import CDNA_START_CODON, CDNA_STOP_CODON

GenomeType = Mapping[str, str]


@dataclass
class Position:
    """A position in the genome."""

    chrom: str
    chrom_start: int
    chrom_stop: int
    is_forward_strand: bool

    def __str__(self) -> str:
        return f"{self.chrom}[{self.chrom_start}:{self.chrom_stop}]"


@dataclass
class Gene:
    name: str


class Transcript:
    """RefGene Transcripts for hg19

    A gene may have multiple transcripts with different combinations of exons.
    """

    def __init__(
        self,
        name: str,
        version: Optional[int],
        gene: str,
        tx_position: Position,
        cds_position: Position,
        is_default=False,
        exons: Optional[list["Exon"]] = None,
    ):
        self.name = name
        self.version = version
        self.gene = Gene(gene)
        self.tx_position = tx_position
        self.cds_position = cds_position
        self.is_default = is_default
        self.exons = exons if exons else []

    @property
    def full_name(self) -> str:
        if self.version is not None:
            return f"{self.name}.{self.version}"
        else:
            return self.name

    @property
    def is_coding(self) -> bool:
        # Coding transcripts have CDS with non-zero length.
        return self.cds_position.chrom_stop - self.cds_position.chrom_start > 0

    @property
    def strand(self):
        return "+" if self.tx_position.is_forward_strand else "-"

    @property
    def coding_exons(self):
        return [exon.get_as_interval(coding_only=True) for exon in self.exons]


@dataclass
class BED6Interval:
    chrom: str
    chrom_start: int
    chrom_end: int
    name: str
    score: str
    strand: str

    def distance(self, offset: int):
        """Return the distance to the interval.

        if offset is inside the exon, distance is zero.
        otherwise, distance is the distance to the nearest edge.

        distance is positive if the exon comes after the offset.
        distance is negative if the exon comes before the offset.
        """

        start = self.chrom_start + 1
        end = self.chrom_end

        if start <= offset <= end:
            return 0

        start_distance = start - offset
        end_distance = offset - end

        if abs(start_distance) < abs(end_distance):
            return start_distance
        else:
            return -end_distance


class Exon:
    def __init__(self, transcript: Transcript, tx_position: Position, exon_number: int):
        self.transcript = transcript
        self.tx_position = tx_position
        self.exon_number = exon_number
        self.name = f"{self.transcript.name}.{self.exon_number}"

    def get_as_interval(self, coding_only=False) -> Optional[BED6Interval]:
        """Returns the coding region for this exon as a BED6Interval.

        This function returns a BED6Interval objects containing  position
        information for this exon. This may be used as input for
        pybedtools.create_interval_from_list() after casting chrom_start
        and chrom_end as strings.

        coding_only: only include exons in the coding region

        """

        exon_start = self.tx_position.chrom_start
        exon_stop = self.tx_position.chrom_stop

        # Get only exon coding region if requested
        if coding_only:
            if (
                exon_stop <= self.transcript.cds_position.chrom_start
                or exon_start >= self.transcript.cds_position.chrom_stop
            ):
                return None
            exon_start = max(exon_start, self.transcript.cds_position.chrom_start)
            exon_stop = min(
                max(exon_stop, self.transcript.cds_position.chrom_start),
                self.transcript.cds_position.chrom_stop,
            )

        return BED6Interval(
            self.tx_position.chrom,
            exon_start,
            exon_stop,
            self.name,
            ".",
            self.strand,
        )

    @property
    def strand(self):
        strand = "+" if self.tx_position.is_forward_strand else "-"
        return strand


class ChromosomeSubset:
    """
    Allow direct access to a subset of the chromosome.
    """

    def __init__(self, name: str, genome=None):
        self.name = name
        self.genome = genome

    def __getitem__(self, key):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        if isinstance(key, slice):
            start, end = (key.start, key.stop)
            start -= self.genome.start
            end -= self.genome.start
            return self.genome.genome[self.genome.seqid][start:end]
        else:
            raise TypeError(f"Expected a slice object but received a {type(key)}.")

    def __repr__(self):
        return f'ChromosomeSubset("{self.name}")'


class GenomeSubset:
    """
    Allow the direct access of a subset of the genome.
    """

    def __init__(self, genome, chrom, start, end, seqid):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seqid = seqid
        self._chroms = {}

    def __getitem__(self, chrom):
        """Return a chromosome by its name."""
        if chrom in self._chroms:
            return self._chroms[chrom]
        else:
            chromosome = ChromosomeSubset(chrom, self)
            self._chroms[chrom] = chromosome
            return chromosome


class CDNACoord:
    """
    A HGVS cDNA-based coordinate.

    A cDNA coordinate can take one of these forms:

    N = nucleotide N in protein coding sequence (e.g. 11A>G)

    -N = nucleotide N 5' of the ATG translation initiation codon (e.g. -4A>G)
         NOTE: so located in the 5'UTR or 5' of the transcription initiation
         site (upstream of the gene, incl. promoter)

    *N = nucleotide N 3' of the translation stop codon (e.g. *6A>G)
         NOTE: so located in the 3'UTR or 3' of the polyA-addition site
         (including downstream of the gene)

    N+M = nucleotide M in the intron after (3' of) position N in the coding DNA
          reference sequence (e.g. 30+4A>G)

    N-M = nucleotide M in the intron before (5' of) position N in the coding
          DNA reference sequence (e.g. 301-2A>G)

    -N+M / -N-M = nucleotide in an intron in the 5'UTR (e.g. -45+4A>G)

    *N+M / *N-M = nucleotide in an intron in the 3'UTR (e.g. *212-2A>G)
    """

    def __init__(self, coord=0, offset=0, landmark=CDNA_START_CODON, string=""):
        """
        coord: main coordinate along cDNA on the same strand as the transcript

        offset: an additional genomic offset from the main coordinate.  This
                allows referencing non-coding (e.g. intronic) positions.
                Offset is also interpreted on the coding strand.

        landmark: ('cdna_start', 'cdna_stop') indicating that 'coord'
                  is relative to one of these landmarks.

        string: a coordinate from an HGVS name.  If given coord, offset, and
                landmark should not be specified.
        """

        if string:
            if coord != 0 or offset != 0 or landmark != CDNA_START_CODON:
                raise ValueError(
                    "coord, offset, and landmark should not "
                    "be given with string argument"
                )

            self.parse(string)
        else:
            self.coord = coord
            self.offset = offset
            self.landmark = landmark

    def parse(self, coord_text: str) -> "CDNACoord":
        """
        Parse a HGVS formatted cDNA coordinate.
        """

        match = re.match(r"(|-|\*)(\d+)((-|\+)(\d+))?", coord_text)
        if not match:
            raise ValueError(f"unknown coordinate format '{coord_text}'")
        coord_prefix, coord, _, offset_prefix, offset = match.groups()

        self.coord = int(coord)
        self.offset = int(offset) if offset else 0

        if offset_prefix == "-":
            self.offset *= -1
        elif offset_prefix == "+" or offset is None:
            pass
        else:
            raise ValueError(f"unknown offset_prefix '{offset_prefix}'")

        if coord_prefix == "":
            self.landmark = CDNA_START_CODON
        elif coord_prefix == "-":
            self.coord *= -1
            self.landmark = CDNA_START_CODON
        elif coord_prefix == "*":
            self.landmark = CDNA_STOP_CODON
        else:
            raise ValueError(f"unknown coord_prefix '{coord_prefix}'")
        return self

    def __str__(self):
        """
        Return a formatted cDNA coordinate
        """
        if self.landmark == CDNA_STOP_CODON:
            coord_prefix = "*"
        else:
            coord_prefix = ""

        if self.offset < 0:
            offset = str(self.offset)
        elif self.offset > 0:
            offset = "+" + str(self.offset)
        else:
            offset = ""

        return f"{coord_prefix}{self.coord}{offset}"

    def __eq__(self, other):
        return (self.coord, self.offset, self.landmark) == (
            other.coord,
            other.offset,
            other.landmark,
        )

    def __repr__(self):
        """
        Returns a string representation of a cDNA coordinate.
        """
        if self.landmark != CDNA_START_CODON:
            return f"CDNACoord({self.coord}, {self.offset}, '{self.landmark}')"
        else:
            return f"CDNACoord({self.coord}, {self.offset})"
