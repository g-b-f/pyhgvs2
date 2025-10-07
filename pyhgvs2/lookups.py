"""
Given X, get Y
"""

from typing import TYPE_CHECKING, Optional, cast

if TYPE_CHECKING:
    from .hgvsg_name import HGVSName

from .constants import CDNA_START_CODON, CDNA_STOP_CODON
from .models import BED6Interval, CDNACoord, Exon, GenomeType, Position, Transcript

# The RefSeq standard for naming contigs/transcripts/proteins:
# http://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly  # nopep8
REFSEQ_PREFIXES = [
    ("AC_", "genomic", "Complete genomic molecule, usually alternate assembly"),
    ("NC_", "genomic", "Complete genomic molecule, usually reference assembly"),
    ("NG_", "genomic", "Incomplete genomic region"),
    ("NT_", "genomic", "Contig or scaffold, clone-based or WGS"),
    ("NW_", "genomic", "Contig or scaffold, primarily WGS"),
    ("NS_", "genomic", "Environmental sequence"),
    ("NZ_", "genomic", "Unfinished WGS"),
    ("NM_", "mRNA", ""),
    ("NR_", "RNA", ""),
    ("XM_", "mRNA", "Predicted model"),
    ("XR_", "RNA", "Predicted model"),
    ("AP_", "Protein", "Annotated on AC_ alternate assembly"),
    ("NP_", "Protein", "Associated with an NM_ or NC_ accession"),
    ("YP_", "Protein", ""),
    ("XP_", "Protein", "Predicted model, associated with an XM_ accession"),
    ("ZP_", "Protein", "Predicted model, annotated on NZ_ genomic records"),
]

REFSEQ_PREFIX_LOOKUP = {prefix: (kind, description) for prefix, kind, description in REFSEQ_PREFIXES}


def get_refseq_type(name: str) -> Optional[str]:
    """
    Return the RefSeq type for a refseq name.
    """
    prefix = name[:3]
    return REFSEQ_PREFIX_LOOKUP.get(prefix, (None, ""))[0]


def get_exons(transcript: Transcript) -> list[Exon]:
    """Yield exons in coding order."""
    transcript_strand = transcript.tx_position.is_forward_strand
    if hasattr(transcript.exons, "select_related"):
        exons = list(transcript.exons.select_related("tx_position"))
    else:
        exons = list(transcript.exons)
    exons.sort(key=lambda exon: exon.tx_position.chrom_start)
    if not transcript_strand:
        exons.reverse()
    return exons


def get_coding_exons(transcript: Transcript):
    """Yield non-empty coding exonic regions in coding order."""
    for exon in get_exons(transcript):
        region = exon.get_as_interval(coding_only=True)
        if region:
            yield region


def get_utr5p_size(transcript: Transcript):
    """Return the size of the 5prime UTR of a transcript."""

    transcript_strand = transcript.tx_position.is_forward_strand
    exons = get_exons(transcript)

    # Find the exon containing the start codon.
    start_codon = (
        transcript.cds_position.chrom_start
        if transcript_strand
        else transcript.cds_position.chrom_stop - 1
    )
    cdna_len = 0
    for exon in exons:
        exon_start = exon.tx_position.chrom_start
        exon_end = exon.tx_position.chrom_stop
        if exon_start <= start_codon < exon_end:
            break
        cdna_len += exon_end - exon_start
    else:
        raise ValueError("transcript contains no exons")

    if transcript_strand:
        return cdna_len + (start_codon - exon_start)
    else:
        return cdna_len + (exon_end - start_codon - 1)


def find_stop_codon(exons: list["Exon"], cds_position: Position):
    """Return the position along the cDNA of the base after the stop codon."""
    if cds_position.is_forward_strand:
        stop_pos = cds_position.chrom_stop
    else:
        stop_pos = cds_position.chrom_start
    cdna_pos = 0
    for exon in exons:
        exon_start = exon.tx_position.chrom_start
        exon_stop = exon.tx_position.chrom_stop

        if exon_start <= stop_pos <= exon_stop:
            if cds_position.is_forward_strand:
                return cdna_pos + stop_pos - exon_start
            else:
                return cdna_pos + exon_stop - stop_pos
        else:
            cdna_pos += exon_stop - exon_start
    raise ValueError("Stop codon is not in any of the exons")


def get_genomic_sequence(genome: GenomeType, chrom: str, start: int, end: int):
    """
    Return a sequence for the genomic region.

    start, end: 1-based, end-inclusive coordinates of the sequence.
    """
    if start > end:
        return ""
    else:
        return str(genome[str(chrom)][start - 1 : end]).upper()


def cdna_to_genomic_coord(transcript: Transcript, coord: CDNACoord):
    """Convert a HGVS cDNA coordinate to a genomic coordinate."""
    transcript_strand = transcript.tx_position.is_forward_strand
    exons = get_exons(transcript)
    utr5p = get_utr5p_size(transcript) if transcript.is_coding else 0

    # compute starting position along spliced transcript.
    if coord.landmark == CDNA_START_CODON:
        if coord.coord > 0:
            pos = utr5p + coord.coord
        else:
            pos = utr5p + coord.coord + 1
    elif coord.landmark == CDNA_STOP_CODON:
        if coord.coord < 0:
            raise ValueError(
                "CDNACoord cannot have a negative coord and landmark CDNA_STOP_CODON"
            )
        pos = find_stop_codon(exons, transcript.cds_position) + coord.coord
        if not transcript.is_coding:
            if not transcript_strand:
                pos = exons[-1].tx_position.chrom_start
                return pos - coord.coord + 1
    else:
        raise ValueError(f'unknown CDNACoord landmark "{coord.landmark}"')

    # 5' flanking sequence.
    if pos < 1:
        if transcript_strand:
            return transcript.tx_position.chrom_start + pos
        else:
            return transcript.tx_position.chrom_stop - pos + 1

    # Walk along transcript until we find an exon that contains pos.
    cdna_start = 1
    cdna_end = 1
    for exon in exons:
        exon_start = exon.tx_position.chrom_start + 1
        exon_end = exon.tx_position.chrom_stop
        cdna_end = cdna_start + (exon_end - exon_start)
        if cdna_start <= pos <= cdna_end:
            break
        cdna_start = cdna_end + 1
    else:
        # 3' flanking sequence
        if transcript_strand:
            return transcript.cds_position.chrom_stop + coord.coord
        else:
            return transcript.cds_position.chrom_start + 1 - coord.coord

    # Compute genomic coordinate using offset.
    if transcript_strand:
        # Plus strand.
        return exon_start + (pos - cdna_start) + coord.offset
    else:
        # Minus strand.
        return exon_end - (pos - cdna_start) - coord.offset


def genomic_to_cdna_coord(
    transcript: Transcript, genomic_coord: int
) -> Optional[CDNACoord]:
    """Convert a genomic coordinate to a cDNA coordinate and offset."""
    exons = cast(
        list[BED6Interval], [exon.get_as_interval() for exon in get_exons(transcript)]
    )

    if len(exons) == 0:
        return None

    strand = transcript.strand

    if strand == "+":
        exons.sort(key=lambda exon: exon.chrom_start)
    else:
        exons.sort(key=lambda exon: -exon.chrom_end)

    distances = [exon.distance(genomic_coord) for exon in exons]
    min_distance_to_exon = min(map(abs, distances))

    coding_offset = 0
    for exon in exons:
        assert exon is not None
        exon_length = exon.chrom_end - exon.chrom_start
        distance = exon.distance(genomic_coord)
        if abs(distance) == min_distance_to_exon:
            if strand == "+":
                exon_start_cds_offset = coding_offset + 1
                exon_end_cds_offset = coding_offset + exon_length
            else:
                exon_start_cds_offset = coding_offset + exon_length
                exon_end_cds_offset = coding_offset + 1
            # This is the exon we want to annotate against.
            if distance == 0:
                # Inside the exon.
                if strand == "+":
                    coord = exon_start_cds_offset + (
                        genomic_coord - (exon.chrom_start + 1)
                    )
                else:
                    coord = exon_end_cds_offset + (exon.chrom_end - genomic_coord)
                cdna_coord = CDNACoord(coord, 0)
            else:
                # Outside the exon.
                if distance > 0:
                    nearest_exonic = exon_start_cds_offset
                else:
                    nearest_exonic = exon_end_cds_offset
                if strand == "+":
                    distance *= -1

                # If outside transcript, don't use offset.
                if (
                    genomic_coord < transcript.tx_position.chrom_start + 1
                    or genomic_coord > transcript.tx_position.chrom_stop
                ):
                    nearest_exonic += distance
                    distance = 0
                cdna_coord = CDNACoord(nearest_exonic, distance)
            break
        coding_offset += exon_length

    # Adjust coordinates for coding transcript.
    if transcript.is_coding:
        # Detect if position before start codon.
        utr5p = get_utr5p_size(transcript) if transcript.is_coding else 0
        cdna_coord.coord -= utr5p
        if cdna_coord.coord <= 0:
            cdna_coord.coord -= 1
        else:
            # Detect if position is after stop_codon.
            exons = get_exons(transcript)
            stop_codon = find_stop_codon(exons, transcript.cds_position)
            stop_codon -= utr5p
            if (
                cdna_coord.coord > stop_codon
                or cdna_coord.coord == stop_codon
                and cdna_coord.offset > 0
            ):
                cdna_coord.coord -= stop_codon
                cdna_coord.landmark = CDNA_STOP_CODON

    else:  # non coding
        if strand == "+":
            # Detect if position is after last exon.
            if genomic_coord > exons[-1].chrom_end:
                cdna_coord.coord = genomic_coord - exons[-1].chrom_end
                cdna_coord.landmark = CDNA_STOP_CODON
            else:
                # Detect if position is before first exon.
                if genomic_coord <= exons[0].chrom_start:
                    cdna_coord.coord -= 1
        else:  # neg strand
            # Detect if position is after last exon.
            if genomic_coord <= exons[-1].chrom_start:
                cdna_coord.coord = exons[-1].chrom_start - genomic_coord + 1
                cdna_coord.landmark = CDNA_STOP_CODON
            else:
                # Detect if position is before first exon.
                if genomic_coord >= exons[0].chrom_end:
                    cdna_coord.coord -= 1

    return cdna_coord


def get_allele(hgvs, genome, transcript=None):
    """Get an allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True
    )
    ref = get_genomic_sequence(genome, chrom, start, end)
    return chrom, start, end, ref, alt


_indel_mutation_types = {"ins", "del", "dup", "delins"}


def get_vcf_allele(
    hgvs: "HGVSName", genome: GenomeType, transcript: Optional[Transcript] = None
) -> tuple[str, int, int, str, str]:
    """Get an VCF-style allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_vcf_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True
    )
    ref = get_genomic_sequence(genome, chrom, start, end)

    if hgvs.mutation_type in _indel_mutation_types:
        # Left-pad alternate allele.
        alt = ref[0] + alt
    return chrom, start, end, ref, alt


def matches_ref_allele(hgvs, genome, transcript=None):
    """Return True if reference allele matches genomic sequence."""
    ref, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True
    )
    chrom, start, end = hgvs.get_coords(transcript)
    genome_ref = get_genomic_sequence(genome, chrom, start, end)
    return genome_ref == ref
