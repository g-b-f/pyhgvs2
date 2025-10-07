from typing import Callable, Literal, Optional

from .hgvsg_name import HGVSName
from .lookups import genomic_to_cdna_coord, get_vcf_allele
from .models import GenomeSubset, GenomeType, Transcript
from .variants import justify_indel, normalize_variant, revcomp


def hgvs_justify_dup(
    chrom: str, offset: int, ref: str, alt: str, genome: GenomeType
) -> tuple[str, int, str, str, str]:
    """
    Determines if allele is a duplication and justifies.

    chrom: Chromosome name.
    offset: 1-index genomic coordinate.
    ref: Reference allele (no padding).
    alt: Alternate allele (no padding).
    genome: pygr compatible genome object.

    Returns duplicated region [start, end] if allele is an insert that
    could be represented as a duplication. Otherwise, returns None.
    """

    if len(ref) == len(alt) == 0:
        # it's a SNP, just return.
        return chrom, offset, ref, alt, ">"

    if len(ref) > 0 and len(alt) > 0:
        # complex indel, don't know how to dup check
        return chrom, offset, ref, alt, "delins"

    if len(ref) > len(alt):
        # deletion -- don't dup check
        return chrom, offset, ref, alt, "del"

    indel_seq = alt
    indel_length = len(indel_seq)

    # Convert offset to 0-index.
    offset -= 1

    # Get genomic sequence around the lesion.
    prev_seq = str(genome[str(chrom)][offset - indel_length : offset]).upper()
    next_seq = str(genome[str(chrom)][offset : offset + indel_length]).upper()

    # Convert offset back to 1-index.
    offset += 1

    if prev_seq == indel_seq:
        offset = offset - indel_length
        mutation_type = "dup"
        ref = indel_seq
        alt = indel_seq * 2
    elif next_seq == indel_seq:
        mutation_type = "dup"
        ref = indel_seq
        alt = indel_seq * 2
    else:
        mutation_type = "ins"

    return chrom, offset, ref, alt, mutation_type


def hgvs_justify_indel(
    chrom: str, offset: int, ref: str, alt: str, strand: str, genome: GenomeType
) -> tuple[str, int, str, str]:
    """
    3' justify an indel according to the HGVS standard.

    Returns (chrom, offset, ref, alt).
    """
    if len(ref) == len(alt) == 0:
        # It's a SNP, just return.
        return chrom, offset, ref, alt

    if len(ref) > 0 and len(alt) > 0:
        # Complex indel, don't know how to justify.
        return chrom, offset, ref, alt

    # Get genomic sequence around the lesion.
    start = max(offset - 100, 0)
    end = offset + 100
    seq = str(genome[str(chrom)][start - 1 : end]).upper()
    cds_offset = offset - start

    # indel -- strip off the ref base to get the actual lesion sequence
    is_insert = len(alt) > 0
    if is_insert:
        indel_seq = alt
        cds_offset_end = cds_offset
    else:
        indel_seq = ref
        cds_offset_end = cds_offset + len(indel_seq)

    # Now 3' justify (vs. cDNA not genome) the offset
    justify: Literal["left", "right"] = "right" if strand == "+" else "left"
    offset, _, indel_seq = justify_indel(
        cds_offset, cds_offset_end, indel_seq, seq, justify
    )
    offset += start

    if is_insert:
        alt = indel_seq
    else:
        ref = indel_seq

    return chrom, offset, ref, alt


def hgvs_normalize_variant(
    chrom: str,
    offset: int,
    ref: str,
    alt: str,
    genome: GenomeType,
    transcript: Optional[Transcript] = None,
):
    """Convert VCF-style variant to HGVS-style."""
    if len(ref) == len(alt) == 1:
        if ref == alt:
            mutation_type = "="
        else:
            mutation_type = ">"
    else:
        # Remove 1bp padding
        offset += 1
        ref = ref[1:]
        alt = alt[1:]

        # 3' justify allele.
        strand = transcript.strand if transcript else "+"
        chrom, offset, ref, alt = hgvs_justify_indel(
            chrom, offset, ref, alt, strand, genome
        )

        # Represent as duplication if possible.
        chrom, offset, ref, alt, mutation_type = hgvs_justify_dup(
            chrom, offset, ref, alt, genome
        )
    return chrom, offset, ref, alt, mutation_type


def parse_hgvs_name(
    hgvs_name: str,
    genome: GenomeType,
    transcript=None,
    get_transcript: Optional[Callable] = None,
    flank_length=30,
    normalize=True,
    lazy=False,
):
    """
    Parse an HGVS name into (chrom, start, end, ref, alt)

    hgvs_name: HGVS name to parse.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to HGVS name.
    normalize: If True, normalize allele according to VCF standard.
    lazy: If True, discard version information from incoming transcript/gene.
    """
    hgvs = HGVSName(hgvs_name)

    # Determine transcript.
    if hgvs.kind == "c" and not transcript:
        if "." in hgvs.transcript and lazy:
            hgvs.transcript, version = hgvs.transcript.split(".")
        elif "." in hgvs.gene and lazy:
            hgvs.gene, version = hgvs.gene.split(".")
        if get_transcript:
            if hgvs.transcript:
                transcript = get_transcript(hgvs.transcript)
            elif hgvs.gene:
                transcript = get_transcript(hgvs.gene)
        if not transcript:
            raise ValueError("transcript is required")

    if transcript and hgvs.transcript in genome:
        # Reference sequence is directly known, use it.
        genome = GenomeSubset(
            genome,
            transcript.tx_position.chrom,
            transcript.tx_position.chrom_start,
            transcript.tx_position.chrom_stop,
            hgvs.transcript,
        )

    chrom, start, end, ref, alt = get_vcf_allele(hgvs, genome, transcript)
    if normalize:
        chrom, start, ref, [alt] = normalize_variant(
            chrom, start, ref, [alt], genome, flank_length=flank_length
        ).variant
    return (chrom, start, ref, alt)


def variant_to_hgvs_name(
    chrom: str,
    offset: int,
    ref: str,
    alt: str,
    genome: GenomeType,
    transcript: Transcript,
    max_allele_length=4,
) -> HGVSName:
    """
    Populate a HGVSName from a genomic coordinate.

    chrom: Chromosome name.
    offset: Genomic offset of allele.
    ref: Reference allele.
    alt: Alternate allele.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to allele.
    max_allele_length: If allele is greater than this use allele length.
    """
    # Convert VCF-style variant to HGVS-style.
    chrom, offset, ref, [alt] = normalize_variant(
        chrom, offset, ref, [alt], genome
    ).variant
    chrom, offset, ref, alt, mutation_type = hgvs_normalize_variant(
        chrom, offset, ref, alt, genome, transcript
    )

    # Populate HGVSName parse tree.
    hgvs = HGVSName()

    # Populate coordinates.
    if not transcript:
        # Use genomic coordinate when no transcript is available.
        hgvs.kind = "g"
        hgvs.start = offset
        hgvs.end = offset + len(ref) - 1
    else:
        # Use cDNA coordinates.
        hgvs.kind = "c"
        if mutation_type == "ins":
            # Always use a range for insertions (HGVS: c.306_307insG)
            offset_start = offset - 1
            offset_end = offset
            if transcript.strand == "-":
                offset_start, offset_end = offset_end, offset_start
            hgvs.cdna_start = genomic_to_cdna_coord(transcript, offset_start)
            hgvs.cdna_end = genomic_to_cdna_coord(transcript, offset_end)
        else:
            is_single_base_indel = (
                mutation_type in ("del", "delins", "dup") and len(ref) == 1
            ) or mutation_type == ">"
            if mutation_type == ">" or is_single_base_indel:
                # Use a single coordinate.
                hgvs.cdna_start = genomic_to_cdna_coord(transcript, offset)
                hgvs.cdna_end = hgvs.cdna_start
            else:
                # Use a range of coordinates.
                offset_start = offset
                offset_end = offset + len(ref) - 1
                if transcript.strand == "-":
                    offset_start, offset_end = offset_end, offset_start
                hgvs.cdna_start = genomic_to_cdna_coord(transcript, offset_start)
                hgvs.cdna_end = genomic_to_cdna_coord(transcript, offset_end)

    # Populate prefix.
    if transcript:
        hgvs.transcript = transcript.full_name
        hgvs.gene = transcript.gene.name

    # Convert alleles to transcript strand.
    if transcript and transcript.strand == "-":
        ref = revcomp(ref)
        alt = revcomp(alt)

    # Convert to allele length if alleles are long.
    ref_len = len(ref)
    alt_len = len(alt)
    if (mutation_type == "dup" and ref_len > max_allele_length) or (
        mutation_type != "dup"
        and (ref_len > max_allele_length or alt_len > max_allele_length)
    ):
        ref = str(ref_len)
        alt = str(alt_len)

    # Populate alleles.
    hgvs.mutation_type = mutation_type
    hgvs.ref_allele = ref
    hgvs.alt_allele = alt

    return hgvs


def format_hgvs_name(
    chrom: str,
    offset: int,
    ref: str,
    alt: str,
    genome: GenomeType,
    transcript: Transcript,
    use_prefix=True,
    use_gene=True,
    use_counsyl=False,
    max_allele_length=4,
) -> str:
    """
    Generate a HGVS name from a genomic coordinate.

    chrom: Chromosome name.
    offset: Genomic offset of allele.
    ref: Reference allele.
    alt: Alternate allele.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to allele.
    use_prefix: Include a transcript/gene/chromosome prefix in HGVS name.
    use_gene: Include gene name in HGVS prefix.
    max_allele_length: If allele is greater than this use allele length.
    """
    hgvs = variant_to_hgvs_name(
        chrom,
        offset,
        ref,
        alt,
        genome,
        transcript,
        max_allele_length=max_allele_length,
    )
    return hgvs.format(
        use_prefix=use_prefix,
        use_gene=use_gene,
    )
