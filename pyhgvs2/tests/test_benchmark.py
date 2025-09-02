"""
Benchmark tests for HGVSName parsing using pytest-benchmark.
"""
import gzip
import pathlib
from typing import List

import pytest

from .. import HGVSName


def load_benchmark_data() -> List[str]:
    """Load HGVS names from the benchmark data file."""
    benchmark_file = pathlib.Path(__file__).parent.parent / "data" / "benchmark.txt.gz"
    
    with gzip.open(benchmark_file, 'rt') as f:
        # Load all HGVS names from the file
        all_hgvs_names = [line.strip() for line in f if line.strip()]
    
    # Filter out HGVS names with unsupported features (e.g., repetition notation with brackets)
    # Focus on supported HGVS variants for benchmarking
    supported_hgvs_names = []
    
    for hgvs_name in all_hgvs_names:
        # Skip HGVS names with repetition notation [n] which aren't supported yet
        if '[' in hgvs_name or ']' in hgvs_name:
            continue
        
        # Test if the HGVS name can be parsed successfully
        try:
            HGVSName(hgvs_name)
            supported_hgvs_names.append(hgvs_name)
        except Exception:
            # Skip HGVS names that can't be parsed
            continue
    
    return supported_hgvs_names


# Load the benchmark data once at module level to avoid loading it repeatedly
BENCHMARK_DATA = load_benchmark_data()


def test_benchmark_hgvs_parse_all(benchmark):
    """Benchmark parsing all supported HGVS names from the benchmark file."""
    def parse_all_hgvs():
        results = []
        for hgvs_name in BENCHMARK_DATA:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results
    
    # Run the benchmark
    result = benchmark(parse_all_hgvs)
    
    # Verify we parsed the expected number of HGVS names
    assert len(result) == len(BENCHMARK_DATA)
    # Note: This will be fewer than 100,000 since we filter out unsupported variants


def test_benchmark_hgvs_parse_sample(benchmark):
    """Benchmark parsing a representative sample of HGVS names."""
    # Use every 10th HGVS name for a faster benchmark (since we already filtered for supported ones)
    sample_data = BENCHMARK_DATA[::10]
    
    def parse_sample_hgvs():
        results = []
        for hgvs_name in sample_data:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results
    
    # Run the benchmark
    result = benchmark(parse_sample_hgvs)
    
    # Verify we parsed the expected number of HGVS names
    assert len(result) == len(sample_data)
    assert len(result) > 0  # Should have at least some variants


def test_benchmark_hgvs_format(benchmark):
    """Benchmark formatting HGVS names."""
    # Pre-parse a sample of HGVS names (use fewer since we have fewer total)
    sample_data = BENCHMARK_DATA[::100] if len(BENCHMARK_DATA) > 100 else BENCHMARK_DATA
    parsed_hgvs = [HGVSName(name) for name in sample_data]
    
    def format_hgvs():
        results = []
        for hgvs in parsed_hgvs:
            formatted = hgvs.format()
            results.append(formatted)
        return results
    
    # Run the benchmark
    result = benchmark(format_hgvs)
    
    # Verify we formatted the expected number of HGVS names
    assert len(result) == len(parsed_hgvs)
    assert len(result) > 0  # Should have at least some variants


def test_benchmark_hgvs_parse_and_format(benchmark):
    """Benchmark the complete parse-and-format cycle."""
    # Use a smaller sample for the complete cycle
    sample_data = BENCHMARK_DATA[::100] if len(BENCHMARK_DATA) > 100 else BENCHMARK_DATA
    
    def parse_and_format():
        results = []
        for hgvs_name in sample_data:
            parsed = HGVSName(hgvs_name)
            formatted = parsed.format()
            results.append((parsed, formatted))
        return results
    
    # Run the benchmark
    result = benchmark(parse_and_format)
    
    # Verify the results
    assert len(result) == len(sample_data)
    
    # Verify that parse-format round trip works correctly
    for i, (parsed, formatted) in enumerate(result):
        # The formatted version should be parseable
        reparsed = HGVSName(formatted)
        assert reparsed.kind == parsed.kind
        assert reparsed.mutation_type == parsed.mutation_type


@pytest.mark.parametrize("variant_type", ["substitution", "deletion", "insertion"])
def test_benchmark_hgvs_by_type(benchmark, variant_type):
    """Benchmark parsing by variant type."""
    # Filter HGVS names by type for more targeted benchmarks
    if variant_type == "substitution":
        filtered_data = [name for name in BENCHMARK_DATA if ">" in name]
    elif variant_type == "deletion":
        filtered_data = [name for name in BENCHMARK_DATA if "del" in name]
    elif variant_type == "insertion":
        filtered_data = [name for name in BENCHMARK_DATA if "ins" in name]
    else:
        filtered_data = BENCHMARK_DATA
    
    # Use a sample for performance
    sample_size = min(1000, len(filtered_data))
    if sample_size > 0:
        step = max(1, len(filtered_data) // sample_size)
        filtered_data = filtered_data[::step][:sample_size]
    
    # Skip if no variants of this type found
    if not filtered_data:
        pytest.skip(f"No {variant_type} variants found in sample")
    
    def parse_variants():
        results = []
        for hgvs_name in filtered_data:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results
    
    # Run the benchmark
    result = benchmark(parse_variants)
    
    # Verify we parsed some variants
    assert len(result) == len(filtered_data)
    assert len(result) > 0
