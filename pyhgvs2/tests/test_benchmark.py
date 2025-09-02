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


def test_benchmark_hgvs_parse_sample(benchmark):
    """Benchmark parsing a representative sample of HGVS names."""
    # Use every 100th HGVS name for a faster benchmark that still catches regressions
    sample_data = BENCHMARK_DATA[::100]
    
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


def test_benchmark_hgvs_parse_and_format(benchmark):
    """Benchmark the complete parse-and-format cycle."""
    # Use a smaller sample for the complete cycle - this catches regressions in both parse and format
    sample_data = BENCHMARK_DATA[::500] if len(BENCHMARK_DATA) > 500 else BENCHMARK_DATA[:200]
    
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
