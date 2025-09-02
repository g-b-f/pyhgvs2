"""
Benchmark tests for HGVSName parsing using pytest-benchmark.
"""
import gzip
import pathlib
from typing import List

import pytest

from .. import HGVSName


benchmark_file = pathlib.Path(__file__).parent.parent / "tests"/"data" / "benchmark_filtered.txt.gz"
with gzip.open(benchmark_file, 'rt') as f:
    BENCHMARK_DATA = [line.strip() for line in f if line.strip()]


def test_benchmark_hgvs_parse_sample(benchmark):

    def parse_sample_hgvs():
        results = []
        for hgvs_name in BENCHMARK_DATA:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results
    
    result = benchmark(parse_sample_hgvs)
    assert len(result) == len(BENCHMARK_DATA)