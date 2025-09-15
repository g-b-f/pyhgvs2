import gzip
from pathlib import Path

from pytest_benchmark.fixture import BenchmarkFixture

from .. import HGVSName

benchmark_file = Path(__file__).parent / "data" / "benchmark_filtered.txt.gz"
with gzip.open(benchmark_file, "rt") as f:
    BENCHMARK_DATA = f.readlines()


def test_benchmark_hgvsName(benchmark: BenchmarkFixture):
    def parse_hgvs():
        results = []
        for hgvs_name in BENCHMARK_DATA:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results

    result = benchmark(parse_hgvs)
    assert len(result) == len(BENCHMARK_DATA)
