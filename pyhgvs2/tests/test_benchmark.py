import gzip
from pathlib import Path

from pytest_benchmark.fixture import BenchmarkFixture

from .. import HGVSName

benchmark_file = Path(__file__).parent / "data" / "benchmark_filtered.txt.gz"
with gzip.open(benchmark_file, "rt") as f:
    BENCHMARK_DATA = [line.strip() for line in f if line.strip()]


def test_benchmark_hgvs_parse_sample(benchmark: BenchmarkFixture):
    def parse_sample_hgvs():
        results = []
        for hgvs_name in BENCHMARK_DATA:
            parsed = HGVSName(hgvs_name)
            results.append(parsed)
        return results

    result = benchmark(parse_sample_hgvs)
    assert len(result) == len(BENCHMARK_DATA)
