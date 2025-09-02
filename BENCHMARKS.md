# HGVSName Performance Benchmarks

This project includes comprehensive performance benchmarks for the `HGVSName` parser using `pytest-benchmark` and real-world HGVS variant data.

## Benchmark Data

The benchmarks use `/pyhgvs2/data/benchmark.txt.gz`, which contains 100,000 real HGVS variant names. The benchmark tests filter out unsupported variant types and currently benchmark ~98,839 supported variants (98.8% coverage).

## Running Benchmarks

### Local Development

1. **Install dependencies:**
   ```bash
   uv sync --dev
   ```

2. **Run all benchmarks:**
   ```bash
   python scripts/run_benchmarks.py
   ```

3. **Run specific benchmark tests:**
   ```bash
   uv run pytest pyhgvs2/tests/test_benchmark.py::test_benchmark_hgvs_parse_sample --benchmark-only
   ```

4. **Run benchmarks with JSON output:**
   ```bash
   python scripts/run_benchmarks.py --save-json --verbose
   ```

5. **Generate benchmark histogram:**
   ```bash
   python scripts/run_benchmarks.py --format histogram
   ```

### CI/CD

Benchmarks automatically run in GitHub Actions on:
- Pull requests
- Manual workflow dispatch

The benchmark job:
- Runs on Python 3.11 (Ubuntu latest)
- Executes all benchmark tests
- Uploads results as artifacts
- Reports performance metrics in the workflow logs

## Benchmark Tests

| Test | Description | Sample Size |
|------|-------------|-------------|
| `test_benchmark_hgvs_parse_all` | Parse all supported HGVS names | ~98,839 variants |
| `test_benchmark_hgvs_parse_sample` | Parse representative sample | ~9,884 variants (every 10th) |
| `test_benchmark_hgvs_format` | Format parsed HGVS names | ~989 variants |
| `test_benchmark_hgvs_parse_and_format` | Complete parse-format cycle | ~989 variants |
| `test_benchmark_hgvs_by_type[substitution]` | Parse substitution variants | Variable (filtered) |
| `test_benchmark_hgvs_by_type[deletion]` | Parse deletion variants | Variable (filtered) |
| `test_benchmark_hgvs_by_type[insertion]` | Parse insertion variants | Variable (filtered) |

## Performance Expectations

Based on current benchmarks (as of September 2025):

- **Parsing speed**: ~16-20 operations per second for full dataset
- **Individual variant parsing**: ~1,975 operations per second
- **Format speed**: ~1,975 operations per second
- **Parse-format cycle**: ~156 operations per second

## Benchmark Configuration

Key benchmark settings (configured in `pyproject.toml`):

```toml
[tool.pytest.benchmark]
group_by = "name"
sort = "mean"
min_rounds = 5
max_time = 2.0
timer = "time.perf_counter"
disable_gc = false
warmup = false
```

## Adding New Benchmarks

To add new benchmark tests:

1. Create test functions in `pyhgvs2/tests/test_benchmark.py`
2. Use the `benchmark` fixture provided by `pytest-benchmark`
3. Follow naming convention: `test_benchmark_*`
4. Include assertions to verify correctness
5. Use appropriate sample sizes for reasonable execution time

Example:
```python
def test_benchmark_new_feature(benchmark):
    """Benchmark a new feature."""
    sample_data = BENCHMARK_DATA[::100]  # Use sample for speed
    
    def run_feature():
        results = []
        for hgvs_name in sample_data:
            # Your benchmark code here
            result = some_operation(hgvs_name)
            results.append(result)
        return results
    
    # Run the benchmark
    result = benchmark(run_feature)
    
    # Verify results
    assert len(result) == len(sample_data)
```

## Troubleshooting

### Common Issues

1. **Benchmark fails with parsing errors:**
   - Check if new HGVS variant types need to be filtered out
   - Update the `load_benchmark_data()` function if needed

2. **Benchmarks take too long:**
   - Reduce sample sizes by increasing the step size (e.g., `[::100]` instead of `[::10]`)
   - Adjust `max_time` in benchmark configuration

3. **CI benchmark job fails:**
   - Check if all dependencies are properly installed
   - Verify the benchmark data file is accessible

### Performance Regression Detection

Monitor benchmark results over time:
- Significant performance decreases (>200%) should trigger investigation
- Compare results between branches/commits
- Use benchmark artifacts to track trends

## Contributing

When contributing code that might affect performance:
1. Run benchmarks locally before submitting PRs
2. Include benchmark results in PR descriptions for significant changes
3. Add new benchmark tests for new parsing features
4. Update this documentation if benchmark structure changes
