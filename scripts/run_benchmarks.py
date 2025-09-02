#!/usr/bin/env python3
"""
Benchmark runner script for HGVSName performance testing.

This script can be used to run benchmarks locally or in CI with various configurations.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def run_benchmarks(
    benchmark_file: str = "pyhgvs2/tests/test_benchmark.py",
    output_format: str = "table",
    save_json: bool = False,
    min_rounds: int = 5,
    verbose: bool = False,
):
    """Run pytest benchmarks with specified options."""
    
    # Ensure we're in the project root directory
    project_root = Path(__file__).parent.parent
    benchmark_path = project_root / benchmark_file
    
    if not benchmark_path.exists():
        print(f"❌ Benchmark file not found: {benchmark_path}")
        sys.exit(1)
    
    cmd = [
        sys.executable, "-m", "pytest",
        str(benchmark_path),
        "--benchmark-only",
        f"--benchmark-min-rounds={min_rounds}",
        "--benchmark-sort=name",
    ]
    
    if verbose:
        cmd.append("--benchmark-verbose")
    
    if save_json:
        cmd.extend(["--benchmark-json=benchmark.json"])
    
    if output_format == "histogram":
        cmd.append("--benchmark-histogram=histogram")
    
    # Run the benchmark
    print(f"Running benchmarks with command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=project_root)
    
    if result.returncode == 0:
        print("\n✅ Benchmarks completed successfully!")
        if save_json:
            print("📊 Results saved to benchmark.json")
        if output_format == "histogram":
            print("📈 Histogram saved to histogram.svg")
    else:
        print("\n❌ Benchmarks failed!")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Run HGVSName performance benchmarks")
    parser.add_argument(
        "--file",
        default="pyhgvs2/tests/test_benchmark.py",
        help="Benchmark test file to run (default: pyhgvs2/tests/test_benchmark.py)"
    )
    parser.add_argument(
        "--format",
        choices=["table", "histogram"],
        default="table",
        help="Output format (default: table)"
    )
    parser.add_argument(
        "--save-json",
        action="store_true",
        help="Save benchmark results to JSON file"
    )
    parser.add_argument(
        "--min-rounds",
        type=int,
        default=5,
        help="Minimum number of benchmark rounds (default: 5)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose benchmark output"
    )
    
    args = parser.parse_args()
    
    run_benchmarks(
        benchmark_file=args.file,
        output_format=args.format,
        save_json=args.save_json,
        min_rounds=args.min_rounds,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
