import gzip
import json
import pathlib

from pyhgvs2 import HGVSName

data_path = pathlib.Path(__file__).parent.parent / "data"
benchmark_file = data_path / "benchmark.txt.gz"
errors = {}
benchmark_data = []
MAX_LINES = 1000

with gzip.open(benchmark_file, "rt") as f:
    for i, line in enumerate(f):
        if i >= MAX_LINES:
            break
        line = line.strip()
        if not line:
            continue
        try:
            HGVSName(line)
            benchmark_data.append(line)
        except Exception as e:
            errors[line] = str(e)

print(f"skipped {len(errors)} lines with errors")
with open(data_path / "parsing_errors.json", "w") as f:
    json.dump(errors, f, indent=4)


with gzip.open(data_path / "benchmark_filtered.txt.gz", "wt") as f:
    for hgvs_name in benchmark_data:
        f.write(hgvs_name + "\n")
