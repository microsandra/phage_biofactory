#!/bin/bash
#$ -N s05_depth_normalization
#$ -cwd
#$ -V
#$ -j y
#$ -o logs/
#$ -t 1-6

# Submit the count array first, then submit the sample array after all counts finish.
#   mkdir -p logs
#   qsub s05_depth_normalization.sh count
#   qsub s05_depth_normalization.sh sample
set -eo pipefail

task=${SGE_TASK_ID:-${2:-}}
if [[ ! $task =~ ^[1-6]$ || ! ${1:-} =~ ^(count|sample)$ ]]; then
    echo "Usage: qsub s05_depth_normalization.sh {count|sample}" >&2
    exit 1
fi

python3 - "$1" "$task" "${CLEAN_DIR:-05_unmapped_reads_microbiome}" \
    "${EQUAL_DEPTH_DIR:-taxonomy_equal_depth_reads}" <<'PY'
import gzip
from itertools import zip_longest
import os
from pathlib import Path
import random
import sys
import tempfile

mode, task, clean_dir, output_dir = sys.argv[1:]
samples = [f"Larva_{i:02d}" for i in range(1, 7)]
sample = samples[int(task) - 1]
clean_dir, output_dir = Path(clean_dir), Path(output_dir)
min_length = 70
seed = 22


def input_paths():
    for suffix in (".fastq", ".fastq.gz"):
        paths = tuple(clean_dir / f"{sample}_clean_{mate}{suffix}" for mate in (1, 2))
        if all(path.is_file() for path in paths):
            return paths
    raise SystemExit(f"{sample}: cannot find both clean FASTQs in {clean_dir}")


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        for number, header in enumerate(handle, 1):
            sequence, separator, quality = (handle.readline() for _ in range(3))
            length = len(sequence.rstrip(b"\r\n"))
            if (not header.startswith(b"@") or not separator.startswith(b"+")
                    or not quality or length != len(quality.rstrip(b"\r\n"))):
                raise SystemExit(f"{path}: invalid FASTQ record {number}")
            read_id = header.split(maxsplit=1)[0]
            if read_id.endswith((b"/1", b"/2")):
                read_id = read_id[:-2]
            yield read_id, length, (header, sequence, separator, quality)


def eligible_pairs():
    for number, (r1, r2) in enumerate(
        zip_longest(*(records(path) for path in input_paths())), 1
    ):
        if r1 is None or r2 is None or r1[0] != r2[0]:
            raise SystemExit(f"{sample}: unmatched R1/R2 at record {number}")
        if min(r1[1], r2[1]) >= min_length:
            yield r1[2], r2[2]


counts_dir = output_dir / "eligible_pair_counts"
if mode == "count":
    counts_dir.mkdir(parents=True, exist_ok=True)
    count = sum(1 for _ in eligible_pairs())
    if not count:
        raise SystemExit(f"{sample}: no eligible read pairs")
    with tempfile.NamedTemporaryFile(mode="w", dir=counts_dir, delete=False) as out:
        out.write(f"{count}\n")
    os.replace(out.name, counts_dir / f"{sample}.txt")
    print(f"{sample}: {count:,} pairs with both mates at least {min_length} bases")
else:
    try:
        counts = {
            name: int((counts_dir / f"{name}.txt").read_text()) for name in samples
        }
    except FileNotFoundError as error:
        raise SystemExit(f"Finish all six count jobs first. Missing {error.filename}")
    target = min(counts.values())
    if target < 1:
        raise SystemExit("The common depth must be positive")
    destination = output_dir / sample
    if destination.exists():
        raise SystemExit(f"Output already exists: {destination}")

    # Choose uniformly in one pass while keeping the mates together.
    remaining, needed = counts[sample], target
    rng = random.Random(seed)
    with tempfile.TemporaryDirectory(prefix=f".{sample}.", dir=output_dir) as temp:
        staging = Path(temp) / sample
        staging.mkdir()
        with gzip.open(staging / f"{sample}_R1.fastq.gz", "wb", compresslevel=1) as out1, \
             gzip.open(staging / f"{sample}_R2.fastq.gz", "wb", compresslevel=1) as out2:
            for r1, r2 in eligible_pairs():
                if remaining <= 0:
                    raise SystemExit(f"{sample}: eligible pair count changed. Repeat the count step.")
                if needed == remaining or (needed and rng.randrange(remaining) < needed):
                    out1.writelines(r1)
                    out2.writelines(r2)
                    needed -= 1
                remaining -= 1
        if remaining or needed:
            raise SystemExit(f"{sample}: eligible pair count changed. Repeat the count step.")
        (staging / "pairs.txt").write_text(f"{target}\n")
        staging.rename(destination)
    print(f"{sample}: kept {target:,} of {counts[sample]:,} eligible pairs (seed {seed})")
PY
