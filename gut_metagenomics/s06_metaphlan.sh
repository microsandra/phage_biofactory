#!/bin/bash
#$ -N s06_metaphlan
#$ -cwd
#$ -V
#$ -j y
#$ -o logs/
#$ -pe openmpi 8
#$ -t 1-6

# Submit after s05_depth_normalization.sh sample finishes.
set -eo pipefail

task=${SGE_TASK_ID:-${1:-}}
[[ $task =~ ^[1-6]$ ]] || { echo "Submit with qsub or pass a task number (1–6)." >&2; exit 1; }
printf -v sample 'Larva_%02d' "$task"

input_dir="${EQUAL_DEPTH_DIR:-taxonomy_equal_depth_reads}/$sample"
out_dir=${MPA_OUT_DIR:-metaphlan_equal_depth}
: "${MPA_DB:?Export MPA_DB as the directory containing the installed MetaPhlAn database.}"
mpa_index=${MPA_INDEX:-mpa_vJan25_CHOCOPhlAnSGB_202503}
r1="$input_dir/${sample}_R1.fastq.gz"
r2="$input_dir/${sample}_R2.fastq.gz"
[[ -s $r1 && -s $r2 && -s $input_dir/pairs.txt ]] || {
    echo "Missing equal-depth inputs for $sample. Run s05_depth_normalization.sh first." >&2
    exit 1
}
expected_pairs=$(< "$input_dir/pairs.txt")

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${MPA_CONDA_ENV:-metaphlan4_env}"
metaphlan --version
echo "$sample: $expected_pairs input pairs, MetaPhlAn database $mpa_index"

mkdir -p "$out_dir"
sample_out="$out_dir/$sample"
mkdir "$sample_out" || { echo "Output already exists: $sample_out" >&2; exit 1; }
trap 'rm -rf -- "$sample_out"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

mapout="$sample_out/$sample.mapout.tsv"
# Submit both FASTQs without internal subsampling.
metaphlan "$r1,$r2" \
    --input_type fastq \
    --db_dir "$MPA_DB" --index "$mpa_index" --offline \
    --nproc "${NSLOTS:-8}" \
    --read_min_len 70 \
    --mapout "$mapout" \
    -t rel_ab_w_read_stats \
    -o "$sample_out/${sample}_profile.txt"

observed_reads=$(tail -n 2 "$mapout" | awk '$1 == "#nreads" { print $2 }')
expected_reads=$((2 * expected_pairs))
if [[ $observed_reads != "$expected_reads" ]]; then
    echo "$sample: MetaPhlAn processed ${observed_reads:-unknown} reads, expected $expected_reads." >&2
    exit 1
fi
gzip -- "$mapout"
trap - EXIT INT TERM
echo "$sample: verified $observed_reads individual reads processed"
