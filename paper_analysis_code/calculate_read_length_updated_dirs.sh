#!/bin/bash
#SBATCH --job-name=calculate_read_length
#SBATCH --output=calculate_read_length_%A_%a.out
#SBATCH --error=calculate_read_length_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-40

module load samtools

# Define base and output directories
BASE_DIR="/scratch/sberger/pmgrc_lr_data/duoNovoInputs"
OUTPUT_DIR="/home/lboukas/PMGRC_longread/duo/duoNovo_outputs/coverage_length_stats"

# List of directories
directories=(
  "24-435383" "24-435531" "24-435607" "24-435774" "24-441740" "24-459581"
  "UCI-014" "UCI-025" "UCI-032" "UCI-047" "UCI-082"
  "24-435405" "24-435551" "24-435640" "24-435781" "24-441781" "UCI-002"
  "UCI-020" "UCI-026" "UCI-048" "UCI-088"
  "24-435430" "24-435560" "24-435659" "24-441509" "24-441811" "UCI-022"
  "UCI-029" "UCI-038" "UCI-049"
  "24-435509" "24-435590" "24-435685" "24-441639" "24-441864" "UCI-010"
  "UCI-024" "UCI-030" "UCI-041" "UCI-050"
)

# Total number of directories
total_dirs=${#directories[@]}
  
  # Function to calculate average read length excluding supplementary and secondary alignments
  calculate_avg_readlength() {
    local bam_file="$1"
    samtools stats -@ "$SLURM_CPUS_PER_TASK" "$bam_file" | \
    grep ^SN | cut -f 2- | \
    awk -F'\t' '$1 == "average length:" {print $2}'
  }
  
  # Determine the directory for this array task
  dir_index=$((SLURM_ARRAY_TASK_ID - 1))
  
  # Check if the task ID is within the range of directories
  if [ "$dir_index" -ge "$total_dirs" ]; then
  echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID. Exiting."
  exit 1
  fi
  
  dir="${directories[$dir_index]}"
  echo "Processing directory: $dir"
  
  # Locate the samples.txt file
  samples_file=$(find "${BASE_DIR}/${dir}" -maxdepth 1 -type f -name "*samples.txt" | head -n 1)
  echo "Using samples file: $samples_file"
  
  if [[ ! -f "$samples_file" ]]; then
  echo "samples.txt not found in directory: $dir" >&2
  exit 1
  fi
  
  output_file="${OUTPUT_DIR}/${dir}_read_length.tsv"
  echo -e "Sample\tRelationship\tAverage Read Length" > "$output_file"
  
  # Read and process each line in samples_file
  while read -r sample relationship id number coverage gvcf_path bam_path; do
  # Convert BAM path to absolute
  bam_path="/scratch/sberger/pmgrc_lr_data/${bam_path#../}"
  
  if [ -f "$bam_path" ]; then
  avg_read_length=$(calculate_avg_readlength "$bam_path")
  echo -e "${id}\t${relationship}\t${avg_read_length}" >> "$output_file"
  echo "Processed: $bam_path"
  else
    echo "File not found: $bam_path" >&2
  fi
  done < "$samples_file"
  
  echo "Results for ${dir} saved to $output_file"
  
