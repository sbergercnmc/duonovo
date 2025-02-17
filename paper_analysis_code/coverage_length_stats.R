### First run this on the cluster to combine coverage and read length results into a single tsv

### For coverage
combined_coverage="combined_coverage_results.tsv" 
echo -e "Directory\tSample\tRelationship\tAverage Coverage" > "$combined_coverage" 
for file in *_coverage.tsv; do     
[[ "$file" == "$combined_coverage" ]] && continue     
dir="${file%%_*}"     
if [[ -s "$file" ]]; then         
tail -n +2 "$file" | awk -v dir="$dir" 'BEGIN {OFS="\t"} {print dir, $0}' >> "$combined_coverage"     
else         
  echo "Warning: Coverage file '$file' is empty or missing." >&2     
fi 
done

### For read length
combined_read_length="combined_read_length_results.tsv" 
echo -e "Directory\tSample\tRelationship\tAverage Read Length" > "$combined_read_length" 

for file in *_read_length.tsv; do
[[ "$file" == "$combined_read_length" ]] && continue
dir="${file%%_*}"
if [[ -s "$file" ]]; then
tail -n +2 "$file" | awk -v dir="$dir" 'BEGIN {OFS="\t"} {print dir, $0}' >> "$combined_read_length"
else
  echo "Warning: Read length file '$file' is empty or missing." >&2
fi
done




### Now run this locally
library(readr)
coverage <- read_tsv("C:\\Users\\lboukas\\combined_coverage_results.tsv")
dup_samples <- coverage$Sample[which(duplicated(coverage$Sample))]
coverage_no_duplicates <- coverage[-which(duplicated(coverage$Sample)), ]
mean(coverage_no_duplicates$`Average Coverage`)

read_length <- read_tsv("C:\\Users\\lboukas\\combined_read_length_results.tsv")
dup_samples <- read_length$Sample[which(duplicated(read_length$Sample))]
read_length_no_duplicates <- read_length[-which(duplicated(read_length$Sample)), ]
mean(read_length_no_duplicates$`Average Read Length`)



