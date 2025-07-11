# Install packages in command line
```
sudo apt install parallel
```
# Data aquisition and pre-processing

The Israeli patients' sequence data was downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906
Metabolomic and metadata were downloaded from the aticle supplementary file: https://www.nature.com/articles/s41467-024-48106-6

File processing was performed in command line as follows:

# Combine the data from all archive files
```
tar -xvf GSE199906_RAW.tar

zcat GSM*.txt.gz | awk '
   BEGIN { FS="\t"; OFS="," }      # Set input/output field separators
  {
    split($1, parts, "|");        # Split target_id by "|"
    print parts[1], parts[6], $2  # Print transcript_id, gene_name, TPM
  }' > extracted_data.csv
```

# Collect sample_id, gene_name, tpm, and transcript_id into one file
```
  echo "sample_id,gene_name,tpm,transcript_id" > final_output_combined.csv
  find . -name "*.txt.gz" -print0 | while IFS= read -r -d $'\0' file; do   zcat "$file" | awk -v file="$file" '
      BEGIN {FS="\t"; OFS=","}
      NR > 1 {
        sample_id = file
        sub(/^.*\//, "", sample_id)
        sub(/\.txt\.gz$/, "", sample_id)

        # Detect format
        if ($1 ~ /\|/) {
          # Pipe-delimited format (GSM691*)
          split($1, parts, "|")
          print sample_id, parts[6], $2, parts[1]
        } else if ($1 ~ /_/) {
          # Underscore-delimited format (GSM599*)
          split($1, parts, "_")
          print sample_id, parts[1], $2, parts[2]
        }
      }
    ' >> final_output_combined.csv; done

  awk -F, 'NR>1 {print $1}' final_output_combined.csv | sort | uniq | wc -l
  head -n 5 final_output_combined.csv && echo "..." && tail -n 5 final_output_combined.csv
```
# Add patients IDs from metadata
```
  awk 'BEGIN {FS=OFS=","}
       NR==FNR {a[$1]=$2 OFS $3 OFS $4; next}
       ($1 in a) {print $0, a[$1]}' Israel_metadata.csv final_output_combined.csv > fc_data.csv
```
