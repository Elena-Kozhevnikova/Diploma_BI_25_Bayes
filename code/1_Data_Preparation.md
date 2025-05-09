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

  cat extracted_data.csv --head
  cat extracted_data.csv

  echo "sample_id,gene_name,tpm,transcript_id" > final_output.csv
  find . -name "*.txt.gz" -print0 | parallel -0 -j $(nproc) "
    zcat {} | awk -v file=\"{}\" '
    BEGIN {FS=\"\t\"; OFS=\",\"}
      {
       # Use full filename (without path) as sample_id
       sample_id = file
        sub(/^.*\//, \"\", sample_id)  # Remove directory path
        sub(/\.txt\.gz$/, \"\", sample_id)  # Remove extension
        split(\$1, parts, \"|\")
        if (length(parts) >= 6 && parts[6] != \"\") {
          print sample_id, parts[6], \$2, parts[1]
        }
      }
    '
  " >> final_output.csv


  echo "Total .txt.gz files found:"
  find . -name "*.txt.gz" | wc -l

  echo "sample_id,gene_name,tpm,transcript_id" > final_output1.csv
  find . -name "*.txt.gz" -print0 | while IFS= read -r -d $'\0' file; do   echo "Processing: $file" >&2;   zcat "$file" | awk -v file="$file" '
      BEGIN {FS="\t"; OFS=","; processed=0}
      {
        # Extract base filename without path/extension
        sample_id = file
        sub(/^.*\//, "", sample_id)
        sub(/\.txt\.gz$/, "", sample_id)

        # Count processed records
        split($1, parts, "|")
        if (length(parts) >= 6 && parts[6] != "") {
          print sample_id, parts[6], $2, parts[1]
          processed++
        }
      }
      END {
        if (processed == 0) {
          print "DEBUG: No valid records in " file > "/dev/stderr"
        }
      }
    ' >> final_output1.csv; qz
```
# Change sanple names
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

  awk -F, 'NR>1 {print $1}' final_output.csv | sort | uniq | wc -l
  awk -F, 'NR>1 {print $1}' final_output_combined.csv | sort | uniq | wc -l
  head -n 5 final_output_combined.csv && echo "..." && tail -n 5 final_output_combined.csv

  awk 'BEGIN {FS=OFS="\t"}
       NR==FNR {a[$1]=$2 OFS $3 OFS $4; next}
       $1 in a {print $0, a[$1]}' metadata.tsv expression.tsv > merged_data.tsv


  cat merged_data.tsv | head
output
  A025    DDX11L1 0.0     ENST00000456328.2       Israel_CD       30_39   male
  A025    DDX11L1 0.0     ENST00000450305.2       Israel_CD       30_39   male
  A025    WASH7P  1.85587 ENST00000488147.1       Israel_CD       30_39   male
  A025    MIR6859-1       0.0     ENST00000619216.1       Israel_CD       30_39   male
  A025    RP11-34P13.3    0.0     ENST00000473358.1       Israel_CD       30_39   male
  A025    RP11-34P13.3    0.0     ENST00000469289.1       Israel_CD       30_39   male
  A025    MIR1302-2       0.0     ENST00000607096.1       Israel_CD       30_39   male
  A025    FAM138A 0.0     ENST00000417324.1       Israel_CD       30_39   male
  A025    FAM138A 0.0209104       ENST00000461467.1       Israel_CD       30_39   male

  echo -e "sample_id\tgene_name\ttpm\ttranscript_id\tPatient_group\tAge\tGender" | cat - merged_data.tsv > headed_data.tsv
```
# Perform DGE in command line:
  awk 'BEGIN {
      FS=OFS="\t";
      print "gene_name\tmean_CD\tmean_Control\tlog2FC\tCD_samples\tControl_samples";
  }
  NR==1 {next} # Skip header
  {
      genes[$2][$5] += $3; # Sum TPM by gene and group
      counts[$2][$5]++;    # Count samples per group
  }
  END {
      for (gene in genes) {
          cd_mean = (counts[gene]["Israel_CD"] > 0) ? genes[gene]["Israel_CD"]/counts[gene]["Israel_CD"] : 0;
          ctrl_mean = (counts[gene]["Israel_control"] > 0) ? genes[gene]["Israel_control"]/counts[gene]["Israel_control"] : 0;

          # Avoid division by zero in fold change calculation
          log2fc = (ctrl_mean > 0) ? log(cd_mean/ctrl_mean)/log(2) : "NA";

          print gene, cd_mean, ctrl_mean, log2fc,
                counts[gene]["Israel_CD"]+0, counts[gene]["Israel_control"]+0;
      }
  }' headed_data.tsv > dge_results.tsv
```
