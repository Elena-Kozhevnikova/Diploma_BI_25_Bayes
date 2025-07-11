# Prepare data for further analysis using ```awk``` in command line

## Data aquisition and pre-processing

The Israeli patients' sequence data was downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906
Metabolomic and metadata were downloaded from the aticle supplementary file: https://www.nature.com/articles/s41467-024-48106-6

## Install packages in command line
```
sudo apt install parallel
```

## Combine the data from all archive files
```
tar -xvf GSE199906_RAW.tar

echo "sample_id,transcript_id,gene_name,TPM" > extracted_data.csv

find . -name "GSM*.txt.gz" -print0 | parallel -0 '
    filename={};
    basename=$(basename "$filename" .txt.gz);
    sample_id="$basename";

    zcat "$filename" | awk -v sid="$sample_id" '\''BEGIN { FS="\t"; OFS="," }
    NR>1 {
        split($1, parts, "_");
        gene_name=parts[1];
        transcript_id=parts[2];
        TPM=$2;
        print sid, transcript_id, gene_name, TPM
    }'\'' >> temp_output_"$basename"
'

cat temp_output_* >> extracted_data.csv
rm temp_output_*
```

## Keep patients IDs in ANNN format
```
awk -F',' 'NR==1 {print; next} {
    match($1, /A[0-9]{3}/);
    if (RSTART != 0) {
        $1 = substr($1, RSTART, RLENGTH);
    }
    print
}' extracted_data.csv > fc_data.csv
```
