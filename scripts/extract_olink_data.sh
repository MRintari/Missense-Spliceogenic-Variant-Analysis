#!/bin/bash

# Set input/output files
vars="/slade/home/mr935/data/biobank_data/high_score_mask.tsv"
output="/slade/home/mr935/data/biobank_data/olink_high_score_mask.tsv"

# Add header to output file
if ! grep -q "^Gene\tVariant\tAlt_Freq\tBeta\tlog10$" "$output"; then
        echo -e "Gene\tVariant\tAlt_Freq\tBeta\tlog10" > "$output"
fi

# Loop through input file
tail -n +2 "$vars" | while IFS=$'\t' read -r gene variant; do
        echo "GENE: $gene VARIANT: $variant"
        # Find necessary file
        gene_file=$(ls /slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results/"${gene}"_rint_burden_wgs/single_variants/"${gene}"_*.regenie 2>/dev/null | head $
        echo "GENE FILE: $gene_file"

        if [[ -f "$gene_file" ]]; then
                # Search for variant in file
                match=$(awk -F' +' -v var="$variant" '($3 == var) {print}' "$gene_file")

                if [[ -n "$match" ]]; then
                        # Extract column values
                        Alt_Freq=$(echo "$match" | awk -F' +' '{print $6}')
                        Beta=$(echo "$match" | awk -F' +' '{print $9}')
                        Log10=$(echo "$match" | awk -F' +' '{print $12}')
                else
                        Alt_Freq="NA"
                        Beta="NA"
                        Log10="NA"
                fi
        else
                Alt_Freq="NA"
                Beta="NA"
                Log10="NA"
        fi

        # Append result to output file
        echo -e "${gene}\t${variant}\t${Alt_Freq}\t${Beta}\t${Log10}" >> "$output"
done
