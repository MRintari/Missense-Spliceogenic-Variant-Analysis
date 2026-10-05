# Define input directory containing VEP files and output file path
VEP_DIR="/slade/home/mr935/data/biobank_data/mane" 
OUTPUT_FILE="/slade/home/mr935/data/biobank_data/LoF_ukb_variants.vep.txt"

# Clear/create the output file
> "$OUTPUT_FILE"

# Flag to ensure header row is only written once
HEADER_COPIED=0

# Loop through all VEP files in the directory
for file in "$VEP_DIR"/*.vep.txt; do
    # Skip if no matching files found
    [ -e "$file" ] || continue
    
    echo "Processing: $file"

    awk -v header_done="$HEADER_COPIED" '
        BEGIN { FS="\t"; OFS="\t" }
        
        # Skip VEP meta-header lines
        /^##/ { next }
        
        # Process the column header line
        /^#/ {
            if (header_done == 0) {
                print $0
            }
            next
        }
        
        # Filter data rows: 14th column == HIGH and 87th column == HC
        $14 == "HIGH" && $87 == "HC" {
            print $0
        }
    ' "$file" >> "$OUTPUT_FILE"

    # Once the first file processes the header, set flag to 1
    HEADER_COPIED=1
done

echo "Done! Filtered variants saved to: $OUTPUT_FILE"
