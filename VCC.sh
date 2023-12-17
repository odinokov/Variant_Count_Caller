#!/usr/bin/bash

# Function to display usage information
usage() {
    cat <<EOF
Usage: $0 -c <CPU> -b <BAM file> -r <Reference genome file> -e <BED file>
Options:
  -c    Number of CPUs to use
  -b    Path to the BAM file
  -r    Path to the reference genome file
  -e    Path to the BED file
  -h    Display this help message
EOF
}

# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist."
        exit 1
    fi
}

# Function to check if required tools are installed
check_dependencies() {
    for cmd in samtools bcftools awk; do
        if ! command -v $cmd &> /dev/null; then
            echo "Error: $cmd is not installed."
            exit 1
        fi
    done
}

# Default CPU calculation
cpu=$(($(grep -c ^processor /proc/cpuinfo) - 1))
if [ $cpu -gt 8 ]; then
    cpu=8
fi

# Parse command-line options
while getopts ":c:b:r:e:h" opt; do
    case $opt in
        c) cpu=$OPTARG ;;
        b) bam=$OPTARG ;;
        r) ref=$OPTARG ;;
        e) bed=$OPTARG ;;
        h) usage; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$cpu" ] || [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$bed" ]; then
    echo "Error: All arguments are required."
    usage
    exit 1
fi

# Check if required files exist
check_file_exists "$bam"
check_file_exists "$ref"
check_file_exists "$bed"

# Check for required dependencies
check_dependencies

# Create a temporary file
temp_bam=$(mktemp)

# Step 1: Filter BAM file
if samtools view -@ "${cpu}" -h -F 3084 -q 30 -L "${bed}" "${bam}" \
   | awk '{if($7=="=") print $0}' \
   | cat <(samtools view -H "${bam}") - \
   | samtools sort -@ "${cpu}" -o "${temp_bam}" - \
   && samtools index -@ "${cpu}" "${temp_bam}"; then

# Step 2: Run bcftools mpileup
   bcftools mpileup -a "INFO/AD" -d 100000 --threads "${cpu}" -f "${ref}" -R "${bed}" -Ov "${temp_bam}" \
    | awk -f <(cat << 'EOF'
        BEGIN {
            FS = "\t";
            OFS = "\t";
            print "Chrom", "Pos", "ID", "Ref", "A_count", "C_count", "G_count", "T_count";
        }

        /^#/ { next; }

        {
            # Initialize counts for A, C, G, T to zero
            counts["A"] = counts["C"] = counts["G"] = counts["T"] = 0;

            match($8, /AD=([^;]+)/, arr);
            split($4 "," $5, alleles, /,|<[^>]*>/);

            for (i = 1; i <= split(arr[1], ad_values, ","); i++) {
                allele = alleles[i];
                if (allele ~ /^[ACGT]$/) counts[allele] = ad_values[i];
            }

            print $1, $2, $3, $4, counts["A"], counts["C"], counts["G"], counts["T"];
        }
EOF
    )
    # Optional: Clean up temporary files
    rm "${temp_bam}" "${temp_bam}.bai"
else
    echo "Error during BAM file processing."
    exit 1
fi