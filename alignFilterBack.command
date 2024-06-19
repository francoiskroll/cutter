#!/bin/bash

shopt -s nullglob

checkPath() {
    local file_path="$1"

    # Check if the file exists
    if [ -e "$file_path" ]; then
        :
    else
        echo "Error. File does not exist: $file_path"
        exit 1
    fi
}

### from user, get path to config file
CFG="$1"
checkPath "$CFG"

### from user, get folder with fastq files
FAQ="$2"
checkPath "$FAQ"

### from user, get folder with fasta references
FAS="$3"
checkPath "$FAS"

################################################################################

# loop through rows of the config file        
while IFS=',' read -r WELL REF
do
    echo ""
    echo "WELL = $WELL"
    echo "REF = $REF"

    ### in fastq folder, find files that start with WELL
    find /path/to/folder -type f -name 'D01*'
fi

done < "$CFG"