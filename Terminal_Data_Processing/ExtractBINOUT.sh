#!/bin/bash

# Set path to the executable
EXTRACTOR="/home/harryjang/Simulations/Campaign_Template_Fracture/Midcourse_Data_Processing/ExtractBINOUT.l2a"

# Define the base directory where MAIN_* folders are located
BASE_DIR=$(pwd)

# Create the directories for collected matsum and glstat files
COLLECTED_MATSUM_DIR="$BASE_DIR/Collected_MATSUM"
COLLECTED_GLSTAT_DIR="$BASE_DIR/Collected_GLSTAT"

mkdir -p "$COLLECTED_MATSUM_DIR"
mkdir -p "$COLLECTED_GLSTAT_DIR"

# Iterate over all MAIN_* folders
for MAIN_FOLDER in MAIN_*; do
    # Ensure it's a directory
    if [[ -d "$MAIN_FOLDER" ]]; then
        echo "Processing folder: $MAIN_FOLDER"
        
        # Define the path to the binout0000 file
        BINOUT_FILE="$MAIN_FOLDER/binout0000"

        # Check if the binout0000 file exists
        if [[ -f "$BINOUT_FILE" ]]; then
            # Change into the MAIN_FOLDER directory before running the extractor
            cd "$MAIN_FOLDER"
            
            # Run the extraction command
            "$EXTRACTOR" binout0000
            
            # Return to the base directory
            cd "$BASE_DIR"

            # Extract the wildcard portion of the folder name
            WILDCARD_SUFFIX="${MAIN_FOLDER#MAIN_}"

            # Check if the matsum file exists and copy it
            if [[ -f "$MAIN_FOLDER/matsum" ]]; then
                NEW_MATSUM_FILE="$COLLECTED_MATSUM_DIR/matsum_${WILDCARD_SUFFIX}"
                cp "$MAIN_FOLDER/matsum" "$NEW_MATSUM_FILE"
                echo "Copied matsum file to: $NEW_MATSUM_FILE"
            else
                echo "Warning: matsum file not found in $MAIN_FOLDER"
            fi

            # Check if the glstat file exists and copy it
            if [[ -f "$MAIN_FOLDER/glstat" ]]; then
                NEW_GLSTAT_FILE="$COLLECTED_GLSTAT_DIR/glstat_${WILDCARD_SUFFIX}"
                cp "$MAIN_FOLDER/glstat" "$NEW_GLSTAT_FILE"
                echo "Copied glstat file to: $NEW_GLSTAT_FILE"
            else
                echo "Warning: glstat file not found in $MAIN_FOLDER"
            fi
        else
            echo "Warning: binout0000 file not found in $MAIN_FOLDER"
        fi
    fi
done

echo "Script completed."
