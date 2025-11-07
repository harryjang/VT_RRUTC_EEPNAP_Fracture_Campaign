#!/bin/bash

# Define the destination folder
DEST_DIR="Collected_Output"

# Create the destination folder if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each MAIN_* folder
for folder in MAIN_*/Analyze/nodout_1/GetClusters/Output_files; do
    # Extract the MAIN_* folder name
    MAIN_FOLDER=$(echo "$folder" | cut -d'/' -f1)

    # Define the target directory within collected_output
    TARGET_DIR="$DEST_DIR/$MAIN_FOLDER"

    # Create the target directory if it doesn't exist
    mkdir -p "$TARGET_DIR"

    # Copy all files except "neighbors.csv" and "SavedSparseCluster.mat"
    find "$folder" -maxdepth 1 -type f ! -name "neighbors.csv" ! -name "saved_sparse_cluster.mat" -exec cp {} "$TARGET_DIR/" \;
done

echo "Files have been collected and organized in $DEST_DIR, excluding neighbors.csv and SavedSparseCluster.mat."
