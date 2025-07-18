#!/bin/bash

# Stop on first error
set -e

# --- Configuration ---
ROOT_DIR=$(pwd)
DATA_DIR="$ROOT_DIR/exampleData"
GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
COMPRESSED_FILE_NAME="hg38.fa.gz"
FINAL_FILE_NAME="hg38.fa"
FULL_COMPRESSED_PATH="$DATA_DIR/$COMPRESSED_FILE_NAME"
FULL_FINAL_PATH="$DATA_DIR/$FINAL_FILE_NAME"

# --- Check for wget ---
if ! command -v wget &> /dev/null
then
    echo "Error: wget is not installed. Please install it to continue."
    exit 1
fi

# --- Create Data Directory ---
echo "--- Preparing Data Directory ---"
mkdir -p "$DATA_DIR"
echo "Directory '$DATA_DIR' is ready."
echo ""


echo "--- Downloading Human Reference Genome (hg38) ---"
echo "This is a large file (~1 GB) and may take a while."

if [ -f "$FULL_FINAL_PATH" ]; then
    echo "File '$FINAL_FILE_NAME' already exists. Skipping download."
else
    if [ -f "$FULL_COMPRESSED_PATH" ]; then
        echo "Compressed file '$COMPRESSED_FILE_NAME' already exists. Skipping download."
    else
        wget -O "$FULL_COMPRESSED_PATH" "$GENOME_URL"
        echo "Download complete."
    fi

    echo ""
    echo "--- Decompressing Genome File ---"
    gunzip -k "$FULL_COMPRESSED_PATH"
    echo "Decompression complete."
fi

echo ""
echo "----------------------------------------------------"
echo "Setup finished successfully!"
echo "The genome file is located at: $FULL_FINAL_PATH"
echo "----------------------------------------------------"