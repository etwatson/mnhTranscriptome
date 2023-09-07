#!/bin/bash

# Script to annotate transcripts using GffCompare

# Input files
ASSEMBLED_TRANSCRIPTS=$1
REFERENCE_ANNOTATION=$2
OUTPUT_PREFIX=$3

# Check if files are provided
if [ -z "$ASSEMBLED_TRANSCRIPTS" ] || [ -z "$REFERENCE_ANNOTATION" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Usage: $0 <assembled_transcripts.gtf> <reference_annotation.gtf> <output_prefix>"
    exit 1
fi

# Run GffCompare
gffcompare -r $REFERENCE_ANNOTATION -G -o $OUTPUT_PREFIX $ASSEMBLED_TRANSCRIPTS

# Check if GffCompare ran successfully
if [ $? -eq 0 ]; then
    echo "GffCompare ran successfully. Check the output files with prefix: $OUTPUT_PREFIX"
else
    echo "Error running GffCompare."
    exit 1
fi

