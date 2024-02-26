#!/bin/bash

# Remove all .nextflow.log* files in the current directory
find . -type f -name '.nextflow.log*' -delete

# Remove the content of the local work directory (modify the path if needed)
rm -rf work

# Remove all null directories
rm -rf null

# Remove all traces
find . -type f -name 'trace.txt*' -delete

echo "Cleanup complete."
