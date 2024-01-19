#!/bin/bash

# Remove all .nextflow.log* files in the current directory
find . -type f -name '.nextflow.log*' -delete

# Remove the content of the local work directory (modify the path if needed)
rm -rf work

echo "Cleanup complete."
