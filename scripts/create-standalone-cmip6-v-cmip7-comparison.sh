#!/bin/bash
#
# Create the standalone CMIP6 vs. CMIP7 comparison

bash scripts/create-pdf-from-single-notebook.sh \
    -s notebooks/cmip-phase-comparison-historical-standalone.py \
    -t "Comparison of CMIP6 and CMIP7 historical datasets" \
    -d "Comparison of the CMIP6 and CMIP7 historical greenhouse gas concentration forcing datasets."
