#!/bin/bash
#
# Create the historical user guide

bash scripts/create-pdf-from-single-notebook.sh \
    -s notebooks/user-guide-historical.py \
    -t "CMIP7 Greenhouse Gas (GHG) Concentration Forcing Historical Dataset" \
    -d "User guide for the CMIP7 greenhouse gas (GHG) concentration forcing historical dataset." \
    --toc
