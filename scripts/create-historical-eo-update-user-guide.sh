#!/bin/bash
#
# Create the historical (EO update) user guide

bash scripts/create-pdf-from-single-notebook.sh \
    -s notebooks/user-guide-historical-EO-update.py \
    -t "Adding Earth Observation data to the generation of CMIP7 Greenhouse Gas (GHG) Concentration Historical Dataset: Data Description and User Guide" \
    -d "User guide for the CMIP7 greenhouse gas (GHG) concentration forcing historical dataset, updated with EO extensions." \
    --toc
