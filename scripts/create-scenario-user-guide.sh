#!/bin/bash
#
# Create the scenario user guide

bash scripts/create-pdf-from-single-notebook.sh \
    -s notebooks/user-guide-scenarios.py \
    -t "CMIP7 Greenhouse Gas (GHG) Concentration Forcing Scenarios Dataset" \
    -d "User guide for the CMIP7 greenhouse gas (GHG) concentration forcing scenarios dataset." \
    --toc
