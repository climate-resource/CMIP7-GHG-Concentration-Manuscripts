#!/bin/bash
#
# Create the historical user guide

function check_output() {

    last_command_status=$?

    if [[ $last_command_status -ne 0 ]]; then
        echo "Running \`$1\` failed"
        exit $last_command_status
    fi

    echo "Running \`$1\` passed"

}

mkdir -p build/cmip-phase-comparison-historical-standalone/
uv run jupytext --to myst notebooks/cmip-phase-comparison-historical-standalone.py --output build/cmip-phase-comparison-historical-standalone/cmip-phase-comparison-historical-standalone.md
check_output "Convert \`.py\` notebook file to myst with jupytext"

uv run python scripts/write-jupyter-book-config.py \
    --base jupyter-book/base-config.yaml \
    --output build/cmip-phase-comparison-historical-standalone/config.yaml \
    --title "Comparison of CMIP6 and CMIP7 historical datasets" \
    --description "Comparison of the CMIP6 and CMIP7 historical greenhouse gas concentration forcing datasets." \
    --references-bib-file references/references.bib \
    --source-file build/cmip-phase-comparison-historical-standalone/cmip-phase-comparison-historical-standalone.md
check_output "Write config file"

uv run jupyter-book build -v \
    build/cmip-phase-comparison-historical-standalone/cmip-phase-comparison-historical-standalone.md \
    --builder latex \
    --config build/cmip-phase-comparison-historical-standalone/config.yaml \
    --path-output build/cmip-phase-comparison-historical-standalone/
check_output "Convert myst to latex with jupyter-book"

uv run python scripts/compile-jupyter-book-latex.py \
    --source build/cmip-phase-comparison-historical-standalone/_build/_page/cmip-phase-comparison-historical-standalone/latex/projectnamenotset.tex \
    --output build/cmip-phase-comparison-historical-standalone/cmip-phase-comparison-historical-standalone.pdf \
    --references-bib-file build/cmip-phase-comparison-historical-standalone/references.bib
check_output "Compile latex to pdf"

echo "Output file is in build/cmip-phase-comparison-historical-standalone/cmip-phase-comparison-historical-standalone.pdf"
