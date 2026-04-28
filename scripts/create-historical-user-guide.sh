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

stem="cmip-phase-comparison-historical-standalone"
build_dir="build/${stem}"
source_md="${build_dir}/${stem}.md"
source_tex="${build_dir}/_build/_page/${stem}/latex/${stem}.tex"
output_pdf="${build_dir}/${stem}.pdf"

mkdir -p "${build_dir}/"
uv run jupytext --to myst "notebooks/${stem}.py" --output "${source_md}"
check_output "Convert \`.py\` notebook file to myst with jupytext"

uv run python scripts/write-jupyter-book-config.py \
    --base jupyter-book/base-config.yaml \
    --output "${build_dir}/config.yaml" \
    --title "Comparison of CMIP6 and CMIP7 historical datasets" \
    --description "Comparison of the CMIP6 and CMIP7 historical greenhouse gas concentration forcing datasets." \
    --references-bib-file references/references.bib \
    --source-file "${source_md}"
check_output "Write config file"

uv run jupyter-book build -v \
    "${source_md}" \
    --builder latex \
    --config "${build_dir}/config.yaml" \
    --path-output "${build_dir}/"
check_output "Convert myst to latex with jupyter-book"

uv run python scripts/compile-jupyter-book-latex.py \
    --source "${source_tex}" \
    --output "${output_pdf}" \
    --references-bib-file "${build_dir}/references.bib"
check_output "Compile latex to pdf"

echo "Output file is in ${output_pdf}"
