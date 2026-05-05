#!/bin/bash
#
# Create a PDF from a single source notebook file
#
# Options:
#
# -s: Source filepath
# -t: Title
# -d: Description

source_filepath=""
title=""
description=""

while getopts "s:t:d:" OPTION; do
    case $OPTION in
    s) source_filepath="${OPTARG}" ;;
    t) title="${OPTARG}" ;;
    d) description="${OPTARG}" ;;
    *)
        echo "usage: $0 [-s source-file] [-d description] [-t title]" >&2
        exit 1
        ;;
    esac
done

if [[ -z $source_filepath ]]; then
    echo "Please provide a source file path with the \`-s\` flag"
    exit 1
fi

if [[ -z $title ]]; then
    echo "Please provide a title with the \`-t\` flag"
    exit 1
fi

if [[ -z $description ]]; then
    echo "Please provide a description with the \`-d\` flag"
    exit 1
fi

function check_output() {

    last_command_status=$?

    if [[ $last_command_status -ne 0 ]]; then
        echo "XXXXXXXXXXXXXXXXXXXXX"
        echo "Running \`$1\` failed"
        echo "XXXXXXXXXXXXXXXXXXXXX"
        exit $last_command_status
    fi

    echo "====================="
    echo "Running \`$1\` passed"
    echo "====================="

}

source_filename=$(basename "${source_filepath}")
stem="${source_filename%.*}"

build_dir="build/${stem}"
source_md="${build_dir}/${stem}.md"
source_tex="${build_dir}/_build/_page/${stem}/latex/${stem}.tex"
output_pdf="${build_dir}/${stem}.pdf"
output_pdfs_grouped_dir="output-pdfs"

mkdir -p "${build_dir}/"
uv run jupytext --to myst "${source_filepath}" --output "${source_md}"
check_output "Convert \`.py\` notebook file to myst with jupytext"

uv run python scripts/write-jupyter-book-config.py \
    --base jupyter-book/base-config.yaml \
    --output "${build_dir}/config.yaml" \
    --title "${title}" \
    --description "${description}" \
    --references-bib-file references/references.bib \
    --source-file "${source_md}" \
    --source-file-raw "${source_filepath}"
check_output "Write config file"

pwd=$(pwd)
# So stupid, --config has to be absolute so sphinx doesn't explode
uv run jupyter-book build -v \
    "${source_md}" \
    --builder latex \
    --config "${pwd}/${build_dir}/config.yaml" \
    --path-output "${build_dir}/"
check_output "Convert myst to latex with jupyter-book"

uv run python scripts/compile-jupyter-book-latex.py \
    --source "${source_tex}" \
    --output "${output_pdf}" \
    --references-bib-file "${build_dir}/references.bib"
check_output "Compile latex to pdf"

echo "Output file is in ${output_pdf}"

mkdir -p "${output_pdfs_grouped_dir}"
cp "${output_pdf}" "${output_pdfs_grouped_dir}"/"$(basename "${output_pdf}")"
check_output "Copy output pdf to ${output_pdfs_grouped_dir}"

echo "Output file is also in ${output_pdfs_grouped_dir}"
