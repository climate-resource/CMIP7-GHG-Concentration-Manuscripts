#!/bin/bash
# Handy trick (full details here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230):
# -e: exit immediately if any command fails
# -u: exit if you reference any unset variable
# -o: pipefail means that a non-zero exit code is returned if any command in the script fails
set -euo pipefail

# Create a PDF with a table of contents from a single source notebook file
#
# Options:
#
# -s: Source filepath
# -t: Title
# -d: Description

source_filepath=""
title=""
description=""
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

while getopts "s:t:d:" OPTION; do
    case $OPTION in
    s) source_filepath="${OPTARG}" ;;
    t) title="${OPTARG}" ;;
    d) description="${OPTARG}" ;;
    *)
        echo "usage: $0 -s source-file -d description -t title" >&2
        exit 1
        ;;
    esac
done

if [[ "${source_filepath}" = /* ]]; then
    source_filepath_abs="${source_filepath}"
else
    source_filepath_abs="${repo_root}/${source_filepath}"
fi

if [[ ! -f "${source_filepath_abs}" ]]; then
    echo "Source file does not exist: ${source_filepath_abs}" >&2
    exit 1
fi

cd "${repo_root}" || exit 1

source_filename=$(basename "${source_filepath_abs}")
stem="${source_filename%.*}"

build_dir="${repo_root}/build/${stem}"
source_md="${build_dir}/${stem}.md"
source_tex="${build_dir}/_build/latex/${stem}.tex"
output_pdf="${build_dir}/${stem}.pdf"
output_pdfs_grouped_dir="${repo_root}/output-pdfs"

mkdir -p "${build_dir}/"
uv run jupytext --to myst "${source_filepath_abs}" --output "${source_md}"

uv run python "${script_dir}/write-jupyter-book-config.py" \
    --base "${repo_root}/jupyter-book/base-config.yaml" \
    --output "${build_dir}/config.yaml" \
    --title "${title}" \
    --description "${description}" \
    --references-bib-file "${repo_root}/references/references.bib" \
    --source-file "${source_md}" \
    --source-file-raw "${source_filepath_abs}"

uv run python "${script_dir}/write-jupyter-book-toc-single-file.py" \
    --output "${build_dir}/toc.yaml" \
    --source-file "${source_md}"

# Jupyter Book/Sphinx handles absolute config and TOC paths more reliably.
uv run jupyter-book build -v \
    "${build_dir}/" \
    --builder latex \
    --config "${build_dir}/config.yaml" \
    --toc "${build_dir}/toc.yaml" \
    --path-output "${build_dir}/"

uv run python "${script_dir}/compile-jupyter-book-latex.py" \
    --source "${source_tex}" \
    --output "${output_pdf}" \
    --references-bib-file "${build_dir}/references.bib"

echo "Output file is in ${output_pdf}"

mkdir -p "${output_pdfs_grouped_dir}"
cp "${output_pdf}" "${output_pdfs_grouped_dir}"/"$(basename "${output_pdf}")"

echo "Output file is also in ${output_pdfs_grouped_dir}"
