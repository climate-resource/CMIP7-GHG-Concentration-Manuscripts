#!/bin/bash
# Handy trick (full details here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230):
# -e: exit immediately if any command fails
# -u: exit if you reference any unset variable
# -o: pipefail means that a non-zero exit code is returned if any command in the script fails
set -euo pipefail

# Create a PDF from a single source notebook file
#
# Options:
#
# -s: Source filepath
# -t: Title
# -d: Description
# --toc: Include a table of contents

source_filepath=""
title=""
description=""
use_toc=false
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

usage() {
    echo "usage: $0 -s source-file -d description -t title [--toc]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -s)
        if [[ $# -lt 2 ]]; then
            echo "Option -s requires an argument" >&2
            usage
            exit 1
        fi
        source_filepath="$2"
        shift 2
        ;;
    -t)
        if [[ $# -lt 2 ]]; then
            echo "Option -t requires an argument" >&2
            usage
            exit 1
        fi
        title="$2"
        shift 2
        ;;
    -d)
        if [[ $# -lt 2 ]]; then
            echo "Option -d requires an argument" >&2
            usage
            exit 1
        fi
        description="$2"
        shift 2
        ;;
    --toc)
        use_toc=true
        shift
        ;;
    *)
        echo "Unknown option: $1" >&2
        usage
        exit 1
        ;;
    esac
done

if [[ -z "${source_filepath}" || -z "${title}" || -z "${description}" ]]; then
    usage
    exit 1
fi

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
output_pdf="${build_dir}/${stem}.pdf"
output_pdfs_grouped_dir="${repo_root}/output-pdfs"
build_source="${source_md}"
build_args=(
    --builder latex
    --config "${build_dir}/config.yaml"
    --path-output "${build_dir}/"
)

if [[ "${use_toc}" == true ]]; then
    source_tex="${build_dir}/_build/latex/${stem}.tex"
    build_source="${build_dir}/"
    build_args+=(
        --toc "${build_dir}/toc.yaml"
    )
else
    source_tex="${build_dir}/_build/_page/${stem}/latex/${stem}.tex"
fi

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

if [[ "${use_toc}" == true ]]; then
    uv run python "${script_dir}/write-jupyter-book-toc-single-file.py" \
        --output "${build_dir}/toc.yaml" \
        --source-file "${source_md}"
fi

# Jupyter Book/Sphinx handles absolute config and TOC paths more reliably.
uv run jupyter-book build -v \
    "${build_source}" \
    "${build_args[@]}"

uv run python "${script_dir}/compile-jupyter-book-latex.py" \
    --source "${source_tex}" \
    --output "${output_pdf}" \
    --references-bib-file "${build_dir}/references.bib"

echo "Output file is in ${output_pdf}"

mkdir -p "${output_pdfs_grouped_dir}"
cp "${output_pdf}" "${output_pdfs_grouped_dir}"/"$(basename "${output_pdf}")"

echo "Output file is also in ${output_pdfs_grouped_dir}"
