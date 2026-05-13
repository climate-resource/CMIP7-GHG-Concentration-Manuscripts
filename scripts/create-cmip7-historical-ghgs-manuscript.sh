#!/bin/bash
# Handy trick (full details here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230):
# -e: exit immediately if any command fails
# -u: exit if you reference any unset variable
# -o: pipefail means that a non-zero exit code is returned if any command in the script fails
set -euo pipefail

# Create CMIP7 historical GHGs manuscript
#
# Our take on a build pipeline for latex-based manuscripts.
# Let's see what we can reuse.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

cd "${repo_root}"

abstract_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/abstract.tex"
latex_metadata_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/metadata.toml"
replacements_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/replacements.yaml"
references_bib="${repo_root}/references/references.bib"

clean_copernicus_template_filename="template_clean.tex"
copernicus_latex_template_dir="${repo_root}/copernicus-latex-package"

build_dir="${repo_root}/build/historical-ghg-forcing-for-cmip7"

output_pdf_dir="${repo_root}/compiled-manuscripts"
output_pdf="${output_pdf_dir}/historical-ghg-forcing-for-cmip7.pdf"

mkdir -p "${output_pdf_dir}/"

# run python stuff to generate inputs
#   - caching in the python (just have user config to set the caching for each step with basic decorators)

# create pdf or dump out to a single text file (that can then be dumped onto google docs, maybe easiest to do this with AI)
#   - caching here would be cool based on changes to the input hashes or content excluding comments, but likely overkill for many steps
uv run python "${script_dir}/compile-gmd-template-based-latex.py" \
    --abstract "${abstract_file}" \
    --replacements "${replacements_file}" \
    --metadata "${latex_metadata_file}" \
    --references-bib-file "${references_bib}" \
    --copernicus-template-dir "${copernicus_latex_template_dir}" \
    --clean-copernicus-template-filename "${clean_copernicus_template_filename}" \
    --build-dir "${build_dir}" \
    --output "${output_pdf}"

echo "Output file is in ${output_pdf}"
