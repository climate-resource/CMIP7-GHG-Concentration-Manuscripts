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

zenodo_bundle_dir="${repo_root}/data/raw/cmip-ghg-concentration-generation/v1.0.0"
# Where we store the original zenodo bundle once downloaded.
# Careful: if you change this, make sure to gitignore the new destination too.
original_run_notebooks_dir="${repo_root}/../CMIP-GHG-Concentration-Generation/output-bundles/v1.0.0/notebooks-executed"
# Where the original run notebooks are.
# We didn't include these in the zenodo archive, stupidly.
# If you need these but don't have them, you have to ask Zeb or someone who does have them.

abstract_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/abstract.tex"
introduction_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/introduction.tex"
output_requirements_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/output-requirements.tex"

methods_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/methods.tex"
co2_methods_figure_file="${repo_root}/figures/historical-ghg-forcing-for-cmip7/co2_methods.pdf"
ch4_methods_figure_file="${repo_root}/figures/historical-ghg-forcing-for-cmip7/ch4_methods.pdf"
n2o_methods_figure_file="${repo_root}/figures/historical-ghg-forcing-for-cmip7/n2o_methods.pdf"
# Every gas processed like SF6. These all share one figure,
# so this is the list of gases to draw it for, not a list of figures.
# It has to match local.historical_ghg_forcing_for_cmip7.SF6_LIKE_GASES;
# a gas which is not in that tuple is rejected before anything is drawn.
sf6_like_gases=(
    c2f6 c3f8 ccl4 cf4
    cfc11 cfc113 cfc114 cfc115 cfc12
    ch2cl2 ch3br ch3ccl3 ch3cl chcl3
    halon1211 halon1301 halon2402
    hcfc141b hcfc142b hcfc22
    hfc125 hfc134a hfc143a hfc152a hfc227ea
    hfc23 hfc236fa hfc245fa hfc32 hfc365mfc hfc4310mee
    nf3 sf6 so2f2
)

sf6_like_methods_figure_files=()
for sf6_like_gas in "${sf6_like_gases[@]}"; do
    sf6_like_methods_figure_files+=(
        "${sf6_like_gas}=${repo_root}/figures/historical-ghg-forcing-for-cmip7/${sf6_like_gas}_methods.pdf"
    )
done

# Every gas processed like C4F10. As with the SF6-like gases,
# these all share one figure, so this is the list of gases to draw it for.
# It has to match local.historical_ghg_forcing_for_cmip7.C4F10_LIKE_GASES.
c4f10_like_gases=(
    c4f10 c5f12 c6f14 c7f16 cc4f8
)

c4f10_like_methods_figure_files=()
for c4f10_like_gas in "${c4f10_like_gases[@]}"; do
    c4f10_like_methods_figure_files+=(
        "${c4f10_like_gas}=${repo_root}/figures/historical-ghg-forcing-for-cmip7/${c4f10_like_gas}_methods.pdf"
    )
done

# C8F18 is in a group of its own, so it is a plain file rather than a list.
c8f18_methods_figure_file="${repo_root}/figures/historical-ghg-forcing-for-cmip7/c8f18_methods.pdf"

# methods_subfile="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/methods-detail.tex"
results_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/results.tex"
code_and_data_availability_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/code-and-data-availability.tex"
author_contribution_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/author-contribution.tex"
competing_interests_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/competing-interests.tex"
acknowledgments_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/acknowledgements.tex"

conclusion_file="${repo_root}/manuscripts/historical-ghg-forcing-for-cmip7/conclusion.tex"
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
# `${a[@]+"${a[@]}"}` rather than `"${a[@]}"` throughout:
# bash 3.2, which is what macOS ships, treats an empty array as unset
# and `set -u` then kills the script, so an empty list of gases
# would fail here rather than simply drawing no SF6-like figures.
sf6_like_args=()
for sf6_like_methods_figure_file in ${sf6_like_methods_figure_files[@]+"${sf6_like_methods_figure_files[@]}"}; do
    sf6_like_args+=(--sf6-like-methods-figure-file "${sf6_like_methods_figure_file}")
done

c4f10_like_args=()
for c4f10_like_methods_figure_file in ${c4f10_like_methods_figure_files[@]+"${c4f10_like_methods_figure_files[@]}"}; do
    c4f10_like_args+=(--c4f10-like-methods-figure-file "${c4f10_like_methods_figure_file}")
done

uv run python "${script_dir}/historical-ghg-forcing-for-cmip7/generate-tex-inputs.py" \
    --co2-methods-figure-file "${co2_methods_figure_file}" \
    --ch4-methods-figure-file "${ch4_methods_figure_file}" \
    --n2o-methods-figure-file "${n2o_methods_figure_file}" \
    --c8f18-methods-figure-file "${c8f18_methods_figure_file}" \
    ${sf6_like_args[@]+"${sf6_like_args[@]}"} \
    ${c4f10_like_args[@]+"${c4f10_like_args[@]}"} \
    --bundle-dir "${zenodo_bundle_dir}" \
    --original-run-notebooks-dir "${original_run_notebooks_dir}"
# --force-rerun \

# create pdf or dump out to a single text file (that can then be dumped onto google docs, maybe easiest to do this with AI)
#   - caching here would be cool based on changes to the input hashes or content excluding comments, but likely overkill for many steps
uv run python "${script_dir}/compile-gmd-template-based-latex.py" \
    --abstract "${abstract_file}" \
    --introduction "${introduction_file}" \
    --section "${output_requirements_file}" \
    --section "${methods_file}" \
    --figure-file "<n2o-methods-figure>=${n2o_methods_figure_file}" \
    --section "${results_file}" \
    --conclusion "${conclusion_file}" \
    --code-and-data-availability "${code_and_data_availability_file}" \
    --author-contribution "${author_contribution_file}" \
    --competing-interests "${competing_interests_file}" \
    --acknowledgements "${acknowledgments_file}" \
    --replacements "${replacements_file}" \
    --metadata "${latex_metadata_file}" \
    --references-bib-file "${references_bib}" \
    --copernicus-template-dir "${copernicus_latex_template_dir}" \
    --clean-copernicus-template-filename "${clean_copernicus_template_filename}" \
    --build-dir "${build_dir}" \
    --output "${output_pdf}"
# --extra "${methods_subfile}" \

echo "Output file is in ${output_pdf}"
