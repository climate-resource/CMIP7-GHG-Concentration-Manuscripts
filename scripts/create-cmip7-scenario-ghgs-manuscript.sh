#!/bin/bash
# Handy trick (full details here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230):
# -e: exit immediately if any command fails
# -u: exit if you reference any unset variable
# -o: pipefail means that a non-zero exit code is returned if any command in the script fails
set -euo pipefail

# Create CMIP7 scenario GHGs manuscript
#
# Based on pytask.
#
# The first run downloads the parts of the scenario run's
# [Zenodo archive](https://zenodo.org/records/21501391) that we need
# and re-runs notebooks from that run in the run's own environment,
# which is installed with [pixi](https://pixi.sh) (so you need pixi installed).
# Installing that environment can be slow
# (for the historical run, the first notebook run after installing it took the better part of an hour),
# but it only happens once.
#
# A task re-runs when any of its dependencies or products changed
# (by content, not modification time) or when the module it is defined in changed.
# Helper modules (e.g. `src/local/scenario_ghg_forcing_for_cmip7/figures.py`)
# are only tracked where a task lists them as a dependency.
# The slow notebook re-runs are in their own module (`task_rerun_notebooks.py`)
# so that editing the other tasks doesn't re-run them.
#
# To see what would run, without running it:
#
# ```sh
# uv run pytask --dry-run
# # Every task's dependencies and products
# uv run pytask collect --nodes
# # The graph of the tasks (needs graphviz, e.g. `brew install graphviz`)
# uv run --with networkx --with pygraphviz pytask dag -o dag.pdf
# ```
#
# To re-run tasks even if they're up to date, use `-f` (force)
# with a selection (otherwise everything re-runs), e.g.
#
# ```sh
# # Re-run the notebooks
# uv run pytask -f -k rerun
# ```
#
# Tasks are selected like tests in pytest:
#
# - `-k` matches expressions against task IDs (module and function name), e.g.
#   `-k task_inverse_emissions_figure`, `-k "task_tex_inputs.py::task_summary_table"`,
#   `-k "figure and not table"`.
#   The whole graph is still collected,
#   so the tasks the selected ones depend on are selected too,
#   but they only run if they are out of date.
#   Tasks which depend on the selected ones are not run,
#   they catch up the next time you run `uv run pytask`
#   (if what they depend on actually changed).
# - `-m` matches expressions against markers (we don't use any yet).
# - Paths, e.g. `uv run pytask tasks/scenario_ghg_forcing_for_cmip7/task_manuscript_scenario_ghg_forcing_for_cmip7.py`,
#   limit which tasks are collected at all.
#   Tasks in other modules are then not part of the graph,
#   so what they produce is used as-is if it exists (even if it is out of date),
#   and the build fails if it doesn't.
#   Prefer `-k` unless you specifically want this.
#   (This is also why `pyproject.toml` doesn't set pytask's `paths`:
#   in pytask 0.6.0, `paths` set there overrides paths given on the command line.)
#
# pytask collects tasks from the directory it is run in,
# so it has to be run from the root of the repository.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

# Careful: a `-k` which matches no task deselects everything
# and pytask still reports success, so check the task name.
uv run pytask -k "task_manuscript_scenario_ghg_forcing_for_cmip7.py::task_compile_pdf"
