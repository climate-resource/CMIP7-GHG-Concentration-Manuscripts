#!/bin/bash
# Handy trick (full details here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230):
# -e: exit immediately if any command fails
# -u: exit if you reference any unset variable
# -o: pipefail means that a non-zero exit code is returned if any command in the script fails
set -euo pipefail

# run python stuff to generate inputs
#   - caching in the python (just have user config to set the caching for each step with basic decorators)
# create pdf
#   - caching here would be cool based on changes to the input hashes or content excluding comments, but maybe overkill?
