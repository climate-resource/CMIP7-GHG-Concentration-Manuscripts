# Editing the manuscripts in Word

The historical paper and scenario paper have separate MM branches and worktrees.
Open the main repository folder for `historical-manuscript-draft-1_MM`, and the
sibling `CMIP7-GHG-Concentration-Manuscripts-scenario_MM` folder for
`scenario-manuscript-draft-1_MM`. Each has its own environment and build outputs.

## Build and export

From the relevant worktree, synchronize its environment and build its paper:

```sh
uv sync --locked
# Historical worktree
bash scripts/create-cmip7-historical-ghgs-manuscript.sh
# Scenario worktree
bash scripts/create-cmip7-scenario-ghgs-manuscript.sh
```

The export needs Pandoc, Poppler (`pdftoppm`), and Python with `python-docx`.
With Codex, use its bundled artifact Python runtime. Otherwise, create a separate
tool environment, for example `uv run --with python-docx python ...`.

```sh
python scripts/manuscript-word-review.py export --manuscript historical
# Or, in the scenario worktree:
python scripts/manuscript-word-review.py export --manuscript scenario
```

The editable copies are `compiled-manuscripts/<manuscript>_MM.docx`.
Track Changes is enabled. Use comments for questions and suggestions; edit prose
normally. Track Changes is helpful but not required for detecting text edits.
Equations, tables, citations, headings and paragraph moves need individual review.
TODOs hidden in LaTeX comments are included in a manuscript TODO section, with
their original filenames and line numbers. Existing TODOs in visible prose remain
in place.

Each export has a unique identifier in its Word properties. Its untouched original,
source snapshots, Git commit, paragraph map and conversion intermediates are kept
in `review-baselines/<manuscript>/<export-id>/`. Keep that folder. Neither it nor
the editable review document is committed to Git. Copy the entire baseline folder
if transferring the review workflow to another computer.

## Return edits to LaTeX

Save the edited document with a new filename and return it to Codex. Do not
overwrite its baseline. Ask Codex to reconcile that review with the corresponding
MM branch. The comparison command prepares the evidence:

```sh
python scripts/manuscript-word-review.py compare path/to/reviewed.docx \
  --output build/word-review/report.json
```

It also writes a readable `.md` report with text edits, tracked revisions,
comment anchors and candidate source sections. It detects changes to source files
since export and marks overlapping files conservatively as conflicts. It refuses
documents with a mismatched export identifier. If needed, pass
`--baseline review-baselines/<manuscript>/<export-id>/manifest.json` explicitly.

This is an **agent-assisted merge**, not an automatic DOCX-to-LaTeX replacement.
Codex compares the original source snapshot, current source and reviewed text;
then edits the original files under `manuscripts/`. Citation keys, labels, macros,
equations and figure-generation code must be retained. Comments become entries
in the review report and `% MM review ...` annotations near the relevant LaTeX.
Ambiguous scientific suggestions remain unresolved. Generated tables or captions
must be changed at their generating source, not patched only in build output.

After transfer, inspect `git diff`, rebuild with scientific checks enabled, and
export a fresh Word copy. Rebuilding checks numeric statements against the data;
an intentional scientific change may require reviewing its corresponding check.

## Verification

```sh
python scripts/test-manuscript-word-review.py --manuscript historical --compile
python scripts/test-manuscript-word-review.py --manuscript scenario --compile
```

These tests use disposable copies: tracked and ordinary edits, an inserted
paragraph, a comment, source conflict detection, and preservation of citations
and equation markup. The final delivery also requires compilation of the edited
test manuscript and visual inspection of the Word render. A comparison report
alone does not establish scientific correctness or conversion fidelity.
