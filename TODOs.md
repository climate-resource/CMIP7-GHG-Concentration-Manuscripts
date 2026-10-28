# TODOs

## scenario GHG forcing for CMIP7 manuscript

### Outstanding

- [ ] make full TODO list

#### To ask co-authors

- [ ] Invite all the data producers to be co-authors (ZN to do, not a manuscript edit)

## historical GHG forcing for CMIP7 manuscript

Compiled from `manuscripts/historical-ghg-forcing-for-cmip7/*`,
`src/local/`, the last build (`build/historical-ghg-forcing-for-cmip7/latex/main.log`, 2026-10-05)
the historical part of `manuscript-planning.md` (moved here, see the end of this file)
and the methods TODOs from `NOTES.md` (moved into section 1).
Line numbers are as of commit `bdf36a4`.
Paths are relative to `manuscripts/historical-ghg-forcing-for-cmip7/` unless stated otherwise.

### Outstanding

The section each item came from is in parentheses (see Done for the context).

#### To ask co-authors

- [ ] (Section 10) Invite all the data producers to be co-authors (ZN to do, not a manuscript edit)
- [ ] (Section 5) Ask co-authors for other (ideally independent) comparison datasets (`results.tex:10-11`).

#### Waiting on others

- [ ] (Section 6) The C3S paragraph at the end of `conclusion.tex` stays commented out until C3S is resolved
      (it also needs the C3S project name).

#### Needs editor or reviewer feedback

- [ ] (Section 3) "Float too large for page by 10pt": this is the CFC-12-like per-gas `table*` (`methods.tex:995-1020`).
      Being too tall is why it was pushed to the end of the document.
      A `\clearpage` at the end of `methods.tex` now flushes it there instead (and one around each appendix section
      keeps the appendix figures in their own appendix), but the warning remains.
      Leaving the warning for now: check with the editor whether they want this long table
      or for us to split it over multiple pages (ZNTODO next to the table in `methods.tex`).

#### Optional

Low priority, not manuscript-blocking (skip the `src/local` tooling items for now).

- [ ] (Section 2) CH_4 Law Dome smoothing figure (noise, windows, regression; `methods.tex:821-823`, "if requested by reviewers").
- [ ] (Section 9) `data_loading.py:140,298,333,337,353,461,494,498,520`: ESGF result checking, loading logic, hard-coded paths.
- [ ] (Section 9) `esgf/` TODOs: uniqueness constraints, process count, mapping organisation, node handling, URL sorting, chunk sizes, retry logic.
