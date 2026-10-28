# TODOs

## scenario GHG forcing for CMIP7 manuscript

Compiled from `manuscripts/scenario-ghg-forcing-for-cmip7/*`
(every `TODO`, `ZNTODO`, `TODO-CHECK`, `zalue-check`, `[CITE ...]`, `X`/`TBD` placeholder and "Placeholder" comment)
plus the few TODOs in `src/local/scenario_ghg_forcing_for_cmip7/` which end up in the manuscript.
Line numbers are as of commit `c12567e`.
Paths are relative to `manuscripts/scenario-ghg-forcing-for-cmip7/` unless stated otherwise.

### Outstanding

- [x] make full TODO list

#### To ask co-authors

- [ ] Invite all the data producers to be co-authors (ZN to do, not a manuscript edit)
- [ ] `results.tex:160-162`: check with co-authors that future CFC-12-eq concentrations are slightly higher
      "because of updated handling of bank emissions" (and get the citations for it).
- [ ] `metadata.toml:32,65`: email for FB (currently `tbd@invalid.com`).
- [x] `metadata.toml:12`: author contributions are placeholders, confirm them.
- [x] `metadata.toml:1`: title is a placeholder (the comment points to `manuscript-planning.md`, which no longer exists).

#### Waiting on others

- [ ] `conclusion.tex:36-39`: C3S paragraph stays commented out until C3S is confirmed (also needs the C3S project name).
- [ ] Scenario user guide needs to be published so we can cite/link it
      (`methods.tex:1099` "[TODO cite scenario user guide]", `conclusion.tex:6` "[TODO link]").
      Source is `notebooks/user-guide-scenarios.py`; the historical one is `nicholls_historical-user-guide_2026`.
- [x] Historical paper has no bib entry yet, so all the "[TODO cite historical paper]" items below wait on that
      (see Citations). Added `nicholls_historical-ghgs_2026` (unpublished, in preparation).

#### Text to write or rework

Introduction

- [ ] `introduction.tex:1-2`: consider adding paragraphs about the need for standardisation, production of forcings, harmonisation etc.
- [ ] `introduction.tex:4-19`: first paragraph "is a bit of a mess and probably trying to say too much", rework.

Methods

- [ ] `methods.tex:112-113`: replace the bracketed "[TODO discussion point about this: ...]" (why the one-box route)
      with a forward ref, the discussion already covers it (`discussion.tex:49-54`).
- [ ] `methods.tex:88`: decide whether the single concentration projection route needs a forward ref to the discussion
      (`discussion.tex:64-70` covers using a single source for ODSs).
- [ ] `methods.tex:274-277`: we also run MAGICC v7.5.3 for every scenario.
      Fig. `methods-ch4` a, b already shows these runs, but the CO_2 Step 1 text doesn't say we run v7.5.3.
      Decide whether to say so there (ZNTODO says don't show the tiny CO_2 feedback difference unless asked).
- [ ] `methods.tex:288-290`: say that using CMIP6 concentrations to 2014 maintains consistency with the AR6 calibration,
      recalibrating on CMIP7 history would be better, but the effect is minor given how similar the histories are
      (partly covered at `methods.tex:293-294` already).
- [ ] `methods.tex:353-356`: check use of "scenario" throughout
      ("scenario" only for the ScenarioMIP scenarios, "experiments"/"datasets" when VL-CF is included).
- [ ] VL-CF regression inconsistencies: state each in the methods text (and check `discussion.tex:72-87` covers them):
    - [ ] CO_2 seasonality change PC vs. NPP (`methods.tex:425-435`, also the commented "that's a bug" note at `methods.tex:444-450`):
          VL-CF regresses over 1981-2015 (slope 1.78 yr / GtC) vs. 1981-2023 for the markers (1.28 yr / GtC),
          i.e. ~40% more seasonality change for the same NPP change.
          Also, the seasonality change PC is anchored in 2023 (not $y_h$) in every experiment, including VL-CF.
    - [ ] CH_4 lat. gradient PC0 (`methods.tex:604-617`): VL-CF regression ~4% steeper, anchored in 2015. Small effect.
    - [ ] N_2O lat. gradient PC0 (`methods.tex:686-697`): VL-CF regression ~40% steeper,
          PC0 ends near -0.14 (VL-CF) vs. -0.55 (VL) with the same emissions after 2100.
          One-EOF gases (`1030`) only differ by anchor, giving a constant offset (e.g. HFC-134a in 2300: VL-CF -1.81, VL 4.54).
- [ ] `methods.tex:709-710`: check the discussion says that the historical one-EOF regressions are against RCMIP SSP2-4.5 emissions,
      while we apply them to the ScenarioMIP emissions (`discussion.tex:92` touches on this).
- [ ] `methods.tex:62-68` and the commented-out VL-CF section (`methods.tex:990-1012`):
      VL-CF equals history up to 2021 and VL from 2022 for the non-MAGICC gases.
      The Fig. `methods-hfc125` caption says so, but nowhere says why. Add a sentence, then delete the commented-out section.
- [ ] `methods.tex:325-326` (notes for the ScenarioMIP paper): check whether "scaling arguments for whatever else is left
      (based on M2020)" and "CO2 AFOLU from IAMs" need to be described here.

Results

- [ ] `results.tex:12`: rename the high/medium/low groups on the figures so they don't clash with the CMIP7 scenario names (H, M, L).
- [ ] `results.tex:13`: plot IGCC concentrations where available (shows which scenario is closest in 2023, 2024 and 2025).

Discussion

- [ ] Two "The sixth issue" paragraphs (`discussion.tex:89,96`), renumber.
- [ ] (Optional) `discussion.tex:4-7`: possible addition, changes in scenarios compared to long-term records.

Appendices

- [ ] `appendices.tex:27-28`: wherever we refer to the CH_4 regression (there and in `sssec:methods-ch4-project-seasonality-lat-gradient`),
      refer to the discussion of inconsistencies between inputs.

Back matter

- [ ] `code-and-data-availability.tex:1-2`: placeholder, "check this".
      Also needs the bundle/Zenodo reference for the VL-CF run,
      which was produced after the 1.1.0 bundle with later code (`methods.tex:995-998`).
- [ ] `competing-interests.tex:1`: placeholder, confirm.
- [ ] `conclusion.tex:1`: placeholder comment points to `manuscript-planning.md` (no longer exists), remove once happy.

#### Numbers to fill in and statements to check

Radiative effect numbers are calculated from our concentrations (not loaded from external ERF totals like IGCC).

- [ ] `abstract.tex:10-20`: all the X/Y/S values (CO_2, CH_4, N_2O ranges, years and scenarios),
      the CMIP6 comparison ranges (CH_4 and N_2O are "TBD", CFC-12 "[higher/lower] by X"), plus value checks for all of them.
- [ ] `results.tex:136`: N_2O radiative effect difference between SSP3-7.0 and H in 2500 ("[TODO quantify]").
- [ ] `results.tex:34-37`: check "Post 2100, the [CO_2] datasets either settle between 300 and 350 ppm".
- [ ] `results.tex:97-99`: check "Post 2100, CH_4 concentrations stabilise between 900 and 1200 ppb".
- [ ] `methods.tex:340-344`: impact of the 2014/2015 jump in MAGICC's inputs
      ("less than [TODO number, I think 0.03 degree C and 0.05 W / m2]").
      Notebooks doing these runs exist somewhere, maybe in emissions-harmonisation-historical; add a (cached) value check.
- [ ] `methods.tex:505-512`: the notebook's check that 2022 monthly values match the historical dataset never runs
      (`.isel(time=slice(-6, 0, None))` selects nothing, presumably meant `slice(None, -6)`).
      Fix, re-check, then quantify "The differences are negligible [TODO quantify and value check]".
- [ ] `methods.tex:561-562`: check CH_4 VL vs. VL-CF "differ by ~70 ppb until 2040".
- [ ] `methods.tex:731-733`: list the gases, scenarios and years where the lat. gradient has to be scaled down further
      (currently the list is just the gases the code allows this for).
- [ ] `methods.tex:777-778`: value check for filling the final year (2500) with linear extrapolation.
- [ ] `methods.tex:815-818`: check "small by 2100" with numbers, relative to each gas' maximum over historical and future.
- [ ] `methods.tex:939-940`: check the assumption that emissions are constant over the year in the one-box model.
- [ ] `methods.tex:943-945`: one-box gases don't guarantee a continuous growth rate at 2022, check how big the kinks are.
      The comment says "appendix figure of the join, see the figure plan at the top", but there's no figure plan any more,
      so decide whether we want that figure.
- [ ] `methods.tex:286`: confirm MAGICC uses the CMIP6 (RCMIP v4.0.0) concentrations up to 2014
      (neither the probabilistic distribution nor `0021_run-magicc.py` seem to override them).
- [ ] `methods.tex:1044-1047`: `GHG_RADIATIVE_EFFICIENCIES` docstring in the scenario code says Table 7.SM.6,
      the historical manuscript checked that 7.SM.7 is right. Check the scenario values match the historical ones.
- [ ] `methods.tex:1114`: check whether M20 adds anything about the vertical dimension (as far as ZN knows, it points back to M17).

#### Value checks

None of the scenario manuscript's checks are live yet: they're all written as `TODO-CHECK`, `TODO zalue-check` or `zalue-check`
(the misspelling stops the checker treating them as real checks).
Each needs implementing in `src/local/scenario_ghg_forcing_for_cmip7/value_checks.py`, then converting to `% value-check: {...}`.

- [ ] Already implemented, just convert the comment:
      `single-projection-magicc-max-erf-difference-sum-l` (`methods.tex:855`),
      `single-projection-magicc-max-erf-difference-sum-range` (`methods.tex:857`).
- [ ] Tag and bounds given, needs implementing:
    - [ ] `native-grid-lat-band-width` (`methods.tex:46,1066`)
    - [ ] `output-first-year` (`output-requirements.tex:20`, `methods.tex:528`), `output-last-year` (`output-requirements.tex:21`)
    - [ ] `magicc-n-ensemble-members` (`methods.tex:268`)
    - [ ] `magicc-switch-from-conc-to-emissions-year` (`methods.tex:279`)
    - [ ] `co2-harmonisation-convergence-time` (`methods.tex:386`), `ch4-...` (`methods.tex:550`), `n2o-...` (`methods.tex:646`)
    - [ ] `scenarios-identical-december-2022` (`methods.tex:501`)
    - [ ] `cfc12-like-wmo-historical-agreement-2022` (`methods.tex:782`)
    - [ ] `n-output-variables` (`methods.tex:1061`)
- [ ] Tag only, needs bounds/spec as well as implementing:
    - [ ] `co2-seasonality-change-npp-regression-years` (`methods.tex:424`)
    - [ ] `co2-lat-gradient-regression-years` (`methods.tex:478`)
    - [ ] `n2o-lat-gradient-regression-years` (`methods.tex:699`)
    - [ ] `cfc12-like-gases` (`methods.tex:752`), `hfc125-like-gases` (`methods.tex:920`)
    - [ ] `one-box-lifetimes` (`methods.tex:810`, against WMO 2022 Table A-5)
    - [ ] `hfc125-largest-erf-one-box-group-2022` (`methods.tex:917`)
    - [ ] `source-ids` (`methods.tex:1080`), `output-time-chunks` (`methods.tex:1091`)
- [ ] No tag yet:
    - [ ] last year PC2 held constant is 2022 (`methods.tex:131`)
    - [ ] years when emissions diverge from history (2023 ScenarioMIP, 2016 VL-CF) (`methods.tex:330`)
    - [ ] CH_4 VL vs. VL-CF difference (`methods.tex:561`, see above)

#### Citations and cross-references

Historical paper (now `nicholls_historical-ghgs_2026`):
`methods.tex:5,9,40,71,149,412,464,473,489,591,601,704,972,974,1017,1040,1063,1113`
(and commented `371`),
`introduction.tex:55,67,81,84,96,108` (`81` is written "[CITE historical paper]"),
`output-requirements.tex:4`, `results.tex:20,159`, `discussion.tex:91,102`, `conclusion.tex:15`, `appendices.tex:31`
and the generated C_8F_18 caption (`src/local/scenario_ghg_forcing_for_cmip7/appendix_figures.py:221`).

- [x] Add the historical paper bib entry and replace all of the above.
      `\citep` where the citation is parenthetical, `\citet` after "noted in"/"discussed in"/"refer to",
      `\citealp` inside existing parentheses. Also updated `nicholls_scenario-ghgs_2026` to the scenario `metadata.toml` title and authors.
- [ ] `introduction.tex:96`: also needs the section number in the historical paper (vertical dimension).
      It was 3.11 (`ssec:methods-vertical`) in the last historical build (2026-10-21),
      but the historical methods have changed since, so re-check before hard-coding it.

Others (bib entry exists for those marked ✓):

- [ ] ScenarioMIP / van Vuuren ✓ (`van-vuuren_scenariomip_2025`): `methods.tex:17,52`, `conclusion.tex:3`,
      `output-requirements.tex:13` ("[cite van Vuuren]"), plus "other citations" at `introduction.tex:21`.
- [ ] Underlying papers for each scenario: `introduction.tex:26`.
- [ ] Riahi: `introduction.tex:33`.
- [ ] SSP paper: `introduction.tex:10`.
- [ ] Paris Agreement: `introduction.tex:30,32`.
- [ ] PolMIP (Annika): `introduction.tex:35` (check citation), `methods.tex:17,18` (PolMIP stuff and dataset),
      `methods.tex:53`, `output-requirements.tex:13` ("[cite Annika?]").
- [ ] Kikstra et al. (harmonisation): `methods.tex:18,307,329` (plus "list of other sources"), `methods.tex:474` (plus underlying datasets),
      `discussion.tex:22` (with Gidden).
- [ ] Nicholls et al. (ScenarioMIP infilling): `methods.tex:309,774`.
- [ ] Sandstad et al. (extensions): `methods.tex:328,825,826`.
- [ ] MAGICC: Meinshausen et al. 2009 (`methods.tex:26,263`),
      Meinshausen et al. 2011 RCP concentrations (`methods.tex:28,263`),
      Meinshausen et al. 2020 ✓ (`methods.tex:263`),
      RCMIP phase 2 (`methods.tex:27,272`),
      AR6-calibrated configuration "[TODO citations]" (`methods.tex:96`),
      feedbacks "[TODO cite papers again?]" (`methods.tex:101`),
      AR6 WG1 Ch. 7 ✓ for probabilistic output (`methods.tex:348`).
- [ ] Sloughter et al. (wetland CH_4 feedback): `methods.tex:273,545`.
- [ ] gradient-aware-harmonisation package: `methods.tex:373`.
- [ ] CDIAC conversion table or equivalent: `methods.tex:808`.
- [ ] Decide whether to cite the WMO data file for CMIP7 (`MixingRatiosCMIP7_20250210.xlsx`) separately: `methods.tex:771`.
- [ ] Bank emissions (CFC-12-eq): `results.tex:162` (see To ask co-authors).

Cross-references:

- [ ] `discussion.tex:36`: ref the section showing VL vs. VL-CF harmonisation differences (CH_4 Step 1, `methods.tex:558-566`, and N_2O Step 1).
- [ ] `methods.tex:853`: forward ref to the note on approximate radiative effect (`ssec:results-radiative-forcing-comparisons-note`).
- [ ] `methods.tex:1101`: repeat the guidance pages link.

#### Typos and small fixes

Not marked as TODOs, spotted while compiling this list.

- [ ] `abstract.tex:12-13,18-19`: CH_4 and N_2O units are ppm, should be ppb.
- [ ] `abstract.tex:14,20`: "CFC12" → "CFC-12".
- [ ] `introduction.tex:104`: extensions run "2101 to the end of 2025", should be 2500.
- [ ] `introduction.tex:96`: "how add" → "how to add".
- [ ] `introduction.tex:79`: "once an agreed communication format has been agreed".
- [ ] `discussion.tex:66`: "Montreal Agreement" → "Montreal Protocol".
- [ ] `discussion.tex:102`: "The second, like also affected the historical dataset" (grammar).
- [ ] `methods.tex:360`: `[as done in][]\citep{meinshausen_ssp-ghgs_2020}` is malformed, should be `\citep[as done in][]{...}`.
- [ ] `methods.tex:562`: `\~70ppb` renders as an accent, use `$\sim$70~ppb`.
- [x] `methods.tex:972`: unclosed parenthesis in the zero-seasonality gas list (closed when adding the historical citation).
- [ ] Spelling: "conistency" (`introduction.tex:15`), "counterfactural" (`introduction.tex:36`), "throughly" (`introduction.tex:78`),
      "experiemnt" (`output-requirements.tex:57`), "timesescale" (`methods.tex:654`), "scenarion" (`methods.tex:1084`).

#### Build code

- [ ] `src/local/scenario_ghg_forcing_for_cmip7/scenarios.py:80,86`: official colour palettes, and a colour for VL-CF.

#### Needs editor or reviewer feedback

- [ ] `introduction.tex:52-53`: the rest of the introduction "isn't really introduction", could become its own section.
- [ ] `output-requirements.tex:74`: maybe mention other limits.
- [ ] `results.tex:8`: deliberately minimal, add things based on feedback.

#### Not manuscript edits

- [ ] Produce corrected versions of VL-CF and upload them to ESGF
      (`discussion.tex:73-75`, possibly in the DestinE context; `discussion.tex:84` commits to "endeavour" to do this).

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
