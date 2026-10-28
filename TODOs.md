# TODOs

## scenario GHG forcing for CMIP7 manuscript

Compiled from `manuscripts/scenario-ghg-forcing-for-cmip7/*`
(every `TODO`, `ZNTODO`, `TODO-CHECK`, `zalue-check`, `[CITE ...]`, `X`/`TBD` placeholder and "Placeholder" comment)
plus the few TODOs in `src/local/scenario_ghg_forcing_for_cmip7/` which end up in the manuscript.
Line numbers are as of commit `c12567e`.
Paths are relative to `manuscripts/scenario-ghg-forcing-for-cmip7/` unless stated otherwise.

### Outstanding

- [x] make full TODO list
- [x] `methods.tex:336-354`: impact of the 2014/2015 jump in MAGICC's inputs
      (was "less than [TODO number, I think 0.03 degree C and 0.05 W / m2]").
      Done 2026-10-29: text gives the largest difference in MAGICC's CO2, CH4 and N2O over 2015-2500 and 1750-2014, with value checks.
    - Notebook `data/raw/scenario-ghg-forcing-for-cmip7/notebooks-modified/0021_run-magicc_cmip7-history-check.py`,
      runs `src/local/scenario_ghg_forcing_for_cmip7/magicc_jump.py`,
      comparison `src/local/scenario_ghg_forcing_for_cmip7/magicc_jump_summary.py`
      (separate so changing the comparison doesn't re-run MAGICC),
      tasks `tasks/scenario_ghg_forcing_for_cmip7/task_magicc_jump.py`
      (all 7 markers + VL-CF, 600 members, `published` and `cmip7-history` set-ups,
      outputs and `summary.json` in the bundle's `manuscript-outputs/magicc-jump`).
    - `published` reproduces the published output exactly (GSAT, total ERF, CO2, CH4, N2O, every experiment).
    - Checked: passing the scenario from 1750 is enough for MAGICC to use the CMIP7 historical emissions
      (its fossil CO2 emissions before 2015 equal the CMIP7 input), no new `*_EMIS.IN` files needed.
    - Largest median difference (cmip7-history minus published):
        - 1750-2014: GSAT 0.026 degC, ERF 0.063 W/m2 (both 1967; same in every experiment).
        - 2015-2100: GSAT 0.079 degC (H, 2100), ERF 0.11 W/m2 (M, 2025).
        - 2015-2500: GSAT 0.088 degC (H, 2297), ERF 0.13 W/m2 (HL, 2146).
    - GHG vs. aerosol ERF (both saved in every run and in `summary.json`):
        - In the 2020s, the ~0.1 W/m2 ERF difference is about half GHGs (~0.05, mostly N2O and CH4)
          and half aerosols (~0.04-0.05, less negative aerosol ERF with the CMIP7 history).
        - Later, the GHG difference grows (up to 0.18 W/m2, HL 2166; N2O alone reaches ~0.1 by 2500)
          while the aerosol difference changes sign (down to -0.1 W/m2), offsetting part of it.
          So the long-term GSAT difference is driven by GHGs.
        - Not investigated: why the aerosol ERF differs after 2015 when both set-ups use the same emissions
          (a guess: MAGICC scales aerosol forcing relative to emissions in a historical reference year).
    - Confirmed (ZN): the post-2015 aerosol difference comes from MAGICC scaling aerosol forcing
      relative to emissions in a historical reference year.
    - Concentrations (raw MAGICC medians, not harmonised, ZN's call), largest difference 2015-2500:
      CO2 6.8 ppm (HL, 2153; ~0.09 W/m2), CH4 99 ppb (HL, 2076; ~0.04 W/m2),
      N2O 18 ppb (VL, 2500; ~0.06 W/m2). To 2100: CO2 3.8 ppm, CH4 99 ppb, N2O 8.6 ppb.
      In the text (2015-2500 and 1750-2014) instead of the GSAT/ERF numbers (ZN, 2026-10-29).
    - The GSAT/ERF numbers (and their value checks) are for the ScenarioMIP paper instead.
    - Value checks: `magicc-jump-runs-reproduce-published`, `magicc-jump-n-runs`,
      `magicc-jump-max-<gsat|erf|erf-ghg|erf-aerosols|co2|ch4|n2o>-difference-<history|to-2100|to-2500>`.
      Original notes, how to do it (fiddly):
    - We need two runs.
      Step 1: a run which reproduces the output we already have (make sure we can reproduce the existing output).
      Step 2: alter the run so that it uses CMIP7 historical emissions and concentrations.
    - To get the concentrations in, we have to write all the `*CONC.IN` files.
    - To get the emissions in, ZN thinks we just need to write the SCEN file with input from 1750 onwards
      and don't need to write new `*_EMIS.IN` files, but we should check.
    - Easiest start: notebook `1200_run-magicc-concentration-driven` from `cmip7-scenariomip-ghg-concentrations`.
      Update it so that it runs all the CMIP7 ScenarioMIP scenarios
      (using the MAGICC emissions input files we already have from the output bundle here)
      in the two ways described above.
    - Running the notebooks will take a while, so this check needs to be cached
      (pytask probably does that for us already, so we likely don't need our own caching).
    - The difference between the two runs gives the maximum difference we could expect from this jump.
      If we had recalibrated MAGICC (as we wanted), the difference would be smaller,
      because we would have updated the model physics to stay close to the AR6 assessed ranges.
    - Searched for an existing notebook which does this (2026-10-28), none found:
      `emissions_harmonization_historical`, `cmip7-scenariomip-ghg-concentrations` and `repos-agcec/magicc`
      (all branches). `1200` only writes `*_CONC.IN` (CMIP7 emissions only from 2015).
      `repos-agcec/magicc/tests/regression/scripts/process_cmip7_scenariomip_markers.py` is our standard setup
      (CMIP6 concentrations to 2014), but its `processing.py` has helpers for writing MAGICC input files.

#### To ask co-authors

- [ ] Invite all the data producers to be co-authors (ZN to do, not a manuscript edit)
- [ ] `results.tex:160-162`: check with co-authors that future CFC-12-eq concentrations are slightly higher
      "because of updated handling of bank emissions" (and get the citations for it).
- [ ] `metadata.toml:32,65`: email for FB (currently `tbd@invalid.com`).
- [x] `metadata.toml:12`: author contributions are placeholders, confirm them.
- [x] `metadata.toml:1`: title is a placeholder (the comment points to `manuscript-planning.md`, which no longer exists).

#### Waiting on others

- [ ] `conclusion.tex:36-39`: C3S paragraph stays commented out until C3S is confirmed (also needs the C3S project name).
- [x] Scenario user guide needs to be published so we can cite/link it
      (`methods.tex:1099` "[TODO cite scenario user guide]", `conclusion.tex:6` "[TODO link]").
      Source is `notebooks/user-guide-scenarios.py`; the historical one is `nicholls_historical-user-guide_2026`.
      Published at https://zenodo.org/records/23719968, added as `nicholls_scenario-user-guide_2026`,
      cited in the methods and linked and cited in the conclusion (as in the historical manuscript).
- [x] Historical paper has no bib entry yet, so all the "[TODO cite historical paper]" items below wait on that
      (see Citations). Added `nicholls_historical-ghgs_2026` (unpublished, in preparation).

#### Text to write or rework

Introduction

- [x] `introduction.tex:1-2`: consider adding paragraphs about the need for standardisation, production of forcings, harmonisation etc.
- [x] `introduction.tex:4-19`: first paragraph "is a bit of a mess and probably trying to say too much", rework.

Methods

- [x] `methods.tex:112-113`: replace the bracketed "[TODO discussion point about this: ...]" (why the one-box route)
      with a forward ref, the discussion already covers it (`discussion.tex:49-54`).
- [x] `methods.tex:88`: decide whether the single concentration projection route needs a forward ref to the discussion
      (`discussion.tex:64-70` covers using a single source for ODSs). Added one.
- [x] `methods.tex:274-277`: we also run MAGICC v7.5.3 for every scenario.
      Fig. `methods-ch4` a, b already shows these runs, but the CO_2 Step 1 text doesn't say we run v7.5.3.
      Decide whether to say so there (ZNTODO says don't show the tiny CO_2 feedback difference unless asked).
      Now says we also run v7.5.3 (no wetland feedback) for comparison, not used in the output, pointing to CH_4 Step 1.
- [x] `methods.tex:288-290`: say that using CMIP6 concentrations to 2014 maintains consistency with the AR6 calibration,
      recalibrating on CMIP7 history would be better, but the effect is minor given how similar the histories are
      (partly covered at `methods.tex:293-294` already).
- [x] `methods.tex:353-356`: check use of "scenario" throughout
      ("scenario" only for the ScenarioMIP scenarios, "experiments"/"datasets" when VL-CF is included).
      Whole manuscript checked. Changed where VL-CF is included (results panel c and the groups, abstract,
      discussion, `methods.tex` "scenarion name", intro "scenario datasets").
      Left: ScenarioMIP/CMIP6 scenarios, scenarios as a concept, emissions "which underpin the scenarios".
      Intro now also says VL-CF's data starts in 2016.
- [x] VL-CF regression inconsistencies: state each in the methods text (and check `discussion.tex:72-87` covers them).
      Done, with value checks (`<key>-vl-cf-regression-slope-ratio`, computed by refitting up to 2015,
      which `vl-cf-regression-is-up-to-harmonisation-year` confirms reproduces VL-CF's CO_2 seasonality regression;
      `<gas>-<vl|vl-cf>-pc0-2300`). The discussion's fifth issue now says which are biggest
      and that VL-CF's CO_2 seasonality change follows history to 2023.
    - [x] CO_2 seasonality change PC vs. NPP (`methods.tex:425-435`, also the commented "that's a bug" note at `methods.tex:444-450`):
          VL-CF regresses over 1981-2015 (slope 1.78 yr / GtC) vs. 1981-2023 for the markers (1.28 yr / GtC),
          i.e. ~40% more seasonality change for the same NPP change.
          Also, the seasonality change PC is anchored in 2023 (not $y_h$) in every experiment, including VL-CF.
    - [x] CH_4 lat. gradient PC0 (`methods.tex:604-617`): VL-CF regression ~4% steeper, anchored in 2015. Small effect.
    - [x] N_2O lat. gradient PC0 (`methods.tex:686-697`): VL-CF regression ~40% steeper,
          PC0 ends near -0.14 (VL-CF) vs. -0.55 (VL) with the same emissions after 2100.
          One-EOF gases (`1030`) only differ by anchor, giving a constant offset (e.g. HFC-134a in 2300: VL-CF -1.81, VL 4.54).
- [x] `methods.tex:709-710`: check the discussion says that the historical one-EOF regressions are against RCMIP SSP2-4.5 emissions,
      while we apply them to the ScenarioMIP emissions (`discussion.tex:92` touches on this).
      It does now (sixth issue), also saying the re-derived regressions use CMIP7, not historical-dataset, emissions.
- [x] `methods.tex:62-68` and the commented-out VL-CF section (`methods.tex:990-1012`):
      VL-CF equals history up to 2021 and VL from 2022 for the non-MAGICC gases.
      The Fig. `methods-hfc125` caption says so, but nowhere says why. Add a sentence, then delete the commented-out section.
      Sentence added in the general approach (their emissions are VL's), the section was already deleted.
- [x] `methods.tex:325-326` (notes for the ScenarioMIP paper): check whether "scaling arguments for whatever else is left
      (based on M2020)" and "CO2 AFOLU from IAMs" need to be described here.
      No: the original run's code takes the emissions team's complete emissions as is (no AFOLU handling),
      so both belong to the ScenarioMIP emissions paper (cited for infilling). Notes left as they are.

Results

- [x] `results.tex:12`: rename the high/medium/low groups on the figures so they don't clash with the CMIP7 scenario names (H, M, L).
      Now "Higher peak", "Medium peak" and "Lower peak" (figure titles, results and emissions figure captions,
      "lower group" in the CH_4 results text).
- [x] `results.tex:13`: plot IGCC concentrations where available (shows which scenario is closest in 2023, 2024 and 2025).
      IGCC v6.4.0 from 2020 onwards, as circles in panel c of the CO_2, CH_4 and N_2O results figures,
      with an inset zooming in on 2020 to IGCC's last year (2025).
      Not shown for CFC-12-eq and HFC-134a-eq, IGCC defines them differently (captions say so).
      The ScenarioMIP scenarios are identical until 2024, so IGCC can't single out a scenario.
      Instead, one sentence per gas in the results compares the ScenarioMIP experiments with IGCC in 2025
      (CO_2 ~0.7 ppm below, CH_4 30-35 ppb above with ~11 ppb already there in 2022, N_2O ~0.2 ppb above),
      with value checks (`igcc-last-year`, `<gas>-scenariomip-minus-igcc-last-igcc-year`, `ch4-history-minus-igcc-2022`).

Discussion

- [x] Two "The sixth issue" paragraphs (`discussion.tex:89,96`), renumber.
- [x] (Optional) `discussion.tex:4-7`: possible addition, changes in scenarios compared to long-term records.

Appendices

- [x] `appendices.tex:27-28`: wherever we refer to the CH_4 regression (there and in `sssec:methods-ch4-project-seasonality-lat-gradient`),
      refer to the discussion of inconsistencies between inputs.
      The discussion's sixth issue now says the CH_4 slope is ~5.6 times the historical one and the fit (1750-2023)
      mixes the ice-core-optimised and emissions-regression years (value checks `ch4-lat-gradient-regression-*`).
      The CH_4 regression figure caption and the CH_4 Step 2 methods text point there.

Back matter

- [x] `code-and-data-availability.tex:1-2`: placeholder, "check this".
      Also needs the bundle/Zenodo reference for the VL-CF run,
      which was produced after the 1.1.0 bundle with later code (`methods.tex:995-998`).
      @Claude check: I think this bundle has VL-CF already
      Checked: 21501391 (2026-07-23) does have VL-CF, but it's a later version of the same Zenodo concept (18690744)
      as 19243739 (v1.1.0, 2026-03-27). The ScenarioMIP files' `doi` attribute is 19243739, VL-CF's is 21501391
      (21501391 also includes the ScenarioMIP files). The text now gives 19243739 for ScenarioMIP and 21501391 for VL-CF.
      The commented-out VL-CF section's ZNTODO (`methods.tex:995-998`) is resolved by this, delete it with that section.
- [x] `competing-interests.tex:1`: placeholder, confirm.
- [x] `conclusion.tex:1`: placeholder comment points to `manuscript-planning.md` (no longer exists), remove once happy.

#### Numbers to fill in and statements to check

Radiative effect numbers are calculated from our concentrations (not loaded from external ERF totals like IGCC).

- [x] `abstract.tex:10-20`: all the X/Y/S values (CO_2, CH_4, N_2O ranges, years and scenarios),
      the CMIP6 comparison ranges (CH_4 and N_2O are "TBD", CFC-12 "[higher/lower] by X"), plus value checks for all of them.
      Done (2022-2500 for both phases), with value checks for every number, year and experiment
      (`<gas>-<cmip7|cmip6>-<min|max>`, `-year`, `-is-<experiment>`, `-range-width-relative-to-cmip6`, `cfc12-*`).
      CMIP7's CFC-12 is above every CMIP6 SSP in every year, ~6% above the CMIP6 mean to 2100.
- [x] `results.tex:136`: N_2O radiative effect difference between SSP3-7.0 and H in 2500 ("[TODO quantify]").
      ~0.8 W/m2 (`n2o-ssp370-minus-h-2500-radiative-effect`).
- [x] `results.tex:34-37`: check "Post 2100, the [CO_2] datasets either settle between 300 and 350 ppm".
      Wasn't quite right: the settled experiments end at 312-363 ppm (VL-CF above 350), text now says 310-365 (`co2-settled-2500-range`).
- [x] `results.tex:97-99`: check "Post 2100, CH_4 concentrations stabilise between 900 and 1200 ppb".
      Holds (952-1065 ppb in 2500, `ch4-stabilised-2500-range`).
- [x] `methods.tex:505-512`: the notebook's check that 2022 monthly values match the historical dataset never runs
      (`.isel(time=slice(-6, 0, None))` selects nothing, presumably meant `slice(None, -6)`).
      Fix, re-check, then quantify "The differences are negligible [TODO quantify and value check]".
      Didn't touch the original notebook; a value check now compares the published 2022 files directly:
      at most 0.0014 W/m2 (CO_2, 0.11 ppm) over every gas, latitude and month
      (`max-2022-monthly-difference-from-history-radiative-effect`, also used for the same claim in the introduction).
- [x] `methods.tex:561-562`: check CH_4 VL vs. VL-CF "differ by ~70 ppb until 2040".
      Wasn't quite right: the difference peaks at ~69 ppb in 2041, rather than holding until 2040. Text now says
      "up to around 70 ppb (around 2040)" (`ch4-vl-minus-vl-cf-max`, `-max-year`).
- [x] `methods.tex:731-733`: list the gases, scenarios and years where the lat. gradient has to be scaled down further
      (currently the list is just the gases the code allows this for).
      Skipped (ZN's call, 2026-10-28). The ZNTODO comment is still in `methods.tex`.
- [x] `methods.tex:777-778`: value check for filling the final year (2500) with linear extrapolation.
      `single-projection-final-year-linearly-extrapolated` (2500 = 2 x 2499 - 2498 for every single projection gas).
- [x] `methods.tex:815-818`: check "small by 2100" with numbers, relative to each gas' maximum over historical and future.
      HCFC-22 8%, HCFC-141b 12%, HCFC-142b 33% of peak in 2100; together 0.006 W/m2. Text now gives both
      (`western-et-al-2024-2100-max-fraction-of-peak`, `western-et-al-2024-2100-radiative-effect`).
- [x] `methods.tex:939-940`: check the assumption that emissions are constant over the year in the one-box model.
      Holds: `0030` holds year y - 1's emissions constant over 12 monthly Euler steps (noted in a comment).
- [x] `methods.tex:943-945`: one-box gases don't guarantee a continuous growth rate at 2022, check how big the kinks are.
      The comment says "appendix figure of the join, see the figure plan at the top", but there's no figure plan any more,
      so decide whether we want that figure.
      Skipped (ZN's call, 2026-10-28). The ZNTODO comment is still in `methods.tex`.
- [x] `methods.tex:286`: confirm MAGICC uses the CMIP6 (RCMIP v4.0.0) concentrations up to 2014
      (neither the probabilistic distribution nor `0021_run-magicc.py` seem to override them).
      Skipped (ZN's call, 2026-10-28). The ZNTODO comment is still in `methods.tex`.
- [x] `methods.tex:1044-1047`: `GHG_RADIATIVE_EFFICIENCIES` docstring in the scenario code says Table 7.SM.6,
      the historical manuscript checked that 7.SM.7 is right. Check the scenario values match the historical ones.
      Values identical to the historical ones for every shared gas (scenario also has Halon-1201).
      Docstring updated to Table 7.SM.7 in `~/Documents/repos/cmip7-scenariomip-ghg-concentrations` (not committed).
- [x] `methods.tex:1114`: check whether M20 adds anything about the vertical dimension (as far as ZN knows, it points back to M17).
      M20 only provides surface mole fractions, no vertical treatment, so nothing to add (noted in a comment).

#### Value checks

All done 2026-10-29: every `TODO-CHECK`, `TODO zalue-check` and `zalue-check` is now a live `% value-check`
(105 checks pass). Where possible the checks derive the statement from the data
(e.g. the gas groups from which interim projection our output matches,
the convergence times from where our data meets MAGICC's median).
Text changed to match the data (please review):

- [x] `methods.tex` general approach, Step 2, and CO2 latitudinal gradient: the second PC is *not* simply
      held constant at its harmonisation-year value. It follows the historical dataset's PCs up to 2023
      (their last year), then is held at its harmonisation-year value.
      For VL-CF, it follows the historical PCs 2016-2023, then drops back to its 2015 value.
      Text now says this. Might deserve a mention in the discussion (artefact, particularly for VL-CF).
- [x] `methods.tex` MAGICC emissions: was "the year in which each experiment's emissions diverge from history
      (2023 for ScenarioMIP, 2016 for VL-CF)", which mixed the last year on history (2023) with the first year off it (2016).
      The ScenarioMIP experiments' emissions are identical up to 2023 and VL-CF's equal VL's up to 2015.
      Now "the last year in which each experiment's emissions follow history (2023 ..., 2015 for VL-CF)".
- [ ] `methods.tex` CO2 (and N2O) latitudinal gradient regression years: the original code compares an array index with years,
      so the regression starts two years before the first PC stops being constant
      (CO2: includes 1750, N2O: includes 1988 and 1989, CH4 unaffected).
      Text now describes this as "an artefact of how the code selects the years"; ZNTODO left for ZN to decide
      whether to call it a bug and/or mention it in the discussion.
- [x] `methods.tex` mean-preserving interpolation: monthly global-means in 2022 are the same across the ScenarioMIP
      experiments to within 1.4e-5 (relative), not exactly. Text now says "(to within 0.01\%)".
      (At native resolution, 2022 differs by up to 0.4% between experiments (CH2Cl2, January),
      from the latitudinal gradient; the text makes no claim about that.)
- [x] Everything else just converted (the statements held):
      lat band width, output first/last years (plus VL-CF's 2016), MAGICC ensemble size and switch year,
      convergence times (15/60/75 years), CFC-12-like and HFC-125-like gas lists,
      CFC-12-like agreement with history in 2022 (max 0.13%, CH3CCl3),
      lifetimes vs. WMO (2022) Table A-5 (transcribed from the annex PDF; SF6 uses the middle of the 850-1280 year range),
      HFC-125 largest one-box radiative effect in 2022 (0.0097 W/m2), number of output variables (46),
      source IDs, output time chunks, single-projection MAGICC differences.
- [x] CH_4 VL vs. VL-CF difference: done earlier (`ch4-vl-minus-vl-cf-max`).

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

- [x] `abstract.tex:12-13,18-19`: CH_4 and N_2O units are ppm, should be ppb.
- [x] `abstract.tex:14,20`: "CFC12" → "CFC-12".
- [x] `introduction.tex:104`: extensions run "2101 to the end of 2025", should be 2500.
- [ ] `introduction.tex:96`: "how add" → "how to add".
- [ ] `introduction.tex:79`: "once an agreed communication format has been agreed".
- [ ] `discussion.tex:66`: "Montreal Agreement" → "Montreal Protocol".
- [ ] `discussion.tex:102`: "The second, like also affected the historical dataset" (grammar).
- [ ] `methods.tex:360`: `[as done in][]\citep{meinshausen_ssp-ghgs_2020}` is malformed, should be `\citep[as done in][]{...}`.
- [x] `methods.tex:562`: `\~70ppb` renders as an accent, use `$\sim$70~ppb`.
- [x] `methods.tex:972`: unclosed parenthesis in the zero-seasonality gas list (closed when adding the historical citation).
- [ ] Spelling: "conistency" (`introduction.tex:15`), "counterfactural" (`introduction.tex:36`), "throughly" (`introduction.tex:78`),
      "experiemnt" (`output-requirements.tex:57`), "timesescale" (`methods.tex:654`), "scenarion" (`methods.tex:1084`).

#### Build code

- [x] `src/local/scenario_ghg_forcing_for_cmip7/scenarios.py:80,86`: official colour palettes, and a colour for VL-CF.
      Done 2026-10-29: official CMIP7 ScenarioMIP palette (from ZN). VL-CF isn't part of ScenarioMIP, so keeps its grey.
      The methods figures didn't depend on `scenarios.py` (so wouldn't pick up colour changes), fixed. All figures regenerated.

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
