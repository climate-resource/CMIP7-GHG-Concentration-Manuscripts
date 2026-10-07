# TODOs: historical GHG forcing for CMIP7 manuscript

Compiled from `manuscripts/historical-ghg-forcing-for-cmip7/*`,
`src/local/`, the last build (`build/historical-ghg-forcing-for-cmip7/latex/main.log`, 2026-10-05)
the historical part of `manuscript-planning.md` (moved here, see the end of this file)
and the methods TODOs from `NOTES.md` (moved into section 1).
Line numbers are as of commit `bdf36a4`.
Paths are relative to `manuscripts/historical-ghg-forcing-for-cmip7/` unless stated otherwise.

## 1. Methods: step labelling consistency

The general approach (`methods.tex:26-154`) defines six steps.
Everything else should use exactly these names:

| Step | Name in general approach |
| ---- | ------------------------ |
| 1 | Collect observations |
| 2 | Bin and interpolate observations |
| 3 | Principal component analysis |
| 4 | Extend |
| 5 | Prepare the components to be combined |
| 6 | Generate native resolution |

- [x] Fix all headers throughout methods (from `NOTES.md`). Done:
    - [x] Step 5 subsubsections renamed from "Process components" to "Prepare the components to be combined",
          matching the general approach and summary list
    - [x] every per-gas step subsubsection labelled `sssec:methods-<gas>-<step>`,
          with `<gas>` one of `n2o`, `co2`, `ch4`, `cfc12-like`, `c4f10-like`
          and `<step>` one of `collect-observations`, `bin-and-interpolate`, `pca`, `extend`,
          `prepare-components`, `generate-native-resolution`; all `\ref`s updated
    - [x] per-gas subsection titles unified as "Generating concentrations for <gas(es)>"
    - [x] C_8F_18 now says what it takes from CMIP6 (global- and hemispheric-, annual-means),
          what it assumes (zero seasonality, linear lat. gradient) and that Steps 5 and 6 are as for N_2O
          (checked against the original run's `1505_c8f18-like_create-pieces-for-gridding` notebook)
    - [x] "see step 3" → "see Step 3", "section \ref" → "Section \ref"
    - [x] stale commented sub-headings (`% \subsubsubsubsection{...}`) turned into bold run-in headings
          (`\textbf{Zonal mean.}` etc.; `\paragraph` is `\subsubsection` in the Copernicus class, so isn't an option).
          Rule: one run-in per component ("Global-, annual-mean.", "Seasonality.", "Latitudinal gradient.")
          wherever a step treats the components separately, in the same order for every gas.
          The order differs between steps because of dependencies
          (in Step 4, the global-mean is calculated from the already-extended latitudinal gradient).
          Where a component doesn't apply to a gas, its run-in is left in commented out,
          with a `% Doesn't apply: <reason>` line above it.
- [x] Sub-steps clash with top-level steps. Done:
    - [x] N_2O Step 3: "The first/second/third/fourth step" → "sub-step"
    - [x] CO_2 Step 3: "The fourth step" → "The last sub-step";
          "residuals from the previous step" → "the residuals from which it was derived"
    - [x] N_2O/CO_2 Step 5: "The last step is to construct..." → "In this step, we construct...";
          CO_2 Step 5's final part → "Finally, we calculate..." (matching N_2O)
    - [x] General Step 5: "The first step is to to put..." → "All components must first be put..."
          (also fixes the "to to" typo)
    - [x] CFC-12-like Step 4 global-mean: "two steps" / "the previous step" → "two sub-steps" / "the first sub-step"
    - Left as is: "this and subsequent steps" (CO_2 Step 3), "processing steps" (CH_4 Step 4)
      and "the next step is to then extend" (C_4F_10-like Step 4), which do refer to top-level steps.
- [x] Per-gas summary list vs the per-gas sections. Done:
    - [x] CFC-12-like: added the Step 2 delta (no data needed in the polar boxes, constant extrapolation into them,
          over multiple boxes for some gases, pointing to the per-gas table; removed the ZNTODO in general Step 2),
          the Step 3 delta (obs-derived seasonality and lat. gradient for every gas, unlike M17)
          and the Step 4 details (regression against RCMIP emissions, quartic gap fill to pre-industrial values).
          Checked against `1301_sf6-like_interpolate-observational-network`:
          every CFC-12-like gas gets the one-box extension, `gases_long_poleward_extension` the multi-box one.
    - [x] C_4F_10-like: Step 3 now says the global-mean is derived there and both are extrapolated from 2018 to 2022;
          Step 4 now says the PC, as well as concentrations, is zero before the Droste data.
    - [x] C_8F_18: now matches the section (CMIP6 global- and hemispheric-means, zero seasonality, linear lat. gradient).
          Kept 2015-2022: the original notebook takes SSP2-4.5 for 2015-2023,
          but the published output ends in 2022 (`..._gm_1750-2022.nc`).
- [x] General Step 4 (`methods.tex:94-105`) says the lat. gradient PCs are extended either by emissions regression
      or by linear/constant extrapolation. It doesn't cover the CH_4 ice-core optimisation (this is already noted in an exception),
      the CO_2 seasonality change PC (temperature-CO_2 composite, this is already noted in an exception)
      or the C_4F_10-like zero assumption (this is a special case of constant i.e. already covered). Generalise or add "except where noted".
- [x] General Step 2 says "The raw observations are binned" (`methods.tex:43`)
      but Step 1 says we use monthly aggregated data rather than raw data. Align terminology.
      (These are two different things, the difference is deliberate)
- [x] Fix the intro section of each per-gas sub-section so they're consistent
      and explain the differences/extra information compared to the base case (from `NOTES.md`). Done:
    - [x] CO_2/CH_4/CFC-12/C_4F_10 call N_2O the "base case" but point to
          `ssec:methods-general-approach` (`methods.tex:496,723,901,1136`); point to `ssec:methods-n2o` or reword.
          Now all "the base case (N_2O, Section \ref{ssec:methods-n2o})".
    - [x] Use the same opening pattern everywhere
          (currently "We begin with", "Next we consider", "Last is") and list which steps differ, matching the summary list.
          Pattern: "<gas> follows the base case (N_2O, ...), except in Step X (<name>), ...",
          then one sentence on why, with each reason tagged by its step.
          N_2O: "follows the general approach without modification, so we treat it as the base case";
          C_8F_18: "does not follow the base case".
          CFC12 intro now also mentions Step 2 (was missing) and the quartic gap fill in Step 4;
          CO_2 intro now also mentions the merged Law Dome - Mauna Loa record (Step 4).
          C_4F_10 intro no longer mentions the zero-before-Droste assumption,
          as it's a special case of constant extension (i.e. Step 4 follows the general approach, matching the summary list).
- [x] C_4F_10-like Step 1 (`methods.tex:1150-1155`) talks about not doing interpolation or PCA (that belongs in Steps 2/3)
      Move that text into Steps 2/3.
- [x] Clean up each per-gas sub-section so it doesn't repeat the base case more than needed (from `NOTES.md`).
      Candidates: CO_2 Step 3 re-describing zonal/global-mean and lat. gradient (`methods.tex:562-582`),
      CH_4 Step 4 smoothing detail vs "very similar to CO_2" (`methods.tex:843-861`),
      Step 5/6 sections that only say "same as N_2O". Done:
    - Kept every per-gas step heading (and run-in heading) so the six-step structure and labels stay intact,
      but cut the text under them.
    - Whole step identical: "This step is the same as for N_2O (Section \ref{sssec:methods-n2o-<step>})."
      Component identical: "As for N_2O (Section ...)." under the run-in.
      Step 2 for CO_2/CH_4: "Binning and interpolation are the same as for N_2O (...)", then the obs. network years.
      All of these now point to the N_2O step, not `ssec:methods-general-approach`
      (CO_2/CH_4/CFC12 Step 2 and CO_2 Step 6 used to).
    - CO_2 Step 3: dropped the restated "As for N_2O, we then take our spatially complete dataset..." opener
      and the re-explanation of retaining two EOFs (a ZNTODO to check the variance explained for CO_2 is left in).
      CO_2 Step 5: dropped the restated "In this step, we construct..." opener.
    - CH_4 Step 4: lat. gradient regression now "As for CO_2"; the Law Dome global-mean paragraph
      now says it's used the same way as Menking et al. for CO_2 and only gives the CH_4 specifics
      (offset size, -67.5° bin), instead of re-describing harmonisation and matching. EPICA paragraph trimmed likewise.
      The Law Dome smoothing detail is CH_4-specific, so left in
      (could move to the appendix with the optional smoothing figure, see Section 2).
    - CFC12 Step 2: now "the same as for N_2O, except that we do not require data in the most northern and southern boxes";
      CFC12/C_4F_10 Step 5: "the same as for N_2O, except that..." for the negative-value check.
    - C_4F_10 Step 4: dropped the "Like for other gases, the next step is to then extend..." filler.
    - Typos: "for for", "PCA analysis", "that to not", "EPIC" → "EPICA".
- [x] Skipped: we're not going to include a table (too terse to be useful to the reader, so we stick with the dot points).
      A redesigned version (keywords per step) replaces the old draft, still commented out, in case others want to reconsider.
      Create the table of differences from the base case (from `NOTES.md`).
      A draft is commented out at `methods.tex:274-295` (`table:methods-differences`),
      but its ZNTODO says that layout doesn't work (too much text per cell).
      Redesign it (e.g. short ticks/keywords per step, details left to the text)
      and decide whether it replaces or complements the per-gas summary list (`methods.tex:161-272`).
- [x] Unify language for "observational network" vs "observation network" and "ice core extension"
      (`methods.tex:347,1125`, both used throughout).
      Also hyphenation: "global-, annual-mean" vs "global- annual-mean". Done (methods, results; discussion already consistent):
    - "observation network" everywhere (was the majority, and what results/discussion use);
      "observational network(s)" and "observing networks" replaced.
      Compounds hyphenated: "observation network-based", "observation network-derived".
    - General Step 2 now defines the `interpolated observation network dataset' (the gap-free result of Step 2)
      and the `interpolated observation network period' (the time it covers,
      which is not the same as the period the observation network itself covers).
      "observational period" → "interpolated observation network period" (CO_2, CH_4 Step 4),
      "before the observation network (pre-1989)" → "before the interpolated observation network period",
      CFC12 table caption/footnote ("Obs. network years", footnote b) use these terms too.
      "observational record" kept where it means the observations in general (general Step 2, CO_2 Step 3),
      but "not already covered in the observational record" (CO_2 Step 4) → "by the interpolated observation network dataset".
    - "ice core extension" was only used once (CH_4 Step 4 PC list), now "the optimisation against ice cores".
    - "global-, annual-mean" everywhere ("global- annual-mean" replaced).
      Checked each "global-mean" on its own: kept where it really is a global-mean
      (global-mean surface air temperature, the derived global-mean monthly/yearly grids, NOAA products in results,
      the "Global-mean source" column header in the CFC12 table),
      changed to "global-, annual-mean" for the relative seasonality denominator (N_2O Step 3),
      the CH_4 ice-core optimisation free parameter (per year, `global_annual_mean_optimised` in `1103_ch4_extend-pcs.py`)
      and the Trudinger composite (CFC12 Step 4).
      UCI's timeseries in the results no longer called "global-mean" (see the UCI TODO in section 5).
    - Removed the two "unify language" TODOs from `methods.tex`.
    - Docstrings, comments and value-check descriptions in `src/local/` and `scripts/` updated to "observation network" too
      (identifiers and file names such as `observational_network_global_annual_mean_file` untouched).

## 2. Appendices

The build doesn't support appendices yet, but most of the figures already exist.

### Infrastructure

- [ ] Add an appendix option to `scripts/compile-gmd-template-based-latex.py`
      (the Copernicus template has a commented `\appendix ... \noappendix` block
      after `\codedataavailability`, `copernicus-latex-package/template_clean.tex:97-115`).
- [ ] Create `manuscripts/historical-ghg-forcing-for-cmip7/appendices.tex`
      and pass it from `scripts/create-cmip7-historical-ghgs-manuscript.sh`.
- [ ] Check appendix figures/tables get numbered A1, A2, ... (`\appendixfigures`/`\appendixtables`, `\noappendix`).
- [ ] Decide what goes in the appendix and what goes in a supplement (roughly 80 extra figures, see below).

### Methods appendix figures (already generated)

- [ ] Put a note as a comment in the latex that we use appendices based on copernicus's
      distinction between appendices ("all material required to understand the essential aspects of the paper")
      and supplementary (not required for understanding the paper i.e. surplus to requirements and "Supplementary material is reserved for items that cannot reasonably be included in the main text or as appendices")
      so that other authors know why we've done this.
      Link to https://publications.copernicus.org/for_authors/manuscript_preparation.html
- [ ] N_2O, CO_2, CH_4, CFC-12: move the draft captions from `methods.tex` into the appendix
      (`methods.tex:327-336,526-539,749-759,960-970`), use the `<gas-methods-appendix-figure>` tags,
      and remove the "update this when we get to appendix figures" ZNTODOs.
- [ ] The other 33 CFC-12-like gases: `<gas>-methods-figure` and `<gas>-methods-appendix-figure`
      are generated and in the manifest but not used anywhere in the text.
- [ ] The other C_4F_10-like gases (C_5F_12, C_6F_14, C_7F_16, cC_4F_8): `<gas>-methods-figure` is generated but not used.
- [ ] Check the CFC-12 appendix draft caption, "e) Regression between the first latitudinal gradient PC and CFC12 emissions",
      against the method (regression against total SSP2-4.5 emissions from RCMIP, `methods.tex:1035-1037`).
- [ ] Optional: CH_4 Law Dome smoothing figure (noise, windows, regression; `methods.tex:821-823`, "if requested by reviewers").

### Results appendix figures (already generated)

- [ ] Results figures exist for all 43 gases plus the equivalent species,
      but only CO_2, CH_4, N_2O, CFC-12, CFC-12-eq and HFC-134a-eq are in the text.
      Add the rest to an appendix (`results.tex:9`) and fill in "[TODO Appendix section and refs]" (`results.tex:286`).
- [ ] Possibly add CFC-11-eq (`scripts/create-cmip7-historical-ghgs-manuscript.sh` notes it isn't discussed yet).

### References from the main text into the appendix

- [ ] N_2O Step 3, "explain almost all of the variance (Figure [ZNTODO figure panel ref])" (`methods.tex:391`), pointing to the variance explained panel.
- [ ] CO_2 Step 3 variance claims (`methods.tex:574-576,597-600`), pointing to the variance explained panels (c, f).
- [ ] CO_2 Step 4 regressions (`methods.tex:607-612,667-688`), pointing to the CO_2 appendix panels e and h (ZNTODO at `methods.tex:688`).
- [ ] CH_4 and CFC-12 Step 4 regressions, pointing to panel e of their appendix figures.
- [ ] General Step 2 interpolation (`methods.tex:40-70`), pointing to panels a and b (best/worst-case interpolation).
- [ ] A pass over the whole methods section looking for other places a figure reference would help (ZNTODO at `methods.tex:688`).
- [ ] A pass over the whole methods section looking for places where we refer to a figure/value for justification for one gas but don't similarly for others (e.g. currently we don't explain why we keep two lat. grad. EOFs for CH4)

## 3. Broken references and build issues

- [x] `\ref{ssec:methods-c8f18-like}` is undefined (`methods.tex:1262`). It should be `ssec:methods-c8f18`.
- [x] `\ref{sssec:methods-co2-pca-derived-components}` is undefined (`methods.tex:572,582`).
      These sentences mean "same as N_2O", so they should point to `sssec:methods-n2o-pca`.
- [ ] CFC-12 results caption refers to `Figure \ref{fig:results-cfc12}a` (itself) for the obs. network (`results.tex:271`).
      It should be `fig:methods-cfc12`.
- [ ] "Float too large for page by 10pt": this is the CFC-12-like per-gas `table*` (`methods.tex:995-1020`).
- [ ] `caption` package "Unknown document class" warning (`main.log:2926`). Check that captions render as intended.
- [ ] `results.tex:109` reads `Figure\ref` (missing `~`).
- [ ] Results and methods sections have no `\label{sec:...}`. Add them if they'll be referenced.
- [x] Rebuild (done 2026-10-06 after the header fixes; no undefined references left).
- [x] "Label(s) may have changed. Rerun" persisted after repeated builds.
      Fixed by a third pdflatex pass after bibtex (`compile_latex` in `scripts/compile-gmd-template-based-latex.py`),
      which now also warns if latex still wants a rerun. LaTeX output goes to stderr.

## 4. Methods: content and correctness

- [ ] General Step 2 and CFC-12-like Step 2 say "at least five data points are required" for spatial interpolation
      (`methods.tex:49,980`), but the original run uses `MIN_POINTS_FOR_SPATIAL_INTERPOLATION = 4`
      with `if n < 4: skip`, i.e. at least four
      (`1001_n2o`, `1101_ch4`, `1201_co2`, `1301_sf6-like` `_interpolate-observational-network` notebooks).
      Fix the text and consider a value-check.
- [ ] N_2O Step 4: the harmonisation year is given as 1984, "the earliest year in our observational network derived global- annual-mean"
      (`methods.tex:437-438`), but the N_2O network runs from 1989 (`methods.tex:346,418`).
      This looks copied from CH_4. Fix it and add a value-check.
- [ ] CH_4 Step 4: EPICA is said to be in the -82.5° bin (`methods.tex:873`) and Law Dome in -67.5°.
      EPICA Dome C is about 75.1°S, on the bin edge. Check which bin the code used.
- [ ] CH_4 Step 4: "[TODO cite EPIC]" typo (`methods.tex:870`); check the first year of EPICA use is 154 (`methods.tex:864`).
- [ ] CO_2 Step 4: "we then extend emissions back to the 1750" (`methods.tex:610`) should be "extend the first PC back to 1750".
- [ ] CO_2 composite equation (`methods.tex:672-677`): the trailing `\\` before `\end{align}` adds an empty numbered line,
      and the `t'` and `c'` lines get their own numbers. Use `\nonumber` or `aligned`.
- [ ] Equation \ref{eq:native-resolution-sum} discussion: explain why each lat. gradient EOF (not just their sum) has zero area-weighted mean
      (`methods.tex:149,379`). The integral as written also doesn't show the area weighting.
- [ ] Relative seasonality: is average(seasonality)/average(global-mean) equal to average(seasonality/global-mean)?
      Convince ourselves, or note the bug (`methods.tex:409-410`).
- [ ] CFC-12 Step 2: the "lat grid bands are 15 degrees" zalue-check at `methods.tex:991` sits next to the 1979-2023 statement.
      It should be a `cfc12-obs-network-start/end` check.
- [ ] CFC-12 Step 4: figure out the Trudinger composite period and add a check (`methods.tex:1064`);
      check the harmonisation decline period (`methods.tex:1069`);
      check pre-industrial sources and pull references from M17 (`methods.tex:1076`);
      check the Velders citation and the 100-year convergence time in the table notes (`methods.tex:1011,1017`).
- [ ] CFC-12 Step 1: introduce the distinction between NOAA's combined HATS product and HATS flask data in the source lists (`methods.tex:938-941`).
- [ ] CO_2/CH_4 Step 1: check in-situ vs flask, and the implications for assuming the same scales (`methods.tex:505,730`).
- [ ] C_4F_10-like: add checks for the Droste first/last years and that Droste is zero in its first year
      (`methods.tex:1193,1205,1221-1233`); fill in "[TODO Droste first year]";
      calculate the ERF share of these gases (`methods.tex:1196-1197`); check the site coordinates (`methods.tex:1161`).
- [ ] CO_2 Step 4: check the 1850 value for the seasonality change PC extension (`methods.tex:183`).
- [ ] CFC-12-like/C_4F_10-like Step 5: consider a check that the negative-value reduction only removes rounding errors (`methods.tex:229,261`).
- [ ] Equivalent species (`methods.tex:1296-1331`): "it not all" should be "if not all";
      "CFC-11 captures" / "CFC12 captures" should be "CFC-11-eq" / "CFC-12-eq";
      double check the AR6 Table 7.SM.6 radiative efficiencies (`methods.tex:1325`; also see the "halving" worry at `results.tex:291-292`).
- [ ] ESGF section (`methods.tex:1333-1354`): explain the frequencies and grid labels; say where full DRS details are;
      cite input4mips-validation and the input4MIPs CVs page/user guide.
- [ ] Vertical dimension section (`methods.tex:1356-1361`) is a stub: summarise the suggested approach and give the M17 section number.
- [ ] Format `numpy.linalg.svd` (and other code names) as `\texttt{}` (`methods.tex:381`).
- [ ] Check M17's names for EOFs and PCs (`methods.tex:383,385`).
- [ ] Convert all remaining `ZNTODO zalue-check` / "add check" comments into real `value-check`s
      (`methods.tex:183,429,434,437,574,598,628,630,635,646,657,785,788,802,804,809,818,835,846,859,864,871,875,987,991,1064,1069`).

## 5. Results

- [ ] Uncommented placeholder text: "Differences over the period A-B are less than X~ppm." (`results.tex:58-59`).
- [ ] Go through all value-checks (`results.tex:41`) and fix the disabled ones:
      `co2-monthly-diff-from-maunoa-loa` (the difference is bigger than stated, `results.tex:49-50`),
      the 1940 dip / 1970 spike checks (`results.tex:80,82`), `cfc12-monthly-diff-from-noaa/agage` (`results.tex:237-238`),
      and the sign of the CFC-12 lat. gradient weakening check (`results.tex:246`).
- [ ] Add the missing checks: CO_2/CH_4/N_2O vs IGCC over the NOAA period (`results.tex:54,117,185`),
      the HFC-134a-eq lat. gradient (`results.tex:350`), the radiative efficiency source (`results.tex:25`).
- [ ] CH_4: update the IGCC comparison to treat the non-NOAA period properly (`results.tex:115`);
      UCI: check what the UCI timeseries is exactly (spatial region, averaging; it's probably not a global-mean)
      and make sure it's described correctly, its latitude bounds and where the data was retrieved from (`results.tex:122-126`);
      "Summit" ice core outlier: what is it, and does the paper explain it (`results.tex:139-141`, `comparison_data.py:486`).
- [ ] Check "IGCC is also based on NOAA data" (`results.tex:60,120,188`) and "the same source is used by IGCC" (`results.tex:244`).
- [ ] Check IGCC citations and processing (`results.tex:388`, `src/local/historical_ghg_forcing_for_cmip7/comparison_data.py:826`).
- [ ] Define AR6 once and cite it (`results.tex:301,357`); check "IGCC uses Hodnebrog" (`results.tex:302,358`).
- [ ] Plain-text gas names that the replacements won't catch: "CO2" (`results.tex:43`), "CH4" (`results.tex:106,132,153`).
- [ ] Heading vs text names: "CFC-12-eq"/"HFC-134a-eq" headings vs "CFC12-eq"/HFC134a in text. Check `replacements.yaml` covers all forms.
- [ ] Specify that piControl should repeat the 1850 values (`output-requirements.tex:117`; probably in results).
- [ ] Ask co-authors for other (ideally independent) comparison datasets (`results.tex:10-11`).

## 6. Abstract, introduction, output requirements, discussion, conclusion, statements

- [ ] Abstract: fill all [TODO X/Y/Z] numbers from `results.tex` (`abstract.tex:6-20`).
      Most already exist in results (e.g. CO_2 4.4 ppm around year 200, 2 ppm since 1850, <0.5 ppm since 1981).
      Ideally drive them from the same value-checks.
- [ ] Use either "radiative forcing" or "radiative effect" consistently (`abstract.tex:12`).
- [ ] Discussion claims all differences are "less than 0.03 W/m2" (`discussion.tex:7-8`),
      but results give about 0.05 W/m2 (CO_2, full period) and about 0.04 W/m2 (N_2O, full period). Reconcile.
- [ ] Check the end year: the abstract says datasets end in 2022, the obs networks run to 2023 (`methods.tex:346,549,768,837`),
      historical ends 2021-12 (`output-requirements.tex:31`).
- [ ] Introduction: cite Vaishali's paper (`introduction.tex:9`); find papers on how forcing generation evolved / the forcings task team
      (`introduction.tex:15-18`); cite other forcing papers and an overview of model inputs (`introduction.tex:22`).
- [ ] Output requirements: Funke et al. solar ref (`output-requirements.tex:16`); scenario paper cite (`:29`);
      check 0.01 W/m2 is about 1 ppm CO_2 (`:60`).
- [ ] Discussion: add the points flagged in methods/results/output requirements (`discussion.tex:4`):
    - more robust interpolation (`methods.tex:54-55`)
    - harmonisation is crude (`methods.tex:439-440,639-640`)
    - Menking CH_4 processing error and not using Menking pre-1850 (`methods.tex:792-797`, `output-requirements.tex:97-99`)
    - calibration scales not converted, and different scales of observations (`methods.tex:930`, `output-requirements.tex:102-106`)
    - constant poleward extrapolation as a weakness (`methods.tex:982`)
    - needing more rapidly updated emissions, or an approach without the regression; the need for a complete dataset back to year 1 (`methods.tex:1038-1047`)
    - fixing upstream processing to avoid small negative values (`methods.tex:1122`)
    - CH_4 seasonal cycle possibly too strong (`results.tex:132-134`)
    - C_8F_18 future iterations (`methods.tex:1257`)
    - annual updates (`abstract.tex:26`) and REF for forcings (`output-requirements.tex:111-115`);
      the "How to do extensions" notes in `NOTES.md` sketch the approach
      (optimise the existing lat. gradient/seasonality change against new network data rather than re-deriving or regressing PCs)
    - vertical resolution (`discussion.tex:60`)
- [ ] Conclusion: user guide link (`conclusion.tex:3-4`); fill the maximum-difference numbers and years (`conclusion.tex:15-18`);
      QUICCA full name / C3S paragraph if confirmed (`conclusion.tex:28-31`).
- [ ] Author contribution is a placeholder (`author-contribution.tex`); maybe auto-generate from `metadata.toml`.
- [ ] Acknowledgements is a placeholder: CMIP IPO, ESA funding (from the data files), plus the planning list in section 10.
- [ ] Competing interests: add a statement about the funding being sought (`competing-interests.tex:1`).
- [ ] Results intro: user guide citation (standalone Zenodo doc to start, `results.tex:16`).

## 7. Citations

Already in `references/references.bib`, just need wiring in:

- [ ] Droste et al. → `acp-20-4787-2020` (`methods.tex:240,250,255,1158,1169,1194,1206,1222,1230,1233`)
- [ ] NOAA CO_2 / CH_4 → `lan_atmospheric_co2_2025`, `lan_atmospheric_ch4_2025` (`methods.tex:729`, `results.tex:47,111,182`)
- [ ] AGAGE → `prinn_history_2000`, `prinn2018history` (`methods.tex:308,729,921`, `results.tex:128,190,239`)
- [ ] NOAA HATS CFC-12 → `noaa_hats_cfc12` (`methods.tex:919`); uncomment the existing `\citep`s
- [x] SSP2-4.5 / M2020 → `meinshausen_shared_2020` (C_8F_18 section)
- [ ] AR6 Ch7 SM → `IPCC_2021_WGI_Ch_7_SM` (`methods.tex:1326`, `results.tex:301,357`)
- [ ] Daniel et al. / WMO 2022 → `wmo_2022_ozone_ch7`? Check (`results.tex:242,260`)

Need adding to the bib:

- [ ] NOAA HATS N_2O, GAGE, ALE (`methods.tex:308-309`)
- [ ] scipy, numpy (`methods.tex:58,381`)
- [ ] Menking et al. 2025 (many places in `methods.tex` and `results.tex`)
- [ ] Law Dome - Mauna Loa merged record (`methods.tex:180,623-651`, `results.tex:79`)
- [ ] NOAA Mauna Loa (`results.tex:48`)
- [ ] PRIMAP-hist, with version (`methods.tex:609,611,783`)
- [ ] HadCRUT (`methods.tex:668,686`)
- [ ] Law Dome CH_4, NEEM, EPICA (`methods.tex:790-877`)
- [ ] UCI CH_4 (`results.tex:122`), CH_4 ice-core comparison paper (`results.tex:136`)
- [ ] IGCC citations (`results.tex:53,116,184`; `forster_indicators_2026` exists, check it's the right one)
- [ ] AR5 WG1 Ch. 8 Appendix 8.A (`results.tex:320`)
- [ ] input4mips-validation, input4MIPs CVs (`methods.tex:1350,1354`)
- [ ] Funke et al. solar, scenario in-prep paper, Vaishali's paper, intro forcing papers (section 6)

## 8. Nomenclature and typos

- [ ] Decide CFC12 vs CFC-12 in the source (abstract asks "[TODO check nomenclature]").
      `replacements.yaml` maps CFC12 → CFC-12 so the output is fine, but the source mixes both.
      Same for HFC134a vs HFC-134a.
- [ ] `replacements.yaml`: double check the mappings for minor species (`:1`);
      add the Halon (1202?) that is in the scenarios but not the historical data (`:28`).
- [ ] Typos: "accomodate" (`methods.tex:500,904`), "that that" (`:460`), "for for" (`:773,1025`),
      "PCA analysis" (`:773,1025,1154`), "to not source" (`:1072`), "it not all" (`:1299`), "Next we consider CH_4" missing full stop (`:722`),
      "differencs" (`results.tex:29`), "minorly" (`results.tex:371`), "radiative focing" (`conclusion.tex:16`), "foward" (`conclusion.tex:27`),
      "a key criteria" (`conclusion.tex:24`), "betwen" (`output-requirements.tex:53`), "than an another" (`output-requirements.tex:86`).

## 9. `src/local`

Manuscript-relevant:

- [ ] `historical_ghg_forcing_for_cmip7/comparison_data.py:486`: check the CH_4 ice-core supplement site metadata (Summit etc.).
- [ ] `historical_ghg_forcing_for_cmip7/comparison_data.py:826`: check the IGCC processing (release v6.4.0).
- [ ] `historical_ghg_forcing_for_cmip7/equivalent_species.py:77`: use pint for unit handling
      (relevant to the radiative efficiency "halving" worry).
- [ ] `historical_ghg_forcing_for_cmip7/results_figure.py:168-175`: `SHOW_OUTPUT_AT_COMPARISON_LATITUDES` is off.
      Decide whether the fairer site comparison is needed (e.g. for the Mauna Loa claim at `results.tex:48-51`).
- [ ] `historical_ghg_forcing_for_cmip7/co2_methods_figure.py:944,982`, `ch4_methods_figure.py:499`: remove the hard-coding.
- [ ] `src/local/historical_ghg_forcing_for_cmip7/__pycache__` has a stale `sf6_like_methods_figure` .pyc (module now `cfc12_like_methods_figure`). Clean up.

General tooling (low priority, not manuscript-blocking, skip all of these for now):

- [ ] `data_loading.py:140,298,333,337,353,461,494,498,520`: ESGF result checking, loading logic, hard-coded paths.
- [ ] `esgf/` TODOs: uniqueness constraints, process count, mapping organisation, node handling, URL sorting, chunk sizes, retry logic.

## 10. From `manuscript-planning.md` (historical manuscript)

Moved from `manuscript-planning.md`. Notes are from reading Meinshausen et al. 2017 (M17):
https://gmd.copernicus.org/articles/10/2057/2017/gmd-10-2057-2017.pdf

- [ ] Invite all the data producers to be co-authors
- Abstract
    - [ ] GHGs are a major driver of past climate change, therefore key for model simulations
    - [ ] Updated GHG concentrations for CMIP7
    - [ ] Once again a composite of multiple input sources, now covering [time range]
    - [ ] Put ERF due to GHGs in; note 'highest ever' or whatever, and the largest increases in recent times
    - [ ] Available for modelling teams to use: "Finally, we describe where to access the data and provide a summary of key data changes compared to CMIP6."
- Content
    - [ ] Mole fraction in dry air vs. mole fraction in the real atmosphere
    - [ ] Different scales (ignored this time, noise in the scheme of things)
    - [ ] Historical experiment ends in 2021 even though this dataset goes to 2022 (and we hope to extend further in future)
    - [ ] Changes compared to CMIP6 (and the reasons)
    - [ ] Comparisons with other datasets (lots of 'time pressure tied our hands' in here)
    - [ ] If you need explanation/justification for methods, refer back to M17
    - [ ] Missing halon in the historical dataset: note the ERF difference is tiny, so not ideal, but not a reason to re-write/re-run
    - [ ] Don't compare the seasonal cycle and latitudinal gradient from CMIP6 ESMs (unless there's way more time than expected)
    - [ ] "Given the negligible radiative forcing from ..., this uncertainty does not affect the overall results."
    - [ ] If time: build a portal for visualising/exploring/accessing results
- Introduction
    - [ ] CMIP context
    - [ ] Unique requirements of CMIP, therefore the goal of this study
- Methods
    - [ ] Flowchart figure for the overall idea, with sub-panels for key variants/details (probably only 5).
          Maybe a table too, or just text given how big the tables are in M17
          (see `tmp-figure-scribbles/hist-methodology-figure-*.jpeg` and the table decision in section 1)
    - [ ] Build out based on notebooks
    - Compared to the CMIP6 workflow
        - [ ] Removed optimisation step
        - [ ] Updated/captured polynomial smoothing (break out a package?)
        - [ ] Removed N_2O interpolation from 1966-1987?
- Results
    - [ ] We don't do uncertainties (future work)
    - [ ] Compare to CMIP6 throughout
    - [ ] Compare to other studies here?
    - CO_2: global-mean, lat. grad. (compare with M17), seasonality
        - [ ] Make plots that show all components of M17 Fig. 9, but don't put them in the paper (outreach product)
        - [ ] Main paper needs a plot with the relevant components: global-mean, global-mean monthly steps, lat. gradient, seasonality
    - CH_4: global-mean, lat. grad., seasonality
        - [ ] Compare with https://www.nature.com/articles/s41586-026-10938-1
    - N_2O: global-mean, lat. grad., seasonality
    - ODSs
        - [ ] Reproduce/check all the pre-industrial choices, extrapolations, hard-coded zero seasonality and lat. gradient etc.
- Data format and recommendations
    - [ ] Use input4MIPs CVs text
    - [ ] Point to the forcings implementation and forcing usage recording docs (i.e. how to record what f1 means), and papers as secondary
- Discussion
    - [ ] Explain the differences from CMIP6 here?
    - Limitations (all relatively small given the use case)
        - [ ] Focussed on CMIP; don't use for inversion studies etc., you need a different product (dup from M17)
        - [ ] No vertical or longitudinal resolution (dup from M17)
        - [ ] Do we want observed concs or concs that reflect background concs (dup from M17)
        - [ ] Hybrid nature of calibration scales (dup from M17)
        - [ ] Uncertainty (dup from M17)
- Acknowledgements
    - [ ] Everyone collecting raw data, particularly for their openness
    - [ ] CMIP organisers and all involved
    - [ ] Forcings TT panel
    - [ ] Direct funding acknowledgements (ESA)
- Data access
    - [ ] ESGF
    - [ ] Zenodo
