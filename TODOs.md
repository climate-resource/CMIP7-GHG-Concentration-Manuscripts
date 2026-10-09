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

The build supports appendices (see Infrastructure), and most of the figures already exist.

### Infrastructure

- [x] Add an appendix option to `scripts/compile-gmd-template-based-latex.py`
      (the Copernicus template has a commented `\appendix ... \noappendix` block
      after `\codedataavailability`, `copernicus-latex-package/template_clean.tex:97-115`).
      Done: `--appendix` (can be given more than once, like `--section`).
      The files are inserted at a new `% @appendix-start@` tag in `template_clean.tex`
      (after `\codedataavailability`, before `\authorcontribution`, where the template's commented block is),
      wrapped in `\appendix ... \noappendix` by the script, so the appendix files only hold `\section`s and content.
- [x] Create `manuscripts/historical-ghg-forcing-for-cmip7/appendices.tex`
      and pass it from `scripts/create-cmip7-historical-ghgs-manuscript.sh`.
      Done: "Appendix A: Detailed methods figures" (`app:methods-figures`)
      and "Appendix B: Detailed results figures" (`app:results-figures`).
      The first figures are in (see below): Appendix A has the N_2O, CO_2, CH_4, CFC-12 methods appendix figures
      and the SF_6 methods and methods appendix figures (A1-A6),
      Appendix B the SF_6 and CFC-11-eq results figures (B1, B2).
- [x] Check appendix figures/tables get numbered A1, A2, ... (`\appendixfigures`/`\appendixtables`, `\noappendix`).
      Done: with the `manuscript` class option `\appendix` resets the figure and table counters per appendix section,
      so figures/tables put inside an appendix `\section` are numbered automatically
      (the template's "Option 1").
      We don't use `\appendixfigures`/`\appendixtables`, despite the template mentioning them:
      the template only asks for them if all floats are put after the reference list ("Option 2"),
      and they break the numbering with Option 1.
      Tested with `\appendixfigures` after each appendix `\section`:
      the Appendix B figures came out as C1, C2 (it steps the section counter on top of `\section`'s step).
      This is explained in a comment at the top of `appendices.tex`.
      Checked with a throwaway build with two extra appendices holding figures and tables:
      captions and `\ref`s came out as Figure C1, C2, D1 and Table C1, D1,
      the headings as "Appendix C: ...", and the main-text figures stayed 1-12.
      Note: the manuscript class option floats tables to the end of the document, appendix tables included.
- [x] Decide what goes in the appendix and what goes in a supplement (roughly 80 extra figures, see below).
      Decision: everything goes in appendices, no supplement.
      Copernicus reserves the supplement for "items that cannot reasonably be included in the main text or as appendices"
      (https://publications.copernicus.org/for_authors/manuscript_preparation.html),
      and figures can reasonably be included as appendices
      (this is the reasoning to put in the latex comment, see the first item under "Methods appendix figures").
      Revisit if the per-gas figures make the PDF unwieldy, e.g. if the editor asks us to move them.

### Methods appendix figures (already generated)

- [x] Put a note as a comment in the latex that we use appendices based on copernicus's
      distinction between appendices ("all material required to understand the essential aspects of the paper")
      and supplementary (not required for understanding the paper i.e. surplus to requirements and "Supplementary material is reserved for items that cannot reasonably be included in the main text or as appendices")
      so that other authors know why we've done this.
      Link to https://publications.copernicus.org/for_authors/manuscript_preparation.html
      Done: at the top of `appendices.tex`, quoting the page (checked 2026-10-07).
      Note: on that page, "All material required to understand the essential aspects of the paper"
      is what should go in the main text, not the appendices.
      Appendices are for "Additional figures, tables, ... which are not critical to support the conclusion of the paper,
      but which provide extra detail and/or support useful for experts in the field
      and whose inclusion in the main text would disrupt the flow", so the comment quotes that instead.
- [x] N_2O, CO_2, CH_4, CFC-12: move the draft captions from `methods.tex` into the appendix
      (`methods.tex:327-336,526-539,749-759,960-970`), use the `@gas-methods-appendix-figure@` tags,
      and remove the "update this when we get to appendix figures" ZNTODOs.
      Done: Figures A1-A4 (`fig:methods-appendix-<gas>`), captions start "Working behind Figure \ref{fig:methods-<gas>}".
      Where the draft captions were, `methods.tex` now has a test cross-ref to the appendix figure
      (see "References from the main text into the appendix").
- [x] The other 33 CFC-12-like gases: `<gas>-methods-figure` and `<gas>-methods-appendix-figure`
      are generated and in the manifest but not used anywhere in the text.
      Done: Appendix A2 "Other gases processed like CFC-12" (Figures A5-A70, `fig:methods-<gas>`, `fig:methods-appendix-<gas>`).
      These (and the other groups of per-gas appendix figures, see below) are generated,
      rather than written out by hand for each gas:
      `src/local/historical_ghg_forcing_for_cmip7/appendix_figures.py` makes the figures,
      whose captions are just "Like Figure \ref{<the CFC-12/C_4F_10 figure of the same kind>}, except for <gas>.".
      `generate-tex-inputs.py` writes one file per group to `build/historical-ghg-forcing-for-cmip7/appendix-figures/`
      and lists them under the manifest's `inline` entry.
      The compilation puts each file's text in place of its `@...-appendix-figures@` tag in `appendices.tex`
      (`--inline-file tag=path`, or the manifest's `inline` entry: "replace this tag with the text of this file").
      The generated tables use the same mechanism (`--table-file` and the manifest's `tables` entry are gone),
      and tags in comments are left alone, both when inlining files and when replacing figure tags,
      so comments can mention tags.
      Tags are written `@tag@` (e.g. `@cfc12-like-per-gas-table@`, `% @body-start@` in the template), not `<tag>`:
      `<`/`>` are redirections in the shell, so `<tag>` had to be quoted on the command line,
      whereas `@` means nothing to bash, zsh or latex (and is the autoconf/CMake substitution convention).
      The manifest's `{tag: path}` entries are passed through as pairs; only CLI values are split on `=`.
      The compilation fails if any `@tag@` is left (outside comments) after the figures and files are put in,
      rather than the tag ending up in the PDF as text.
      Which gases go in each group follows from the gases the build script draws figures for,
      minus any figure whose tag the manuscript files (`--manuscript-file`, the body sections and appendices)
      already include outside a comment, e.g. CFC-12's and C_4F_10's.
      The build script now keeps the body sections and appendices in arrays,
      so the input generation and the compilation are given the same files.
      Each gas starts on a new page: without the `\clearpage`s latex fails with "Too many unprocessed floats"
      (and we can't add packages like `morefloats` for Copernicus).
      Appendix A is now split into subsections: A1 "N_2O, CO_2, CH_4 and CFC-12" (hand-written, Figures A1-A4),
      A2 (CFC-12-like) and A3 (C_4F_10-like).
      Appendix headings write "CFC-12", not "CFC12":
      the replacement's `\nobreakdash` breaks the build in an appendix heading.
- [x] The other C_4F_10-like gases (C_5F_12, C_6F_14, C_7F_16, cC_4F_8): `<gas>-methods-figure` is generated but not used.
      Done: Appendix A3 "Other gases processed like C_4F_10" (Figures A71-A74), generated as above.
- [x] Check the CFC-12 appendix draft caption, "e) Regression between the first latitudinal gradient PC and CFC12 emissions",
      against the method (regression against total SSP2-4.5 emissions from RCMIP, `methods.tex:1035-1037`).
      The caption is now in `appendices.tex` (with a ZNTODO), and the SF_6 appendix caption copies it, so check that too.
      Done: the figure's emissions are the original run's `historical_emissions.csv`,
      which `notebooks/010y_compile-historical-emissions/0109_compile-complete-dataset.py` in the original run
      builds from the RCMIP v5.1.0 SSP2-4.5 World emissions, so it matches the method.
      Panel e) now says "Regression between the first latitudinal gradient PC and total CFC12 emissions."
      (the text has the detail). The generated captions point to the CFC-12 one, so don't repeat it.
- [ ] Optional: CH_4 Law Dome smoothing figure (noise, windows, regression; `methods.tex:821-823`, "if requested by reviewers").

### Results appendix figures (already generated)

- [x] Results figures exist for all 43 gases plus the equivalent species,
      but only CO_2, CH_4, N_2O, CFC-12, CFC-12-eq and HFC-134a-eq are in the text.
      Add the rest to an appendix (`results.tex:9`) and fill in "[TODO Appendix section and refs]" (`results.tex:286`).
      Done: Appendix B is split into B1 "Gases processed like CFC-12" (Figures B1-B33, generated),
      B2 "Gases processed like C_4F_10" (B34-B38: C_4F_10 hand-written, as its results aren't in the main text
      and its caption has to say the input data are Droste's, the others generated and pointing to it)
      and B3 "C_8F_18 and CFC-11-eq" (B39, B40, hand-written: their captions are one-offs,
      neither has input data to show and IGCC has no CFC-11-eq).
      The ZNTODO at `results.tex:9` is removed and the placeholder now reads "please see Appendix \ref{app:results-figures}".
      The C_4F_10 results caption says "[TODO cite Droste et al] values" like the methods (see the citations section).
      The results text doesn't discuss any of these gases (yet).
- [x] Possibly add CFC-11-eq (`scripts/create-cmip7-historical-ghgs-manuscript.sh` notes it isn't discussed yet).
      Done: added to the results gases in the build script, figure generated, in Appendix B (Figure B40, `fig:results-cfc11eq`).
      Its caption notes that IGCC has no CFC-11-eq.
      The results text only points to it (no discussion of it yet).

### References from the main text into the appendix

- [x] Review the test cross-refs, each marked with
      "% TODO: review this cross-ref, this is just here to test the linking between main and appendices e.g. numbering."
      They are: a sentence after each of the N_2O, CO_2, CH_4 and CFC-12 methods figures pointing to its appendix figure
      (the N_2O one also points to Appendix A), one after the CFC-12 one pointing to the SF_6 figures (now A67, A68),
      one at the end of the CFC-12 results section pointing to the SF_6 results (now B32)
      and one in the equivalence datasets intro pointing to CFC-11-eq (now B40, Appendix B).
      Replace them with the real references below (or keep them, reworded).
      Done: all removed. The CFC-12 methods intro points to Appendix A2 (other CFC-12-like gases)
      and C_4F_10 Step 1 to Appendix A3 (other C_4F_10-like gases), so each group of appendix figures is pointed to.
      The SF_6 sentences are gone (covered by these pointers and the "species by species" pointer to Appendix B).
      The CFC-11-eq sentence says it covers the gases in CFC12-eq and HFC134a-eq except CFC12 (not their sum).
- [x] N_2O Step 3, "explain almost all of the variance (Figure [ZNTODO figure panel ref])" (`methods.tex:391`), pointing to the variance explained panel.
- [x] CO_2 Step 3 variance claims (`methods.tex:574-576,597-600`), pointing to the variance explained panels (c, f).
      Both the latitudinal gradient (f) and seasonality change (c) claims point to their panel.
- [x] CO_2 Step 4 regressions (`methods.tex:607-612,667-688`), pointing to the CO_2 appendix panels e and h (ZNTODO at `methods.tex:688`).
- [x] CH_4 and CFC-12 Step 4 regressions, pointing to panel e of their appendix figures.
- [x] General Step 2 interpolation (`methods.tex:40-70`), pointing to panels a and b (best/worst-case interpolation).
      Done in each per-gas Step 2 (N_2O, CO_2, CH_4, CFC-12) rather than in the general approach,
      which comes before any figure.
- [x] A pass over the whole methods section looking for other places a figure reference would help (ZNTODO at `methods.tex:688`).
      Done: every panel of every main-text methods figure (and the PC panels of the appendix figures)
      is now pointed to from the text that describes it:
      Step 2 binning → obs. counts (panel b), Step 3 → global-mean, seasonality, lat. gradient EOF and PC panels,
      Step 4 → extended global-mean and extended PC panels (N_2O, CO_2, CH_4, CFC-12, C_4F_10, C_8F_18).
      The CFC-12 appendix caption now says "PC", not "PCs" (CFC-12-like gases keep only one).
- [x] A pass over the whole methods section looking for places where we refer to a figure/value for justification for one gas but don't similarly for others (e.g. currently we don't explain why we keep two lat. grad. EOFs for CH4)
      Done:
    - CH_4 Step 3: now says it keeps two lat. gradient EOFs because they explain almost all of the variance (Figure A3c).
    - CFC-12-like Step 3 said "the same as for N_2O", but these gases keep one lat. gradient EOF, not two
      (`lat_gradient_n_eofs_to_use: 1` in the original run's config for all 34).
      Step 3 now says so and justifies it (Figure A4c, and the other gases' panels in Appendix A2):
      the first EOF explains more than 90% of the variance for every gas except HFC-236fa (82%),
      from `data/raw/cmip-ghg-concentration-generation/v1.0.0/manuscript-outputs/*_lat-gradient-variance-explained.csv`.
      The CFC-12 intro and the per-gas summary list now include Step 3 as a difference from the base case.
    - CH_4 Step 4: the second PC being kept constant now points to N_2O (where M17 is the justification).
    - Each new variance claim has a `% ZNTODO: add check of this claim` (like N_2O's and CO_2's, see section 4).
    - Checked and fine: interpolation panels and regression panels are pointed to for every gas that has them;
      CO_2 doesn't need N_2O's year 5 → year 1 extrapolation note (Menking et al. CO_2 starts in year 1);
      harmonisation offsets are given for N_2O, CO_2 and CH_4 (the CFC-12-like Trudinger ones are per gas, so not given).
    - Not resolved: CO_2 Step 4 gives no reason for taking the second lat. gradient PC linearly to zero in 1850
      (N_2O and CH_4 keep it constant). Find the reason (M17?) and add it.

## 3. Broken references and build issues

- [x] `\ref{ssec:methods-c8f18-like}` is undefined (`methods.tex:1262`). It should be `ssec:methods-c8f18`.
- [x] `\ref{sssec:methods-co2-pca-derived-components}` is undefined (`methods.tex:572,582`).
      These sentences mean "same as N_2O", so they should point to `sssec:methods-n2o-pca`.
- [x] CFC-12 results caption refers to `Figure \ref{fig:results-cfc12}a` (itself) for the obs. network (`results.tex:271`).
      It should be `fig:methods-cfc12`.
- [ ] "Float too large for page by 10pt": this is the CFC-12-like per-gas `table*` (`methods.tex:995-1020`).
      Being too tall is why it was pushed to the end of the document.
      A `\clearpage` at the end of `methods.tex` now flushes it there instead (and one around each appendix section
      keeps the appendix figures in their own appendix), but the warning remains.
      Leaving the warning for now: check with the editor whether they want this long table
      or for us to split it over multiple pages (ZNTODO next to the table in `methods.tex`).
- [x] `caption` package "Unknown document class" warning (`main.log:2926`). Check that captions render as intended.
      Harmless: `copernicus.cls` loads `caption` itself, then immediately restores its own `\@makecaption`.
      Checked in the PDF (2026-10-07): captions render in the Copernicus style ("**Figure 7.** ...").
- [x] `results.tex:109` reads `Figure\ref` (missing `~`).
      Now `Figure~\ref`. All `Figure(s)`/`Table`/`Section(s)`/`Equation`/`Appendix`/`and` `\ref`s
      in the manuscript (and the generated appendix captions) now use `~` too.
- [x] Results and methods sections have no `\label{sec:...}`. Add them if they'll be referenced.
      Added `sec:methods` and `sec:results` (not referenced yet).
- [x] Rebuild (done 2026-10-06 after the header fixes; no undefined references left).
- [x] "Label(s) may have changed. Rerun" persisted after repeated builds.
      Fixed by a third pdflatex pass after bibtex (`compile_latex` in `scripts/compile-gmd-template-based-latex.py`),
      which now also warns if latex still wants a rerun. LaTeX output goes to stderr.

## 4. Methods: content and correctness

- [x] General Step 2 and CFC-12-like Step 2 say "at least five data points are required" for spatial interpolation
      (`methods.tex:49,980`), but the original run uses `MIN_POINTS_FOR_SPATIAL_INTERPOLATION = 4`
      with `if n < 4: skip`, i.e. at least four
      (`1001_n2o`, `1101_ch4`, `1201_co2`, `1301_sf6-like` `_interpolate-observational-network` notebooks).
      Fix the text and consider a value-check.
- [x] N_2O Step 4: the harmonisation year is given as 1984, "the earliest year in our observational network derived global- annual-mean"
      (`methods.tex:437-438`), but the N_2O network runs from 1989 (`methods.tex:346,418`).
      This looks copied from CH_4. Fix it and add a value-check.
- [x] CH_4 Step 4: EPICA is said to be in the -82.5° bin (`methods.tex:873`) and Law Dome in -67.5°.
      EPICA Dome C is about 75.1°S, on the bin edge. Check which bin the code used.
    - [x] the code (`1104_ch4_extend-global-annual-mean`) uses the nearest bin centre.
          EPICA is at -75.0025°, 0.0025° south of the bin edge, so -82.5° is right.
          The text now says so
    - [x] value-checks `ch4-epica-lat`, `ch4-epica-lat-bin` and `ch4-law-dome-lat-bin`
          (data loaded via `zenodo_missing.py`)
- [x] CH_4 Step 4: "[TODO cite EPIC]" typo (`methods.tex:870`); check the first year of EPICA use is 154 (`methods.tex:864`).
    - [x] typo fixed
    - [x] check the first year of EPICA use is 154.
          The smoothed Law Dome data starts in 154 and EPICA is used for years 1 to 153
          (`np.arange(1, law_dome_start_year)` in `1104_ch4_extend-global-annual-mean`).
          The text said "year 1 to year 154", now "year 1 to year 153".
          Value-checks `ch4-law-dome-start-year` and `ch4-first-year`
- [x] CO_2 Step 4: "we then extend emissions back to the 1750" (`methods.tex:610`) should be "extend the first PC back to 1750".
- [x] CO_2 composite equation (`methods.tex:672-677`): the trailing `\\` before `\end{align}` adds an empty numbered line,
      and the `t'` and `c'` lines get their own numbers. Use `\nonumber` or `aligned`.
      The template says to use align, so we assume this is journal standard.
- [x] Equation \ref{eq:native-resolution-sum} discussion: explain why each lat. gradient EOF (not just their sum) has zero area-weighted mean
      (`methods.tex:149,379`). The integral as written also doesn't show the area weighting.
    - [x] explained in N_2O Step 3: each year's residuals have zero area-weighted mean
          and each EOF with a non-zero singular value is a linear combination of them;
          Step 6 points there
    - [x] integral shows the area weighting, now normalised ($\int a(l) dl = 1$)
    - [x] value-check `lat-gradient-eofs-area-weighted-mean` over all 42 gases' `allyears` EOFs
          (actual value about 6e-14, relative to the EOF's magnitude)
- [x] Relative seasonality: is average(seasonality)/average(global-mean) equal to average(seasonality/global-mean)?
      Convince ourselves, or note the bug (`methods.tex:409-410`).
    - [x] they aren't equal, but it isn't a bug.
          The code's ratio of averages is the average of each year's ratio, weighted by each year's global-, annual-mean.
          Multiplied back out (Step 5), it reproduces the observed average seasonality over the observation network period exactly
          (when the global-, annual-mean over that period is the observation network's, as for N_2O and CH_4).
          The text now says this; value-check `n2o-seasonality-reproduces-observed-average`
    - [x] recomputed both for the 36 gases with relative seasonality (in the bundle's pixi env; reproduces the saved files exactly).
          Max difference relative to the seasonality's max magnitude:
          N_2O 0.5%, CH_4 0.6%, CFC-12 19%, SF_6 37%, HFC-134a 74%, HCFC-141b 180%
    - [x] decide whether the sensitivity of fast-growing gases' seasonality to this choice
          belongs in the discussion of limitations
- [x] CFC-12 Step 2: the "lat grid bands are 15 degrees" zalue-check at `methods.tex:991` sits next to the 1979-2023 statement.
      It should be a `cfc12-obs-network-start/end` check.
    - [x] 1979-2023 statement now has real `cfc12-obs-network-start/end` value-checks
    - [x] "15 degrees" now has value-check `native-grid-lat-band-width`
          (also checks the bands are uniform and cover the globe).
          Every other "15\textdegree" mention in the manuscript (all in `methods.tex`, incl. the four figure captions
          and the Step 2 binning grid, which has the same latitude bands) now has the same tag
- [x] CFC-12 Step 4: figure out the Trudinger composite period and add a check (`methods.tex:1064`);
      check the harmonisation decline period (`methods.tex:1069`);
      check pre-industrial sources and pull references from M17 (`methods.tex:1076`);
      check the Velders citation and the 100-year convergence time in the table notes (`methods.tex:1011,1017`).
    - [x] Trudinger composite period: harmonised Trudinger from 1901 to the year before the network starts
          (1901-2007 for CF_4 and C_2F_6, 1901-2005 for C_3F_8), network from then on
          (also in 2008-2014, which Trudinger covers too; the old text implied otherwise).
          Value-checks `trudinger-start-year` and `trudinger-end-year`
    - [x] decline period: 100 years, from `n_transition_years=100` in `1304_sf6-like_create-global-annual-mean` (`trudinger-harmonisation-transition-years`),
          also tagged in table note b. Text now says the offset declines over the 100 years *before* the harmonisation year
    - [x] pre-industrial values: all match M17 (Sects. 3.4 and 3.5) except CCl_4 (0 vs M17's 0.025 ppt,
          in line with Walker et al. 2000 per M17). Text now says which gases are non-zero and why
          (value-check `cfc12-like-non-zero-pre-industrial-gases`), with M17's underlying references as `[TODO cite ...]`:
          Velders and Daniel 2014, Worton et al. 2006, Aucott et al. 1999, Trudinger et al. 2004, Muhle et al. 2010, Walker et al. 2000
    - [x] pre-industrial *years* aren't given in M17.
          They are close to, but not the same as, the year M17's (CMIP6) concentrations start to rise
          (e.g. CHCl_3 1940 vs ~1925, HFC-23 1950 vs ~1929, HCFC-141b 1950 vs ~1989).
          Decide how to describe where they come from (the table's "Pre-industrial: source" reads as if M17 gives them).
          Done: text says they don't come from M17, are close to the original sources but some don't match exactly
          (time pressure), with minimal impact; the caption says the source is for the value and the years are ours
    - [x] Velders citation: `velders_hfcs_2022` is the right paper (it goes with the Zenodo dataset the run used).
          Table note a's 1988 and 1980 now have value-checks
    - [x] Velders' data start in 1990 (raw and processed; several HFCs are zero in 1990),
          so the 1980 pre-industrial year is our choice, consistent with Velders, rather than a value Velders gives.
          Decide whether the table / note should say so.
          Done in the Step 4 text: 1980 is a pragmatic choice to make the fit work, which does no harm
          because Velders is zero in 1990 for all of these HFCs except HFC-143a (0.5 ppt).
          Value-checks `velders-first-year`, `velders-first-year-zero-except-hfc143a`, `velders-first-year-hfc143a`
- [x] CFC-12 Step 1: introduce the distinction between NOAA's combined HATS product and HATS flask data in the source lists (`methods.tex:938-941`).
    - [x] the generated list already labels them `NOAA HATS combined` (CCl_4, CFC-11, CFC-113, CFC12, SF_6) and `NOAA HATS flask`.
          Added prose before the list (from the draft comment) explaining the two and that these are the labels;
          it doesn't repeat the gas names, so no value-check needed
- [x] CO_2/CH_4 Step 1: check in-situ vs flask, and the implications for assuming the same scales (`methods.tex:505,730`).
    - [x] the original run uses NOAA surface flask and in-situ data for CO_2;
          NOAA surface flask and in-situ, AGAGE GC-MD, GAGE and ALE for CH_4. It applies no scale conversions
          (unlike M17's x1.0003 for AGAGE CH_4). The visible `[ZNTODO ...]`/`[TODO ...]` notes are replaced by what's used
    - [x] M17-based notes on scales added as comments in each gas' Step 1
          (N_2O compatible; CO_2 all NOAA; CH_4 x1.0003; halocarbons quote from M17, other factors 0.99-1.2)
    - [x] general approach, Step 1: not converting / not exploring scales is a limitation, but fine for our purposes (cites M17 Sect. 6.4)
    - [x] discussion: new "fifth" limitation paragraph on calibration scales
- [x] CO_2 Step 1: check more carefully that NOAA's CO_2 surface flask and in-situ data are on the same calibration scale
      (we think so, see the comment in CO_2 Step 1).
    - [x] both are on WMO X2019: the in-situ files' headers say `scale : CO2_X2019`,
          the flask files don't say but NOAA's flask README (section 6) does. Comment in CO_2 Step 1 updated
- [x] C_4F_10-like: add checks for the Droste first/last years and that Droste is zero in its first year
      (`methods.tex:1193,1205,1221-1233`); fill in "[TODO Droste first year]";
      calculate the ERF share of these gases (`methods.tex:1196-1197`); check the site coordinates (`methods.tex:1161`).
    - [x] Droste covers 1934-2018 for all five gases at both sites; output runs to 2022.
          "[TODO Droste first year]" filled in (1934). Value-checks `droste-first-year`, `droste-last-year`, `c4f10-like-last-year`
          (also on the general approach's "(2018) to 2022")
    - [x] Droste isn't exactly zero in 1934, but effectively zero (at most 4.8e-6 ppt; the original run requires < 1e-5).
          Text now says "effectively zero (less than 1e-5 ppt)". Value-check `droste-max-first-year-value`
    - [x] ERF: these gases together are 8.8e-4 W/m^2 in 2022 (from our concentrations: change since 1750 x AR6 radiative efficiency).
          Text says "less than 0.001 W/m^2". Value-check `c4f10-like-erf-2022`.
          No share of total ERF: we'd need an external total, and we decided to stick to estimates from our own values
    - [x] site latitudes from the data: Cape Grim -40.6833, Tacolneston 52.5127.
          Text had 41S, 52N (52.5 rounds to 53), now 40.7S, 52.5N with value-checks.
          Longitudes (145E, 1E) aren't in the processed data (from the Droste paper, comment in the tex), so unchecked
- [x] CO_2 Step 4: check the 1850 value for the seasonality change PC extension (`methods.tex:183`).
    - [x] the PC is constant up to and including 1850 and varies from 1851
          (`1205_co2_extend-seasonality-change-pcs` holds the composite constant before HadCRUT starts).
          Value-check `co2-seasonality-change-regression-start-year`, in the general approach and CO_2 Step 4
- [x] CFC-12-like/C_4F_10-like Step 5: consider a check that the negative-value reduction only removes rounding errors (`methods.tex:229,261`).
    - [x] it doesn't only remove rounding errors (`1305_sf6-like_...` / `1405_c4f10-like_create-pieces-for-gridding`).
          The latitudinal gradient is scaled so its most negative value is half the global-mean (zero where that is zero)
          in months where the two would give negative values, and the seasonality is capped at 35% of the global-mean.
          Affects nearly every CFC12-like gas in the years after its pre-industrial year (at most 44 years after, except HFC-152a),
          HFC-152a's latitudinal gradient in some months 1991-2023, HFC-236fa's seasonality every year 1996-2023,
          and the C_4F_10-like gases' latitudinal gradient 1934-2001.
          Global-, annual-means are unchanged (one factor per month / per year).
    - [x] text in the general approach, the section intros and both Step 5s rewritten to say this.
          Value-checks detect the scaled-down months/years from the output files (match the executed notebooks exactly)
          and copy the 0.5 and 35% constants from the notebooks
- [x] Equivalent species (`methods.tex:1296-1331`): "it not all" should be "if not all";
      "CFC-11 captures" / "CFC12 captures" should be "CFC-11-eq" / "CFC-12-eq";
      double check the AR6 Table 7.SM.6 radiative efficiencies (`methods.tex:1325`; also see the "halving" worry at `results.tex:291-292`).
    - [x] "it not all" → "if not all"
    - [x] "CFC-11 captures" → "CFC-11-eq captures"
    - [x] "CFC12-eq captures": `replacements.yaml` renders it as CFC-12-eq; the source spelling is the Section 8 decision
    - [x] "HFC134-eq captures" → "HFC134a-eq"
    - [x] AR6 radiative efficiencies: all 44 we use match AR6 WG1 Ch. 7 SM Table **7.SM.7** (checked against the PDF).
          Table 7.SM.6 is the carbon cycle response function, so the text (and the original notebook's comment) cited the wrong table;
          the text now says 7.SM.7
    - [x] the "halving" in `results.tex` is real: CH_3Cl is 0.01 W/m^2/ppb in AR5 (CMIP6) and 0.005 in AR6.
          The visible `[ZNTODO ...]` there is now a comment saying so
- [x] ESGF section (`methods.tex:1333-1354`): explain the frequencies and grid labels; say where full DRS details are;
      cite input4mips-validation and the input4MIPs CVs page/user guide.
    - [x] frequencies and grid labels explained, with the five combinations we provide
    - [x] say where full DRS details are: now points to the input4MIPs CVs page and user guide
          (citation placeholder, see Section 7). Checked the five frequency/grid-label combinations
          against `4010_write-input4mips-files`
    - [x] citations are in Section 7 (input4mips-validation, input4MIPs CVs)
- [x] Vertical dimension section (`methods.tex:1356-1361`) is a stub: summarise the suggested approach and give the M17 section number.
- [x] Format `numpy.linalg.svd` (and other code names) as `\texttt{}` (`methods.tex:381`).
    - [x] `numpy.linalg.svd`
    - [x] "scipy's griddata function" in general Step 2 → `\texttt{scipy.interpolate.griddata}` (checked it's what the code uses)
- [x] Check M17's names for EOFs and PCs (`methods.tex:383,385`).
- [x] Convert all remaining `ZNTODO zalue-check` / "add check" comments into real `value-check`s
      (`methods.tex:183,429,434,437,574,598,628,630,635,646,657,785,788,802,804,809,818,835,846,859,864,871,875,987,991,1064,1069`).
    - [x] all `methods.tex` ones converted (the `results.tex` ones are for Section 5).
          Constants copied from the original notebooks/config: min. 4 points for interpolation,
          100-year harmonisation (now one `HARMONISATION_TRANSITION_YEARS`, used by every harmonisation in the run),
          Mauna Loa start 1959, CH_4 Law Dome smoothing settings (from the bundle config).
          Computed from the bundle: EOF variance fractions, harmonisation offsets, CH_4 ice-core optimisation years,
          the CO_2 output matching harmonised Menking, NEEM match
    - [x] text corrected where the numbers were off:
          N_2O Menking offset "around 2.0 ppb" → 1.7 ppb (1.75);
          CO_2 Mauna Loa offset "around 1 ppm" → 1.5 ppm (1.50);
          CO_2 seasonality change "subsequent EOFs only explain 10% or less" → "each only explain around 10% or less" (second EOF 10.5%);
          NEEM match "within 0.5%" → "within 0.1%" (the original run's own tolerance; actual max 0.008%)
    - [x] the "both most northern and southern box" requirement in general Step 2 isn't an explicit rule
          (it follows from linear griddata not extrapolating and gappy months being dropped); comment added saying so

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
      and make sure it's described correctly, its latitude bounds and where the data was retrieved from (`results.tex:122-126`).
      ("Summit" ice core outlier: done, see Section 7.)
- [x] Check "IGCC is also based on NOAA data" (`results.tex:60,120,188`) and "the same source is used by IGCC" (`results.tex:244`).
      Done (2026-10-09), checked against IGCC's `notebooks/01_trace-gas-global-mean.py` (v6.4.0) and Sect. 3 of Forster et al. (2026).
      IGCC is the AR6 timeseries, extended with: NOAA's global-mean alone for CO_2 (X2019 scale, like our inputs);
      the average of NOAA and AGAGE (from 2019 on) for CH_4, N_2O, SF_6 and the halogenated gases both networks report;
      one network (some with an offset) or extrapolation (linear, or from the literature) for the rest.
      "Vimont et al. (2022)" is the ozone-depleting substances section of BAMS' State of the Climate in 2022
      (https://doi.org/10.1175/BAMS-D-23-0090.1), one of IGCC's sources for extrapolating minor gases.
    - [x] CO_2, CH_4, N_2O: `[TODO check]`s replaced by what IGCC uses for each gas
    - [x] CFC12: "the same source is used by the IGCC" was wrong. IGCC doesn't use WMO (2022),
          it averages NOAA and AGAGE, so is also a middle-ground between the two. Text now says so
    - [x] discussion (uncertainties paragraph): IGCC's reported uncertainties (0.4 ppm CO_2, 3.3 ppb CH_4, 0.4 ppb N_2O,
          carried over from AR6) look too small given the networks' lack of global coverage,
          i.e. they miss the leading (methodological) uncertainty. ZN to review the wording (ZNTODO in `discussion.tex`)
    - [ ] CH_4: the "our dataset also includes AGAGE" explanation for the ~10 ppb difference from NOAA
          applies less to IGCC, which includes AGAGE too (ZNTODO in `results.tex`)
- [x] Check IGCC citations and processing (`results.tex:388`, `src/local/historical_ghg_forcing_for_cmip7/comparison_data.py:826`).
      Citations confirmed by ZN (Section 7); `get_igcc_comparison`'s docstring now describes IGCC's sources.
- [x] Define AR6 once and cite it (`results.tex:301,357`); check "IGCC uses Hodnebrog" (`results.tex:302,358`).
      AR6 is defined on first use (Section 7). IGCC's `radeff` cites Hodnebrog et al. (2020); `[TODO check]`s removed.
- [ ] Plain-text gas names that the replacements won't catch: "CO2" (`results.tex:43`), "CH4" (`results.tex:106,132,153`).
- [ ] Heading vs text names: "CFC-12-eq"/"HFC-134a-eq" headings vs "CFC12-eq"/HFC134a in text. Check `replacements.yaml` covers all forms.
- [ ] Specify that piControl should repeat the 1850 values (`output-requirements.tex:117`; probably in results).
- [ ] Ask co-authors for other (ideally independent) comparison datasets (`results.tex:10-11`).

## 6. Abstract, introduction, output requirements, discussion, conclusion, statements

- [ ] Discussion: add the inconsistencies behind the negative-value reductions in CFC12-like/C_4F_10-like Step 5
      (pre-industrial values vs. the emissions-driven latitudinal gradient; HFC-152a's latitudinal gradient scaled down 1991-2023;
      HFC-236fa's seasonality capped every year 1996-2023; see the `ZNTODO` in CFC12-like Step 5 and the Section 4 notes)
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
- [ ] Acknowledgements: the funding text is in (CMIP IPO/ESA, ESA CCI contract).
      Still to add: the rest of the planning list in section 10 (data collectors, CMIP organisers, Forcings TT panel).
- [ ] Competing interests: add a statement about the funding being sought (`competing-interests.tex:1`).
- [ ] Results intro: user guide citation (standalone Zenodo doc to start, `results.tex:16`).

## 7. Citations

Source of truth for the input-data references: the original run's
`../CMIP-GHG-Concentration-Generation/output-bundles/v1.0.0/data/processed/dependencies.db`
(`source` table: reference text and DOI; `dependencies` table: gas to source).
It is what was written into the output files, so it has the same gaps (e.g. no PRIMAP).

Wired in (2026-10-08; build has no undefined citations):

- [x] Droste et al. → `droste_pfcs_2020` (methods, C_4F_10 appendix caption)
- [x] NOAA CO_2 / CH_4 network inputs → `lan_co2-flask_2024`, `lan_ch4-flask_2024` (flask),
      `thoning_co2-in-situ_2024`, `thoning_ch4-in-situ_2024` (in-situ), in CO_2 and CH_4 Step 1.
      These are the versions the run used; the bib's unused 2025-04-26 flask entries were replaced.
- [x] NOAA comparison timeseries in the results are NOAA's Trends files, not the flask data:
      `lan_co2-trends_2026` (global and Mauna Loa, plus `keeling_co2-exchanges_2001` for the Scripps data before May 1974)
      and `lan_ch4-n2o-sf6-trends_2026`
- [x] AGAGE → `prinn_ale-gage-agage_2000`, `prinn_agage_2018`; GAGE and ALE → `prinn_ale-gage-agage_2000`
- [x] NOAA HATS: N_2O → `dutton_hats-n2o_2022` (new); CFC12 and the other combined products → existing `dutton_hats-*_2022`;
      flask → `montzka_ods_1999`
- [x] SSP2-4.5 / M2020 → `meinshausen_ssp-ghgs_2020` (C_8F_18 section)
- [x] AR6 Ch7 SM → `IPCC_AR6-WG1-Ch-7-SM_2021` (methods equivalent species; results, where AR6 is now defined on first use)
- [x] Daniel et al. = WMO 2022 Ch. 7 → `wmo_ozone-ch7_2022` (`results.tex`, CFC12)
- [x] scipy, numpy → `virtanen_scipy_2020`, `harris_numpy_2020`
- [x] Menking et al. 2025 (in prep.) → `menking_law-dome_2025`
- [x] Law Dome - Mauna Loa merged record (Scripps' `spline_merged_ice_core_yearly.csv`)
      → `keeling_co2-exchanges_2001`, `rubino_law-dome-dataset_2019`, as its file header asks
- [x] PRIMAP-hist v2.5.1 → `gutschow_primap-hist-dataset_2024`, `gutschow_primap-hist_2016`
- [x] HadCRUT5 (5.0.2.0) → `morice_hadcrut5_2021`, `met-office_hadcrut5-data_2025`
- [x] Law Dome CH_4 → `rubino_law-dome-dataset_2019` (v3); NEEM → `rhodes-brook_neem-ch4-dataset_2019`, `rhodes_neem-ch4_2013`;
      EPICA → `epica_edml-ch4_2006`
- [x] UCI CH_4 → `simpson_ethane_2012`, retrieved from the supplement of `saunois_methane-budget_2025`;
      latitude range (approximately 71N to 47S) from UCI's network README, `blake_uci-network-readme_2010`
      (https://data.ornldaac.earthdata.nasa.gov/public/nacp/NACP_GHG_Data_Compilation/comp/README_irvinelatnet_flasks.txt;
      the "46S" previously in the `get_uci_ch4_comparison` docstring was unsourced)
      The network's individual flask samples for 2000-2008 are available (free NASA Earthdata login) at
      https://data.ornldaac.earthdata.nasa.gov/protected/nacp/NACP_GHG_Data_Compilation/data/irvinelatnet_flasks.zip
      (https://doi.org/10.3334/ORNLDAAC/1206), should we ever want to process the network ourselves.
      Not used: 2000-2008 only, and individual samples rather than a global-mean (link is also in the docstring).
- [x] CH_4 ice-core comparison → `lamantia_tropical-ch4_2026`, plus the records it compiles
      (`mitchell_ch4-constraints_2013`, `rubino_law-dome-records_2019`, `rhodes_neem-ch4_2013`)
- [x] IGCC → `forster_igcc_2026`, `smith_igcc-forcings_2026` (confirmed correct by ZN)
- [x] AR5 WG1 Ch. 8 Appendix 8.A → `IPCC_AR5-WG1-Ch-8_2013`
- [x] input4mips-validation, input4MIPs CVs → `input4mips-validation_docs_2026`, `input4mips-cvs_docs_2026` (websites, accessed 2026-10-08)
- [x] WMO 2012: not needed, `discussion.tex` now cites M17 alone
- [x] M17's sources for the non-zero pre-industrial values, cited in `tab:cfc12-like-per-gas` instead of M17
      via `PRE_INDUSTRIAL_SOURCE_CITATION_OVERRIDES` in `cfc12_like_tables.py`
      (footnote c: values are M17's estimates based on the cited sources):
      CH_3Cl, CH_3Br → `velders-daniel_ods-uncertainty_2014`; CHCl_3 → `worton_chloroform-firn_2006`, `aucott_chloroform-emissions_1999`; CH_2Cl_2 → `trudinger_halocarbons-firn_2004`;
      CF_4 → `trudinger_pfcs_2016`, `muhle_pfcs_2010`; CCl_4 (zero in the run, not M17's 0.025 ppt) → `walker_cfc-histories_2000`.
      The other zeros are M17's assumption of no natural sources (M17 Sects. 3.4, 3.5 cite nothing further), so M17 stays.
- [x] Funke et al. solar → `funke_solar_2024`; scenario in-prep paper → `nicholls_scenario-ghgs_2026`;
      Vaishali's paper → `naik_input-data-releases_2025`; forcing process over time / task team → `durack_forcing_2025`

Still to do:

- [x] BUG: "Summit" in the CH_4 ice-core supplement is Huascarán Summit Core A (9.122S, 77.605W, 6768 m; Lamantia et al., 2026),
      not Summit, Greenland. `CH4_ICE_CORE_SUPPLEMENT_SITES` in `comparison_data.py` put it at 72.58N.
      Done (2026-10-09): the site is now at 9.122S and labelled "Huascarán" (the latitude is from the paper's abstract,
      it is not in the data), the CH_4 results figure is redrawn
      and the `results.tex` text now says it is a tropical record which, before 1750, sits above our dataset
      and the polar ice cores, with the paper's interpretation (higher equatorial emissions)
      and why our dataset can't reproduce it (gradient constrained only by Law Dome and NEEM).
      We decided not to quantify the difference from our dataset (no value-check).
      The other sites' latitudes were checked too: WAIS Divide 79.47S (https://nsidc.org/data/nsidc-0493/versions/1),
      GISP2 72.58N (https://catalog.data.gov/dataset/noaa-wds-paleoclimatology-gisp2-ice-core-112kyr-methane-concentration-data;
      the code had 72.60), Law Dome 66.73S and NEEM 77.45N (match our input data), Mauna Loa 19.54N (`MAUNA_LOA_LATITUDE`).
- [x] `results.tex` CH_4 section: the UCI sentences still have two `[TODO check ...]`
      (what the UCI timeseries represents exactly; whether the sentence on its latitude range is right).
      Done (2026-10-09). The latitude range `[TODO check]` was already gone (range is from UCI's network README).
      The timeseries is UCI's estimate of the global-mean: about 80 samples from about 45 Pacific basin locations,
      collected four times a year (March, June, September, December), are averaged within each of 16 latitudinal bands
      (each holding an equal volume of air), the values of bands 15 and 16 are interpolated and the 16 band averages are averaged.
      Source: UCI's slides from the 2011 NOAA GML annual meeting
      (https://gml.noaa.gov/publications/annual_meetings/2011/slides/70-110415-A.pdf), found via a web search summary;
      Simpson et al. (2012) is paywalled, so not checked directly.
      The text now calls it "an estimate of the global-mean derived from samples collected four times a year".
      Still open in Section 5: the "UCI global-mean" figure label and anything else under the UCI item there.
- [x] `introduction.tex`: ZN to review the Durack et al. (2025) sentence and choose the other CMIP7 forcing papers
      to cite alongside `funke_solar_2024` (and an overview of model inputs over time, perhaps Durack et al., 2018, Eos, unchecked).
      Done (2026-10-09): ZN reviewed the sentence. Durack et al. (2018, Eos, "Toward Standardized Data Sets for Climate Model
      Experimentation", checked at https://doi.org/10.1029/2018EO101751) added as `durack_input4mips_2018`, cited with `durack_forcing_2025`.
      Other forcing papers: left as `funke_solar_2024` only, with ZN's note in the tex to add others (volcanic, CEDS)
      if needed or requested by reviewers.
- [x] NOAA Trends entries (`lan_co2-trends_2026`, `lan_ch4-n2o-sf6-trends_2026`): the note says "Data files created 5 September 2026";
      swap for NOAA's version string (e.g. "Version 2026-09") once confirmed on the Trends pages.
      Done (2026-10-09): both Trends pages (https://gml.noaa.gov/ccgg/trends/gl_data.html, https://gml.noaa.gov/ccgg/trends_ch4/)
      ask for "Version 2026-09", which matches our files' creation date (5 September 2026). Notes updated.
- [ ] User guide citation (`results.tex:16`, `methods.tex` output format section: `[TODO ref user guide]`), see Section 6.
      Blocked (2026-10-09): there is nothing to cite yet. The user guide is built from this repo
      (`scripts/create-historical-user-guide.sh`) but has no Zenodo record/DOI.
      Once it has one, add the bib entry and replace the four placeholders
      (`results.tex` intro, `methods.tex` output format section x2, `conclusion.tex`).
- [x] Check "IGCC uses Hodnebrog" and "the same source is used by the IGCC" (`[TODO check]`s in `results.tex`), see Section 5.
- [x] Clean-up: bib entries marked "Is this reference correct or needed?" / "rename to ..." (`references/references.bib`).
      Done (2026-10-09), no such markers left:
    - correct and cited, marker replaced by a note on the source: `prinn_ale-gage-agage_2000`,
      the five `dutton_hats-*_2022` (DOIs and versions match `dependencies.db`;
      "Chloroflurocarbon" in three titles is the spelling `dependencies.db` records, left as is)
    - fixed: `prinn_agage_2018` was the Discussions version's details with the final paper's DOI,
      now Earth Syst. Sci. Data, 10, 985-1018
    - kept, not cited yet: `rigby_methane-growth_2008`, `rigby_methane-oxidation_2017` (for the per-gas table below)
    - kept, not cited anywhere (doesn't hurt to keep them), marker removed: `hermanson_wmo-climate-update_2022`,
      `meinshausen_magicc6-part-1_2011`, `meinshausen_magicc6-part-2_2011`, `wcrp-cmip_cvs-mip-tables_2025`, `zenodo_zenodo_2025`
    - not touched: five more entries are not cited anywhere but weren't marked
      (`esgf_docs_2025`, `montreal-protocol_final-act_1989`, `unidata_netcdf_2024`, `van-vuuren_scenariomip_2025`, `IPCC_AR6-WG1-Ch-7_2021`)
- [x] Last: a table which, for each output gas, lists all the input-data references.
      Query `dependencies.db` (it includes the per-gas AGAGE papers, 46 of them, e.g. CH_4: Prinn 2018, Rigby 2008 and 2017),
      then add what it misses by hand (e.g. PRIMAP for CO_2 and CH_4).
      Decide whether the per-gas AGAGE papers are cited only in this table (the text cites Prinn et al. 2000 and 2018).
      Done (2026-10-09): new Appendix C "References for the generation of each gas" (`app:input-data-references`), Tables C1-C4
      (CO_2/CH_4/N_2O; the gases processed like CFC12, split over two tables; the gases processed like C_4F_10 and C_8F_18).
      General approach Step 1 points to it.
    - Generated by `src/local/historical_ghg_forcing_for_cmip7/input_data_references_tables.py`
      from the bundle's `data/processed/dependencies.db`, inlined at `@input-data-references-tables@` in `appendices.tex`
      (written to `tables/historical-ghg-forcing-for-cmip7/input_data_references.tex`; delete it or use `--force-rerun` to regenerate).
    - One row per gas, columns "NOAA", "AGAGE" and "Other sources" (by the source's name in `dependencies.db`).
      Every reference is listed in every row it applies to, even if every gas in the table uses it (e.g. M17).
      A column no gas in the group uses is dropped (the last table only has "Other sources").
      The first table is placed `[hbt!]`, like the first figure of the other appendices, so it follows the appendix heading.
    - The per-gas AGAGE papers are cited only in these tables (the text keeps Prinn et al. 2000 and 2018).
      53 bib entries added for them and for the records Menking et al. compile (CO_2, N_2O),
      generated from Crossref's records of the DOIs in `dependencies.db`.
      Pages for `lunt_ccl4_2018` and `prokopiou_n2o-isotopes_2018` were not in Crossref:
      Prokopiou's are from `dependencies.db`; Lunt's are left out (its reference in `dependencies.db` has none either), add them if found.
    - Added by hand (not in `dependencies.db`): PRIMAP-hist for CO_2 and CH_4; the HadCRUT5 data citation next to the paper;
      the NEEM paper next to its dataset; Rubino et al. next to Keeling et al. for the merged Law Dome - Mauna Loa record.
      Left out: the original run's placeholder for this manuscript ("Nicholls et al., 2025 (in-prep)").
    - Fixes to what `dependencies.db` records (comments in `SOURCE_BIBKEYS`): malformed DOIs for Nevison et al. (2004, 2007);
      "Roeckmann et al., 2006" is the 2003 paper; Ghosh et al. (2023) is recorded with its dataset's DOI, we cite the paper;
      "Fang et al. 2018" is 2019 in print.
    - Equivalent species have no rows: the appendix text says their input data are those of the gases they are calculated from.
    - Checked: every citation in each gas' methods section is in that gas' row, except CH_4's mention of Menking et al.
      (a comparison of method, CH_4 doesn't use their data).
    - [x] `dependencies.db` records RCMIP (Nicholls et al., 2020) for the gases processed like C_4F_10, but it isn't used, so it is left out of their rows
          (`SOURCES_RECORDED_BUT_NOT_USED`). The original run's `1404_c4f10-like_derive-latitudinal-gradient` loads the RCMIP emissions
          and regresses the latitudinal gradient PC against them, but the code which would apply the regression is commented out:
          the PC is held at zero before the Droste et al. data and linearly extrapolated after it.
          `1405_c4f10-like_create-pieces-for-gridding` doesn't use emissions.
    - [x] Same-author, same-year citations can come out of order, e.g. "Thompson et al. (2013, 2014a, c, b)",
          because `copernicus.cls` turns natbib's sorting off (`\def\NAT@sort{0}`). Noted in a comment next to
          `\bibliographystyle` in `copernicus-latex-package/template_clean.tex`. Not fixed
          (the a/b/c order is bibtex's, which the table generator doesn't know).

## 8. Nomenclature and typos

- [ ] Check hyphens and dashes throughout (all `.tex` files, plus the generated captions/tables in `src/local/`).
      We mix `-` and `--`: e.g. ranges are written both ways ("2015-2022" in the C_8F_18 summary vs "d)--f)" in methods),
      and a spaced hyphen is used as a dash ("agree in 1984 - the earliest year", "1981 - their first overlap year")
      and as a separator ("Law Dome - Mauna Loa record", "temperature - CO_2 concentration composite").
      Apply the Copernicus guidance consistently.
      It isn't in the template files in `copernicus-latex-package/`,
      so check the Copernicus manuscript preparation / English guidelines pages
      (https://publications.copernicus.org/for_authors/manuscript_preparation.html) for the rules.
      Leave maths minus signs (equations) and hyphenated compounds (e.g. "global-, annual-mean") alone.
      Consider whether `replacements.yaml` or a check in the build can enforce it.

- [ ] Check that references to equations follow the Copernicus style ("Equation" vs "Eq." vs "Eqn",
      and whether the number goes in parentheses).
      The template files in `copernicus-latex-package/` only have equation examples, not the referencing rule,
      so check the manuscript preparation / English guidelines pages (link in the item above).
      Our understanding (to verify) is "Eq. (1)" / "Eqs. (1) and (2)" in running text
      and "Equation (1)" at the start of a sentence.
      All six of ours are currently `Equation~\ref{...}` mid-sentence, without parentheses
      (`methods.tex:430,506,694,1083,1107,1268`).
      If the rule is confirmed, check whether the same applies to "Figure"/"Fig." and "Section"/"Sect."
      and consider enforcing it in the build (e.g. a macro or a check).

- [ ] Decide CFC12 vs CFC-12 in the source (abstract asks "[TODO check nomenclature]").
      `replacements.yaml` maps CFC12 → CFC-12 so the output is fine, but the source mixes both.
      Same for HFC134a vs HFC-134a.
- [x] `replacements.yaml` is applied in file order (`apply_replacements` in `scripts/compile-gmd-template-based-latex.py`),
      so a key that starts with an earlier key never matches:
      e.g. `HFC134a` is replaced before `HFC134a-eq`, `CFC12` before `CFC12-eq`, `CFC-11` before `CFC-11-eq`,
      so the `-eq` names come out with an ordinary (breakable) hyphen.
      Apply longest keys first, or reorder the file.
      Done: reordered (also `CFC-113/114/115` before `CFC-11` and `HFC-236fa` before `HFC-23`,
      which happened to render the same anyway), with a comment at the top of the file.
      The only change to the rendered text is that the `-eq` names now use `\nobreakdash`
      (The `HFC134a-eq`/`HFC134aeq` entries also had `-\nobreakdash-`, i.e. an en dash; fixed.)
- [ ] `replacements.yaml`: double check the mappings for minor species (`:4`);
      add the Halon (1202?) that is in the scenarios but not the historical data (`:31`).
- [ ] Typos: "accomodate" (`methods.tex:500,904`), "that that" (`:460`), "for for" (`:773,1025`),
      "PCA analysis" (`:773,1025,1154`), "to not source" (`:1072`), "it not all" (`:1299`), "Next we consider CH_4" missing full stop (`:722`),
      "differencs" (`results.tex:29`), "minorly" (`results.tex:371`), "radiative focing" (`conclusion.tex:16`), "foward" (`conclusion.tex:27`),
      "a key criteria" (`conclusion.tex:24`), "betwen" (`output-requirements.tex:53`), "than an another" (`output-requirements.tex:86`).

## 9. `src/local`

Manuscript-relevant:

- [x] `historical_ghg_forcing_for_cmip7/comparison_data.py:486`: check the CH_4 ice-core supplement site metadata (Summit etc.).
      Done, see the Summit/Huascarán item in Section 7.
- [x] `historical_ghg_forcing_for_cmip7/comparison_data.py:826`: check the IGCC processing (release v6.4.0).
      Done, see Section 5.
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
    - [x] GHGs are a major driver of past climate change, therefore key for model simulations. Not needed
    - [x] Updated GHG concentrations for CMIP7
    - [x] Once again a composite of multiple input sources, now covering [time range]. Not needed
    - [x] Put ERF due to GHGs in; note 'highest ever' or whatever, and the largest increases in recent times. Not needed.
    - [x] Available for modelling teams to use: "Finally, we describe where to access the data and provide a summary of key data changes compared to CMIP6." Not needed/already done
- Content
    - [ ] Mole fraction in dry air vs. mole fraction in the real atmosphere
    - [x] Different scales (ignored this time, noise in the scheme of things)
    - [ ] Historical experiment ends in 2021 even though this dataset goes to 2022 (and we hope to extend further in future).
          Go through the whole manuscript and fix this up wherever end years are mentioned:
          the output requirements say historical ends 2021-12 and the abstract says the datasets end in 2022,
          but nothing connects the two (see also the "Check the end year" item in section 6).
    - [x] Changes compared to CMIP6 (and the reasons)
    - [x] Comparisons with other datasets (lots of 'time pressure tied our hands' in here)
    - [x] If you need explanation/justification for methods, refer back to M17
    - [x] Missing halon in the historical dataset: note the ERF difference is tiny, so not ideal, but not a reason to re-write/re-run. Not a thing - only an issue for MAGICC (we don't produce it for scenarios either)
    - [x] Don't compare the seasonal cycle and latitudinal gradient from CMIP6 ESMs (unless there's way more time than expected)
    - [x] "Given the negligible radiative forcing from ..., this uncertainty does not affect the overall results."
    - [x] If time: build a portal for visualising/exploring/accessing results. Not required for this paper
- Introduction
    - [x] CMIP context
    - [x] Unique requirements of CMIP, therefore the goal of this study
- Methods
    - [x] Flowchart figure for the overall idea, with sub-panels for key variants/details (probably only 5). Done differently.
          Maybe a table too, or just text given how big the tables are in M17
          (see `tmp-figure-scribbles/hist-methodology-figure-*.jpeg` and the table decision in section 1)
    - [x] Build out based on notebooks
    - [ ] Say whether the concentrations are meant to be observed concentrations
          or concentrations that reflect background conditions (M17 discusses this; nothing in our text does yet).
    - Compared to the CMIP6 workflow
        - [x] Removed optimisation step. Not removed: we kept it for CH_4
              (the latitudinal gradient PC is optimised against ice cores, CH_4 Step 4), so there is nothing to note.
        - [x] Updated/captured polynomial smoothing (break out a package?)
        - [x] Removed N_2O interpolation from 1966-1987? (Not noted, but also fine, we note the updated ice core source)
- Results
    - [x] We don't do uncertainties (future work)
    - [x] Compare to CMIP6 throughout
    - [x] Compare to other studies here?
    - CO_2: global-mean, lat. grad. (compare with M17), seasonality
        - [x] Make plots that show all components of M17 Fig. 9, but don't put them in the paper (outreach product)
        - [x] Main paper needs a plot with the relevant components: global-mean, global-mean monthly steps, lat. gradient, seasonality
    - CH_4: global-mean, lat. grad., seasonality
        - [x] Compare with https://www.nature.com/articles/s41586-026-10938-1
    - N_2O: global-mean, lat. grad., seasonality
    - ODSs
        - [x] Reproduce/check all the pre-industrial choices, extrapolations, hard-coded zero seasonality and lat. gradient etc.
- Data format and recommendations
    - [x] Use input4MIPs CVs text
    - [ ] Point to the forcings implementation and forcing usage recording docs (i.e. how to record what f1 means), and papers as secondary
- Discussion
    - [x] Explain the differences from CMIP6 here? Done in results instead
    - Limitations (all relatively small given the use case)
        - [ ] Focussed on CMIP; don't use for inversion studies etc., you need a different product (dup from M17).
              Partly there: the output requirements say inversion studies would need our choices considered more carefully
              (`output-requirements.tex`, end of the greenhouse gas specific choices).
              Put these limitations on use in the discussion too.
        - [ ] No vertical or longitudinal resolution (dup from M17)
        - [ ] Do we want observed concs or concs that reflect background concs (dup from M17)
        - [x] Hybrid nature of calibration scales (dup from M17). Done: the fifth limitation in the discussion
        - [x] Uncertainty (dup from M17). Done: the third limitation in the discussion
- Acknowledgements
    - [ ] Everyone collecting raw data, particularly for their openness
    - [ ] CMIP organisers and all involved
    - [ ] Forcings TT panel
    - [x] Direct funding acknowledgements (ESA). Done: CMIP IPO/ESA and the ESA CCI contract are in `acknowledgements.tex`
- Data access
    - [ ] ESGF: not in `code-and-data-availability.tex` (only described in the methods' output format section)
    - [x] Zenodo. Done: the Zenodo DOI (and the GitHub repositories) are in `code-and-data-availability.tex`
