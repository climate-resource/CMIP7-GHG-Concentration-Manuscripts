Historical manuscript:

Notes from reading Meinshausen at al 2017 (M17): https://gmd.copernicus.org/articles/10/2057/2017/gmd-10-2057-2017.pdf

- invite all the data producers to be co-authors
- abstract
    - GHGs major driver of past climate change, therefore key for model simulations
    - updated GHG concentrations for CMIP7
    - once again a composite of multiple input sources, now covering [time range]
    - put ERF due to GHGs in, note 'highest ever' or whatever and largest increases in recent times
    - available for modelling teams to use: "Finally, we describe where to access the data and provide a summary of key data changes compared to CMIP6."

- content
    - mole fraction in dry air vs. mole fraction in real atmopshere hooha/BS
    - different scales
    - historical experiment to end in 2021 even though this dataset goes to 2022 (and we hope to extend further in future)
    - changes compared to CMIP6 (and reasons therefore)
    - comparisons with other datasets
        - lots of 'time pressure tied our hands' in here
    - if you need explanation/justification for methods, refer back to Meinshausen et al. 2017
    - missing halon in historical dataset
        - note ERF difference as tiny, therefore not ideal, but not a reason to re-write/re-run
    - don't compare the seasonal cycle and latitudinal gradient from CMIP6 ESMs (unless I have way more time than I expected/I'm really wanting to test some CMIP data handling setup)
    - "Given the negligible radiative forcing from ..., this uncertainty does not affect the overall results."
    - get placeholder link for providing results from Jared so we can put link in the paper then build when it suits us

    - introduction
        - CMIP context
        - unique requirements of CMIP therefore goal of this study
    - methods
        - flowchart figure for overall idea. sub-panels for key variants/details (I think there's only 5). Maybe have a table too (although maybe just stick to text given how big the tables are in the 2017 paper)
        - build out based on notebooks
        - compared to CMIP6 workflow
            - removed optimisation step
            - updated/captured polynomial smoothing (TODO: break out package ?)
            - removed N2O interpolation from 1966-1987?
    - results
        - we don't do uncertainties (future work)
        - compare to CMIP6 throughout here
        - compare to other studies here ?
        - CO2
            - global-mean
            - lat. grad. (compare with M17)
            - seasonality
            - make plots that show all components of M17 Fig. 9, but don't put such figures in the paper (they can go in some outreach product).
            - wait and see what the text needs before deciding what plots to put in
                - main paper will need plot with the relevant components: global-mean, global-mean monthly steps, lat. gradient, seasonality
        - CH4
            - global-mean
            - lat. grad.
            - seasonality
        - N2O
            - global-mean
            - lat. grad.
            - seasonality
        - ODSs
            - have to reproduce/check all the pre-industrial choices, extrapolations, hard-coding zero seasonality and lat. gradient etc.
    - data format and recommendations
        - use input4MIPs CVs text
        - point to forcings implementation and forcing usage recording (i.e. how to record what f1 means) docs (and papers, but as secondary)
    - discussion
        - explain difference from CMIP6 here ?
        - limitations (all relatively small given use case)
            - focussed on CMIP, don't use for inversion studies etc., you need a different product (dup from M17)
            - no vertical or longitudinal resolution (dup from M17)
            - do we want observed concs or concs that reflect background concs (dup from M17)
            - calibration scales hybrid nature (dup from M17)
            - uncertainty (dup from M17)
    - acknowledgements
        - everyone collecting raw data, particularly openness
        - CMIP organisers and all involved
        - Forcings TT panel
        - direct funding acknowledgements (ESA)

Scenarios manuscript:

Notes from reading Meinshausen at al 2020 (M20): https://gmd.copernicus.org/articles/13/3571/2020/gmd-13-3571-2020.pdf

- invite all input data people to be co-authors
- abstract
    - ScenarioMIP overarching exercise for scenario projections, exploring range of future climates to inform mitigation decisions and adaptation planning
    - Here we provide the GHG concentrations to use in both the 'main' experiment and their extensions, based on a smooth harmonisation with the historical timeseries
    - some key numbers: range of CO2 concentrations and CH4 concentrations ? GHG ERF range too ?
    - Speed required, which means tradeoffs: "Results will inevitably differ from ESM outputs for same input emissions, but the need for speed means that we cannot use ESMs to produce these boundary conditions instead."
    - available for modelling teams to use.

- introduction
    - ScenarioMIP exists, part of wider infrastructure/effort
    - provide data for the ScenarioMIP scenarios, which range from X to Y via (one para description of scenarios, then reference overview paper for the rest of the info, tier information etc.)
    - in terms of technical details/requirements, within ScenarioMIP, range of climate models participate. Increasingly CO2 and CH4 emissions-driven capability, but still requirement for consistent GHG concentrations for other spcecies and to support models without this capability
    - want consistency with underlying emissions/scenarios (to facilitate later cross WG analysis), hence need to produce concentrations with a model that is emissions-driven. Model also needs to be fast, hence use SCMs, specifically MAGICC for its demonstrated ESM-emulation ability and use in IPCC reports.
    - also require harmonisation with historical dataset to avoid spurious jumps in the transition from history to scenario. Build on historical data.

- methods [break down more like historical. Describe general idea, then specific implementation for various gases]
    - MAGICC (same model used for CMIP5 and CMIP6). In this case, we use an interim version of MAGICC, MAGICC v7.6.0a3. Given time constraints etc., provides estimate aligned with AR6, but also including key update to methane natural feedbacks (see Trevor's paper).
    - WMO 2022 and Western for some sources
    - one box model for some gases for clarity and simplicity and consistency with underlying lifetime assessments
    - roll straight off the end of historical, to ensure consistency with scenario emissions (which are harmonised to 2023 values) [find a better way to word and explain this]
    - lean on M20 for methods re extending seasonality and lat. gradients
    - refer to other papers for emissions harmonisation and infilling (use emissions-driven for everything for running MAGICC, slight inconsistency from prescribing the concentrations that are generated in other ways, but ok, minor and these are sensible boundary conditions, not the only answer)
    - extensions based on extended emissions (refer to other paper, again slight inconsistency as above, ok)
    - harmonisation of concentrations: gradient preserving

- results
    - general observations (?) e.g. where co2 is up/down and ranges, same for ch4, n2o, ODSs, other f-gases
    - radiative forcing
    - don't repeat tables of output
    - comparison with CMIP6

- discussion
    - reasons for differences with CMIP6
    - don't discuss recent obs, save that for extension paper
    - changes in context of long-term records
    - no ESM simluations for lat. gradients and seasonality
    - harmonisation still not perfect

- limitations (can also blend with discussion)
    - limited number of scenarios, relative to 'full range' (use non-markers and maybe AR6 db for 'full range')
    - sequential scenario generation process i.e. best-estimate after CMIP7 will certainly differ from the 'pre CMIP7 ESM runs' we provide here. Can't be avoided, don't @ us

- climate paper abstract
    - range of warming (and ERF ?) outcomes and comparison to CMIP6
