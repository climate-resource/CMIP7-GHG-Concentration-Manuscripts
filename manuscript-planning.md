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
