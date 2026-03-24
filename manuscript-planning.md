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
    - don't compare the seasonal cycle and latitudinal gradient from CMIP6 ESMs (unless I have way more time than I expected/I'm really wanting to test some CMIP data handling setup)
    - "Given the negligible radiative forcing from ..., this uncertainty does not affect the overall results."
    - get placeholder link from Jared so we can put link in the paper then build when it suits us

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

- can basically follow Malte's structure?
