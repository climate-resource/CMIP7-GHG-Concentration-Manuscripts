- create the table of differences from the base case
- fix all headers throughout methods
- fix intro sections for each sub-section in methods so they're consistent and explain diffs/extra information compared to base case etc.
- clean up each sub-section in methods so they don't repeat the base case more than is needed

## How to do extensions

- ingest updated network data
- calculate new global- annual-means
- optimise existing latitudinal gradient (and seasonality change for CO_2) to get best fit with (spatially interpolated?) observational network values
    - most sensible way to extend PCs
    - avoids the jump of using a different lat. gradient or seasonality over time
    - also avoids extending PCs using crude regressions
    - have to use relative seasonality calculation for gases other than CO_2, oh well just one less degree of freedom in the optimisation
- calculate full grids
- check no unexpected/spurious jumps
- process, format etc. then publish
