I would like you to make a function to plot the methods figures for the C4F10-like gases.
Base this on the CH4 figure.
The key changes to make can be found in `methods.tex`.
I put my summary of the changes below too.
If you find contradictions between the two, please alert me and clarify what should be done.

Differences between CH4 figure and SF6-like figures:

1. The observational network panel needs to be relabelled as inputs and just plot the input timeseries, one line for each latitude in the input
1. Get rid of the obs. counts and obs. locations and interpolation panels
1. make panel c be the assumed latitudinal gradient
1. panel d is then the derived lat. gradient PC
1. panel e is then the derived global-mean
1. The global-mean and lat. gradient extensions are then both: Droste et al based, assumed constant before the start of Droste and linear extrapolation after the end of Droste

To make this work, you might need to re-run notebooks or grab output. Hopefully it is clear where/how to do this. If you have any doubts, please ask me.

Once you have done this, then please do the plot for C8F18.
As all we do here is use CMIP6 output, all you need to show is:

1. extended global-mean
1. extened lat. gradient PC
1. Lat. gradient EOF
1. Seasonality
1. monthly spatial means, yearly spatial means and native resolution (same as for other gases)
