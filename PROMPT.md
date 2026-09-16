I would like you to make a function to plot the methods figures for the SF6-like gases.
Base this on the CH4 figure.
The key changes to make can be found in `methods.tex`.
I put my summary of the changes below too.
If you find contradictions between the two, please alert me and clarify what should be done.

Differences between CH4 figure and SF6-like figures:

1. The observational network value also needs to include any override values that are used for the global-mean
1. (the interpolation figure does not need to change, the fact we allow polar extension will be obvious in the 'fewest' inputs panel I believe)
1. The extended global-mean panel needs to show the 'obs. global-mean' as an input in the case where the obs. global-mean is overwritten by something else. The inclusion of the Trudinger et al 2016 data should be able to be handled in the same way as we handle the Menking et al data for e.g. for N2O I believe
1. The lat. gradient regression panel is a regression against total emissions, not geological emissions
1. there will only be one pc in the extended lat. gradient panel
