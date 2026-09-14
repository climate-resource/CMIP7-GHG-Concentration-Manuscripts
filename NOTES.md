Redo the figures

I want to make a new version of the current ch4_methods and n2o_methods figure.
I want you to do this from scratch, but you are welcome to copy paste the existing code as needed.
The hard part about these figures seems to be that there are a few constraints:

1. the flying carpet has to be square
1. the maps have a fixed aspect ratio
1. the colour bars are a pain

All the other panels are flexible in terms of their aspect ratio I believe.
Hopefully this gives us the flexibility to basically place the fixed panels first, then fit everything else in around them.

There's a couple of things which I think would work best, but am ultimately flexible about

1. the observational network values panel goes in the top left
1. the figure roughly goes in order of the data processing steps
1. the flying carpet is in the bottom row (bottom right might be based to handle aspect ratio stuff, bottom left is also fine)
1. putting all the maps on the same row will probably work best, so we only pay their fixed aspect ratio price once
1. you might have to make the figure very big in order to get the number of columns you need and essentially trick matplotlib's internal margin book keeping. That is fine, we will let latex scale the figure back down in the paper, so the actual size doesn't matter that much, just the ratios between things
