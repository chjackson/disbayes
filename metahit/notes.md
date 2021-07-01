







Hi Belen,

Thanks.  It seems that the ONS gives incidence for regions, whereas we want incidence for city regions, which are contained within regions?   So the simplest thing we could do, I think, is to assume that the proportion of colorectal cancer cases that are colon cancer is constant within regions.   Then we could scale the city-region colorectal incidence counts from the GBD by this proportion.

Only trouble might be if the proportion is noisy due to small counts - we might then have to obtain this proportion over a larger area e.g. whole of England. 

There may be potential for more sophisticated statistlcal modeIling for this kind of problem if we think it is worth it - I'll have a look when I pick up the paper again. 

Chris


We have incidence for England for colon and rectum cancer [ and english region it looks like ] by 5 year age groups

* Have CRC for city regions from GBD.  Can derive counts 
* Have colon and rectum separately for broader regions [ that contain the city regions ].  Counts are published 

What about mortality, where does it come from ?  Assume same for crc from GBD 

Need to make assumption like relative prop of colon vs rectum is constant within regions

current version of the model 
hierarchical models 
confident get much value from the hierarchical models 

Chris 



Hi Belen,

Thanks for this - it makes sense - agree the GBD don't explain it very well!

I've attached a R Markdown document with code to process the GBD data.  The main things it does are

* determine "effective sample sizes" behind estimates using ci2num. 

* aggregating from local authority to city region

* disaggregating from five-year to one-year age groups

mostly implemented using dplyr.   Running ci2num is the slow part - takes about 40 minutes on my computer but could be parallelised - everything else is fast.

I think that gives a dataset that can be used for running disbayes to estimate case fatality for the city regions.

The compiled HTML document shows the code and some example calculations.   I hope this makes sense and is helpful.    I didn't completely follow everything that your code was doing - are there any other major complex computations that I missed? 

Chris


