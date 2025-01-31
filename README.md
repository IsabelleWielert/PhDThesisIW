# PhDThesisIW
# Analysis of Twitching Motility 
This programme analyses the twitching motility of N. gonorrhoeae on surfaces, as described in detail by Wielert I, Kraus-Römer S, Volkmann TE, Craig L, Higgins PG, et al. (2025) in their study entitled "Pilin antigenic variants impact gonococcal lifestyle and antibiotic tolerance by modulating interbacterial forces". PLOS Biology 23(1): e3003022. https://doi.org/10.1371/journal.pbio.3003022. 

The videos can be opened via the graphical user interface GC_tracking/GC_tracking.mat, there you can directly adjust the tracking parameters and evaluate the tracks. The next step is to further process the tracks via the program Randomwalkmodel_processingoftracks.m, in which you need to read in the Data_all files of the tracking. These will be fitted to a random-walk model and the velocity and persistence time of the twitching bacteria will be determined. 

