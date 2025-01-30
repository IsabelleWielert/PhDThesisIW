# PhDThesisIW
# Analysis of fluorescently labelled type IV pili (T4P)
All Programs used for analysis during my PhD thesis. Here, the prorgams to analyse the pilus number, production rates, lengths, retraction velocities and elongation velocities are deposited.
These programs are developed to analyse pilus dynamics and the number of T4P when took videos of fluorescently labelled T4P (Methods in Kraus-Römer S, Wielert I, Rathmann I, Grossbach J and Maier B (2022) External Stresses Affect Gonococcal Type 4 Pilus Dynamics. Front. Microbiol. 13:839711. doi: 10.3389/fmicb.2022.839711). This should be a short discription how to use the programs and handle your data for analysis. 

To analyse the pilus number and production rate ofthe cell use the T4P_Contour program.m. This is a program, in which you need to click yourself through the analysis, in the end a file is created of a matrix, which is a kymograph of the contour intensity profile. Then, you need to use the T4PNumberProcessedAnalysis.m to determine the pilus number and production rate. 


The dynamics of the T4P is analysed via the T4P_dynamics.m program. This is a semiautomated program in which the videos are processed and the dynamics are determined in an graphical user interface (GUI). This programs creat a file in which, lengths of the pili, duration of the tracks and velocities are stored. Next, to further process the data and get the mean, you need to use the DynamicsProcessedAnalysis.m program which is collecting all .mat files from analysis before and determine the mean. 
