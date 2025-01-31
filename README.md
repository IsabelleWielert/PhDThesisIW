# PhDThesisIW
# Analysis of fluorescently labelled type IV pili (T4P)
The programs for analysing pilus numbers, production rates, lengths, retraction and extension velocities are deposited here.
These programs were developed to analyse pilus dynamics and the number of T4P when videos of fluorescently labelled T4P were recorded (methods in Kraus-Römer S, Wielert I, Rathmann I, Grossbach J and Maier B (2022) External Stresses Affect Gonococcal Type 4 Pilus Dynamics. Front. Microbiol. 13:839711. doi: 10.3389/fmicb.2022.839711). This should be a short description of how to use the programs and how to handle your data for analysis. 


In order to analyse the pilus number and production rate of the cell, it is necessary to utilise the T4P_Contour.m program. This program requires the user to manually navigate through the analysis process. Consequently, a file is created which comprises a matrix and which is a kymograph of the contour intensity profile (.mat format). Thereafter, the T4PNumberProcessedAnalysis.m program must be employed in order to determine the pilus number and production rate. 


The dynamics of the T4P are analysed via the T4P_dynamics.m program, a semiautomated program that processes videos and determines dynamics in a graphical user interface (GUI). This program generates a file in which lengths of the pili, duration of the tracks and velocities are stored. Subsequently, to further process the data and determine the mean, the DynamicsProcessedAnalysis.m program is utilised, which collates all .mat files from the analysis and calculates the mean. 
