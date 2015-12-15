# mvmmm: Multi-view mixed membership models

This if the implementation in R for the posterior inference algorithm for the model presented in 

'Seppo Virtanen, Mattias Rost, Alistair Morrison, Matthew Chalmers and Mark Girolami. 
Uncovering smartphone usage patterns with multi-view mixed membership models. 
Accepted for publication in STAT, Wiley, 2015.' 

Please cite the paper if you use this software (Full citation information will appear later when available).

The folder contains source files for inferring the model and making predictions. You must install and call the R packages "Rcpp" and "Matrix" in order to run the code. The main function is mvmmm included in the modelBatch.R file, which you need to source. Compile the C++ auxiliary functions (hybridVBGibbs.cpp, hybridVBGibbsTopics.cpp and perpRcpp.cpp) using sourceCpp function included in the Rcpp package. Before calling the function you need to call opts <- getDefaults() (and pass the opts list for mvmmm). You may potentially want to modify the ouput arguments of the list. See the function mvmmm in the R file to consult the function inputs and ouputs.
