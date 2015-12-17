# Prediction of locus presence/absence
This set of python, R, and BASH scripts are meant to create an absence-presence matrix (APM) from the pyRAD outputs (.loci and .stats file) for GBS data. The APM (L x L; L = # loci) is then used to make a matrix of pairwise shared loci (n x n; n = # samples), which is then used with a number of locus assembly statistics to build a predictive model. 
