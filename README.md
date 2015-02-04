# clab-scalable-cca

This package hosts a collection of CCA implementations in MATLAB and takes advantage of MATLAB's
fast decomposition routines (SVD, eigendecomposition, QR, Cholesky, etc.).  

See `script_scale.m` to run the experiments with timing and correlation coefficient comparisons. 

Set the parameters of the random matrices through lines 3-6, and switch between sparse and dense
matrices by commenting in lines 9-14 (for dense matrices) and commenting out lines 17-21 or vice 
versa (for sparse matrices).  