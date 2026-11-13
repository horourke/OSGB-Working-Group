# 041_modvar


The current directory contains the functions used to estimate and select tuning parameters for the MOD-VAR procedure. Further descriptions are provided in the comments of the functions.

- <code>./cv.modvar.R</code>: file that contains the function  <em>cv.modvar()</em> used to perform our proposed procedure via cross validation. 

- <code>./bic.modvar.R</code>: file that contains the function  <em>bic.modvar()</em> used to perform our proposed procedure via Bayesian information criteria (BIC). 

- <code>./ada.modvar.R</code>: file that contains the function  <em>ada.modvar()</em> used to perform our proposed procedure with an added adaptive layer. 

- <code>./matrix_fista.cpp/</code>: functions that perform the optimization, written in C++. Rcpp is used to incorporate them into our R environment.

- <code>./auxfunct.R</code>: contains the code for generating the cross-validation sets (by-subject or rolling window), and additional functions for performing model selection via BIC.


- <code>./adaweights.R</code>: function that calculates the adaptive weights to be used by our proposed adaptive MOD-VAR procedure.
