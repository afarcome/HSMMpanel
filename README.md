
Code for "Exact score and information matrix for panel hidden semi-Markov models"

This project reports R code implementing the methods  in Farcomeni, A. "Exact score and information matrix for panel hidden semi-Markov models"

The following libraries shall be installed in R before using code from this project:

library(snipEM)
library(parallel)
library(snow)
library(numDeriv)
library(mhsmm)
library(compiler)

The R package snipEM is not currently available on CRAN, a source .tar.gz and .zip with binaries is attached.

While the model is very general, this sample code focuses on two cases: two binary outcomes, and two continuous outcomes.
Also the sojourn distribution is fixed as a discrete gamma. 

Code for the case of two binary outcomes is in file _functionsBin.R; while
code for the case of two Gaussian outcomes is in file _functionsCont.R.

File exampleBin.R provides an example with two binary outcomes.
File exampleCont.R provides an example with two continuous outcomes.

Since examples can be time consuming, results can be found in workspaces,
namely resBin.RData for the case of two binary outcomes; and resCont.RData for the case of two continuous outcomes. 

The main functions are called estHSMMpanelBin and estHSMMpanelCont.

In both cases inputs are:

Y: n by T by p array, with n sample size, T time replicates, p number of outcomes 
X: n by T by q array, with q number of predictors 
k: number of latent masses 
cl: output of function makeCluster, as the function uses parallel computing 
init:
Pinit: initial values for the transition matrix (it must be matrix(c(0,1,1,0),2,2) when k=2)
se: (default: FALSE)

Outputs are:

pars: the MLE (note that some parameters are on the log or logit scale)  
lik: the value of the log-likelihood at convergence 
score: the score (NA if se=FALSE)
info: the information matrix (NA if se=FALSE)
se: estimated standard errors (NA if se=FALSE)



