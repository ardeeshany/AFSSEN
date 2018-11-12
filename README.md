---
output:
  pdf_document: default
---
# AFSSEN
**Adaptive Function-on-Scalar Smoothing Elastic Net**

AFSSEN is a methodology that simultaneously select significant predicotrs and produce smooth estimates of their parameters in a function-on-scalar linear modelwith sub-Gaussian errors and high-dimensional predictors. 

<br>
</br>

### Documentation

For installing this package, use

`
devtools::install_github("ardeeshany/AFSSEN")
`
<br>
</br>

#### AFSSEN.R

We have option to control sparsity and smoothness separately with using two penalty parameters $\lambda_H$ and $\lambda_K$. We aim to estimate a smooth version of $\bf{\beta}$ to minimize the following target function.

<br>
</br>

<center>
![equation](http://latex.codecogs.com/gif.latex?L_%7B%5Clambda%7D%28%5Cbeta%29%3D%5Cdfrac%7B1%7D%7B2N%7D%20%5C%7C%5Cbf%7BY%7D-%5Cbf%7BX%7D%5Cbf%7B%5Cbeta%7D%5C%7C_%7B%5Cmathbb%7BH%7D%7D%5E2&plus;%5Cfrac%7B%5Clambda_%7BK%7D%7D%7B2%7D%20%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5C%7CL%28%5Cbeta_%7Bi%7D%29%5C%7C_%7B%5Cmathbb%7BK%7D%7D%5E2%20&plus;%20%5Clambda_%7BH%7D%20%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5Ctilde%7Bw%7D_%7Bi%7D%20%5C%7C%20%5Cbeta_%7Bi%7D%5C%7C_%7B%5Cmathbb%7BH%7D%7D)
</center>

<br>
</br>


The following **AFFSEN()** function helps us to estimate the smooth $\bf{\beta}$ and find the significant predictors:

<center>
![](https://github.com/ardeeshany/AFSSEN/blob/master/inst/doc/AFSSEN.png)
</center>


