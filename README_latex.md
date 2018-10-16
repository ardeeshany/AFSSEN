---
output:
  html_document: default
---
# AFSSEN
**Adaptive Function-on-Scalar Smoothing Elastic Net**

AFSSEN is a methodology that simultaneously select significant predicotrs and produce smooth estimates of their parameters in a function-on-scalar linear modelwith sub-Gaussian errors and high-dimensional predictors. 

<br>
</br>

#### AFSSEN.R

We have option to control sparsity and smoothness separately with using two penalty parameters $\lambda_H$ and $\lambda_K$. We aim to estimate a smooth version of $\bf{\beta}$ to minimize the following target function.

<br>
</br>

<center>
![equation](http://latex.codecogs.com/gif.latex?L_%7B%5Clambda%7D%28%5Cbeta%29%3D%5Cdfrac%7B1%7D%7B2N%7D%20%5C%7C%5Cbf%7BY%7D-%5Cbf%7BX%7D%5Cbf%7B%5Cbeta%7D%5C%7C_%7B%5Cmathbb%7BH%7D%7D%5E2&plus;%5Cfrac%7B%5Clambda_%7BK%7D%7D%7B2%7D%20%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5C%7CL%28%5Cbeta_%7Bi%7D%29%5C%7C_%7B%5Cmathbb%7BK%7D%7D%5E2%20&plus;%20%5Clambda_%7BH%7D%20%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5Ctilde%7Bw%7D_%7Bi%7D%20%5C%7C%20%5Cbeta_%7Bi%7D%5C%7C_%7B%5Cmathbb%7BH%7D%7D)
</center
>
<br>
</br>


The following **AFFSEN()** function helps us to estimate the smooth $\bf{\beta}$ and find the significant predictors:

<center>
![](AFSSEN.png)
</center>

<br>
</br>

[Here](Readme.pdf) is the AFSSEN function details.

and here is the functionality of each parameters:

<center>
<p style="text-align: center;">
| Parameter            | functionality              | Details                |
| ---------------------| -------------------------- |----------------------- |
| `X`             |Numerical design matrix |It should be a N*I matrix (N= #observations ; I =#predictors) |
`Y`               |Matrix of pointwise evaluation for observations on T_domain| It should be a N*m matrix where|
`T_domain`|Time domain for evaluation of Y and generating kernel|Default : T_domain = seq(0,1,m=50)|
`type_kernel`|Type of kernel|‘exponential’, ‘gaussian’, ‘sobolev’|
`param_kernel`|Kernel parameter|In all types, the time domain is seq(0,1,50)|
`thres`|<span style="color:red"> Stopping criteria:</span> $\bf{\beta}$ increment threshold| $\| \bf{\beta}^T -  \bf{\beta}^{T-1}\|_{\mathbb{H}} < thres$|
`number_non_zeros` |<span style="color:red"> Stopping Criteria:</span> **Kill switch**; number of nonzero predictors||
`ratio_lambda_H` | $\dfrac{\lambda_{Hmax}}{\lambda_{Hmin}}$||
`number_lambda_H` |Generate number of log-equally spaced in $[\lambda_{Hmin},\lambda_{Hmax}]$||
`num_lambda_H_NONad`|Number of $\lambda_H$ in non-adaptive step||
`lambda_H` |You have option to insert a vector of $\lambda_H$ | If you want to make the log-equally spaced by `ratio_lambda_H` and `number_lambda_H`, set **lambda_H=numeric()** |
`lambda_K` | Vector of $\lambda_K$||
`early_CV` |0 or 1 : applying the “early_CV_thres” stopping criteria or not||
`early_CV_thres` | <span style="color:red"> Stopping Criteria: </span> Breaking point in CV plot | $\dfrac{\| CV(h-1,k) - CV(h,k)\|}{CV(h-1,k)} < early\_CV\_thres$|
`max_ite_nadp`| <span style="color:red"> Stopping Criteria: </span> Maximum iteration of coordinate descent alg. in non-adaptive step||
`max_ite_adp`| **Stopping Criteria:** Maximum iteration of coordinate descent alg. in adaptive step||
`max_ite_final` | **Stopping Criteria:** Maximum iteration of coordinate descent algorithm for the optimum $\lambda_H$ and $\lambda_K$||
`target_inc` | <span style="color:red"> Stopping Criteria: </span> 0 or 1 : if target function is increased, stop!||
`proportion_training_set`| Proportion of training set for estimation in non-adaptive step ||
`fold_ad` | Number of fold for using CV in adaptive steps to find optimum $\lambda_H$ and $\lambda_K$ and then estimation ||
</p>
</center>







