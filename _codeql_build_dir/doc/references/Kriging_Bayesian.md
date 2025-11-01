## Bayesian framework

In the Bayesian framework, we assume that 

$\beta\sim\mathcal{N}(\beta_0,S)$

We obtain 

$\beta|Z\sim\mathcal{N}(\mu_c,\Sigma_c)$


$\Sigma_c = (X^t\Sigma^{-1}X+S^{-1})^{-1}$

and

$\mu_c=\Sigma_c(S^{-1}\beta_0+X^t\Sigma^{-1}Z)$

We obtain the Bayesian quantities:
- the estimator

$Z^{Bayes}_0 =\lambda^t_{SK}Z + (X_0 - \lambda_{SK}^tX)\mu_c$

- the variance of the estimator

$\textrm{Var}(Z^{Bayes}_0) = \textrm{Var}(Z^{SK}_0)-\lambda_{SK}^tX\Sigma_c X^t\lambda_{SK}+X_0\Sigma_c X_0^t$

- the variance of the estimation error

$\sigma_{Bayes}^2 
=\sigma_{SK}^2+(X_0-\lambda_{SK}^tX)\Sigma_c(X_0^t-X^t\lambda_{SK})
$

