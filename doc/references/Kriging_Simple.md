# Kriging 

Let suppose that 

$Z(s_i) = X_i\beta + Y(s_i)$

where $Y$ is a second order stationary random field with mean 0 and covariance function $C$.

If $Z$ is a vector of observations, we denote 
$Z = X\beta + Y$ with $\Sigma$ the covariance of Y

## Simple Kriging 
If $\beta$ is known, we can obtain the simple kriging 

$Z_0^{SK} = X_0\beta + \Sigma_0^t\Sigma^{-1}(Z-X\beta) = X_0\beta + \lambda_{SK}^t(Z-X\beta)$

with:

- the simple kriging weights

$\lambda_{SK}=\Sigma^{-1}\Sigma_0$ 

- the variance of the estimator

$\textrm{Var}(Z_0^{SK})=\lambda_{SK}^t\Sigma\lambda_{SK}=\lambda_{SK}^t\Sigma_0$

- the estimation variance

$\sigma_{SK}^2 = \textrm{Var}(Z_0-Z_0^{SK}) = \sigma_0^2-\Sigma_0^t\Sigma^{-1}\Sigma_0=\sigma_0^2-\lambda_{SK}^t\Sigma_0$

### In matrix notation:

Simple Kriging System

$
      \begin{bmatrix}
	\Sigma
      \end{bmatrix}
      \times
      \begin{bmatrix}
	\lambda_{SK}
      \end{bmatrix}
      =
      \begin{bmatrix}
        \Sigma_0
      \end{bmatrix}
$

Estimation

$    Z_0^{SK} =
     \begin{bmatrix}
	Z
     \end{bmatrix}^t
     \times
     \begin{bmatrix}
	\lambda_{SK}
     \end{bmatrix}
     + m ({ 1 - \sum{\lambda_{SK}}} )
$

Variance of Estimation error

$
   \sigma_{SK}^2 = \sigma_0^2 -
   \begin{bmatrix}
     \lambda_{SK}
   \end{bmatrix}^t
   \times
   \begin{bmatrix}
     \Sigma_0
   \end{bmatrix}
$

Variance of Estimator

$
   \textrm{Var}(Z_0^{SK}) =
   \begin{bmatrix}
     \lambda_{SK}
   \end{bmatrix}^t
   \times
   \begin{bmatrix}
     \Sigma
   \end{bmatrix}
   \times
   \begin{bmatrix}
     \lambda_{SK}
   \end{bmatrix} =
   \begin{bmatrix}
     \lambda_{SK}
   \end{bmatrix}^t
   \times
   \begin{bmatrix}
     \Sigma_0
   \end{bmatrix}
$
