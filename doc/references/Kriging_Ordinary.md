# Kriging 

Let suppose that 

$Z(s_i) = X_i\beta + Y(s_i)$

where $Y$ is a second order stationary random field with mean 0 and covariance function $C$.

If $Z$ is a vector of observations, we denote 
$Z = X\beta + Y$ with $\Sigma$ the covariance of Y


## Universal kriging

If $\beta$ is unknown, we can estimate it by 

$\hat\beta =  \Sigma_c X^t\Sigma^{-1}Z$ 

Introducing the notation

$\Sigma_c =  (X^t\Sigma^{-1}X)^{-1} $

then

$\hat\beta = \Sigma_c X^t\Sigma^{-1}Z$ 

$\textrm{Var}(\hat\beta)=\Sigma_c$

The Universal kriging is obtained by first computing $\hat\beta$ and then pluging $\hat\beta$  in the simple kriging procedure.

$Z^{UK}_0 = X_0\hat\beta + \Sigma_0^t\Sigma^{-1}(Z-X\hat\beta)= \Sigma_0^t\Sigma^{-1}Z + (X_0 - \Sigma_0^t\Sigma^{-1}X)\hat\beta$

We can rewrite everything with respect to $Z$

$Z^{UK}_0 =  (\Sigma_0^t\Sigma^{-1} + (X_0 - \Sigma_0^t\Sigma^{-1}X)\Sigma_c X^t\Sigma^{-1})Z \\
=(\lambda_{SK}^t+(X_0-\lambda_{SK}^tX) \Sigma_c X^t\Sigma^{-1})Z\\
=\lambda_{UK}^tZ$ 

with

- the Universal Kriging Weights

$\lambda_{UK}=\lambda_{SK}+\Sigma^{-1}X \Sigma_c(X_0^t-X^t\lambda_{SK})$

- the Lagrange coefficients

$\mu_{UK}=\Sigma_c (X_0 - \lambda_{SK}^tX)^t$

- the variance of the estimator is

$\textrm{Var}(Z^{UK}_0) = \lambda_{UK}^t\Sigma\lambda_{UK} \\
=\textrm{Var}(Z^{SK}_0) +2\lambda_{SK}^tX \Sigma_c \Sigma_c (X_0^t-X^t\lambda_{SK})+(X_0-\lambda_{SK}^tX)\Sigma_c X^t\Sigma^{-1}X\Sigma_c (X_0^t-X^t\lambda_{SK})\\
=\textrm{Var}(Z^{SK}_0) +2\lambda_{SK}^tX\Sigma_c (X_0^t-X^t\lambda_{SK})+(X_0-\lambda_{SK}^tX)\Sigma_c (X_0^t-X^t\lambda_{SK})\\
=\textrm{Var}(Z^{SK}_0)+(\lambda_{SK}^tX+X_0)\Sigma_c (X_0^t-X^t\lambda_{SK})\\
=\textrm{Var}(Z^{SK}_0)-\lambda_{SK}^tX\Sigma_c X^t\lambda_{SK}+X_0 \Sigma_c X_0^t$

- the estimation variance

$\sigma_{UK}^2 = \sigma_0^2 - 2\textrm{Cov}(Z_0,Z^{UK}_0)+ \textrm{Var}(Z^{UK}_0)\\
= \sigma_0^2 -2\Sigma_0^t\lambda_{UK}+\textrm{Var}(Z^{UK}_0)\\
= \sigma_0^2 -2\Sigma_0^t(\lambda_{SK}+\Sigma^{-1}X \Sigma_c(X_0^t-X^t\lambda_{SK}))+\textrm{Var}(Z^{SK}_0)-\lambda_{SK}^tX \Sigma_c X^t\lambda_{SK}+X_0 \Sigma_c X_0^t\\
=  \sigma_0^2 -\Sigma_0^t\lambda_{SK} -2\Sigma_0^t\Sigma^{-1}X \Sigma_c (X_0^t-X^t\lambda_{SK})-\lambda_{SK}^tX \Sigma_c X^t\lambda_{SK}+X_0 \Sigma_c X_0^t\\
=\sigma_{SK}^2-2\lambda_{SK}^tX \Sigma_c (X_0^t-X^t\lambda_{SK})-\lambda_{SK}^tX \Sigma_c X^t\lambda_{SK}+X_0 \Sigma_c X_0^t\\
=\sigma_{SK}^2+(X_0-\lambda_{SK}^tX) \Sigma_c (X_0^t-X^t\lambda_{SK})
$

### In matrix notation:

Universal Kriging System

$
      \begin{bmatrix}
	\Sigma & X \\
         X^t   & 0
      \end{bmatrix}
      \times
      \begin{bmatrix}
	\lambda_{UK} \\
	-\mu
      \end{bmatrix}
      =
      \begin{bmatrix}
        \Sigma_0 \\
	X_0^t
      \end{bmatrix}
$

Estimation

$    Z_0^{UK} =
     \begin{bmatrix}
	Z \\
	0
     \end{bmatrix}^t
     \times
     \begin{bmatrix}
	\lambda_{UK} \\
	-\mu
     \end{bmatrix}
$

Variance of estimation error

$
   \sigma_{UK}^2 = \sigma_0^2 -
   \begin{bmatrix}
     \lambda_{UK} \\
     -\mu
   \end{bmatrix}^t
   \times
   \begin{bmatrix}
     \Sigma_0 \\
     X_0^t
   \end{bmatrix}
$

Variance of estimator

$
   \textrm{Var}(Z^{UK}_0) =
     \begin{bmatrix}
     \lambda_{UK}
   \end{bmatrix}^t
   \times
   \begin{bmatrix}
     \Sigma
   \end{bmatrix}
   \times
   \begin{bmatrix}
     \lambda_{UK}
   \end{bmatrix}
$
