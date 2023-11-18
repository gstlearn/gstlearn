# Expectation, Variance, Covariance

## Expectation

### Definition

If $X$ is a random variable with output space $H$ and with  density $f(x)$.

Then, $$E[X]=\int_H x f(x) dx$$

### Linearity

#### Sum of random variables

If $X$ and $Y$ are two random variables with respective outputs spaces $H$ and $K$, with respective densities $f$ and $g$, with a joint density $f(x,y)$ and with a finite expectation. 

We have $$E[X+Y]=E[X]+E[Y]$$

##### Idea of the proof

Let consider the function $q$ defined by $q(x,y)=x+y$

$E[q(X,Y)] = \int_{H}\int_{K} q(x,y) f(x,y)dx dy$

$\phantom{E[Z]} = \int_{H}\int_{K} (x+y) f(x,y)dx dy$

$\phantom{E[Z]} = \int_{H}\int_{K} x f(x,y)dx dy + \int_{H}\int_{K} y f(x,y)dx dy$

$\phantom{E[Z]} = \int_{H} x\int_{K} f(x,y)dy dx + \int_{K} y\int_{H}  f(x,y)dx dy$

$\phantom{E[Z]} = \int_{H} xf(x) dx + \int_{K} y(y)g(y)dy$

$\phantom{E[Z]} = E[X]+E[Y]$

#### Product by a constant

We also have

$$E[aX]=aE[X]$$

### Positivity

If $X$ is a positive random variable i.e $$P(X\geq 0)=1$$ then 

$$E[X]\geq 0$$

Where $X$ is a positive random variable, if $E[X]=0$, then 

$$P(X=0)=1$$

### Constant

$$E[a] = a$$

## Covariance and Variance

### Definition

$$\textrm{Cov}(X,Y) = E[(X-E[X])(Y-E[Y])]$$

### Properties 

#### Symmetry

$$\textrm{Cov}(X,Y)=\textrm{Cov}(Y,X)$$

##### Other expression

Sometimes more convenient

$$\textrm{Cov}(X,Y) = E[XY] - E[X]E[Y]$$

In particular

$$\textrm{Var}(X)  = \textrm{Cov}(X,X) = E[X^2] - E[X]^2$$

#### Linearity

$$\textrm{Cov}(aX+bY,Z) = a\textrm{Cov}(X,Z)+b\textrm{Cov}(Y,Z)$$

#### Variance and covariance

$$\textrm{Cov}(X,X) = E[(X-E[X])(X-E[X])]=\textrm{Var}(X)$$

Note that the variance is always positive (as the expectation of a square of a random variable).

#### Covariance between a variable and a constant

$$\textrm{Cov}(X,a) = 0$$

Consequence :

$$\textrm{Var}(a)=0$$

Reciprocally, if a random variable has a variance equal to 0, then the variable is constant.

#### Variance of a linear combination

$$\textrm{Var}\left(\sum_{i=1}^n\lambda_i Z_i\right) = \sum_{i=1}^n\sum_{j=1}^n \lambda_i\lambda_j\textrm{Cov}(Z_i,Z_j)$$

##### Applications

$$\textrm{Var}(aX) = a^2 \textrm{Var}(X)$$

$$\textrm{Var}(aX+bY)=a^2\textrm{Cov}(X,X)+2ab\textrm{Cov}(X,Y)+b^2\textrm{Cov}(Y,Y)$$

$$\textrm{Var}(aX-bY)=a^2\textrm{Cov}(X,X)-2ab\textrm{Cov}(X,Y)+b^2\textrm{Cov}(Y,Y)$$

$$\textrm{Var}(X+a) = \textrm{Var}(X)$$


### Covariance matrix

When we have a set of random variables $Z_1,\dots,Z_n$.

For each pair $(k,l)$, if we denote $$c_{kl} = \textrm{Cov}(Z_k,Z_l)$$

We can store the $c_{kl}$'s in a matrix $$\Sigma = \left[\begin{array}{ccc}c_{11} &\dots & c_{1n}\\
c_{21} & \dots & c_{2n}\\
\vdots & \ddots & \vdots\\
c_{n1} & \dots & c_{nn}\end{array}\right]$$

$\Sigma$ is named the covariance matrix of the random vector $$Z=\left[\begin{array}{c}Z_1\\ \vdots\\ Z_n\end{array}\right]$$

Note that we can rewrite 

$$\textrm{Var}\left(\sum_{i=1}^n\lambda_iZ_i\right) = \lambda^T \Sigma \lambda$$

where $$\lambda =\left[\begin{array}{c}\lambda_1\\ \vdots\\\lambda_n\end{array}\right]$$

and $^T$ designates the transposition

$$\lambda^T =\left[\begin{array}{ccc}\lambda_1& \dots & \lambda_n\end{array}\right]$$

Since a variance is always positive, the variance of any linear combination as to be positive. Therefore, a covariance matrix is always (semi-)positive definite, i.e

For each $\lambda$ $$\lambda^T \Sigma \lambda\geq 0$$

#### Cross-covariance matrix

Let consider two random vectors $X=(X_1,\dots,X_n)$ and $Y=(Y_1,\dots,Y_p)$.

We can consider the cross-covariance matrix $\textrm{Cov}(X,Y)$ where element corresponding to the row $i$ and the column $j$ is $\textrm{Cov}(X_i,Y_j)$

If $A$ and $B$ are some matrices (of constants)

$$\textrm{Cov}(AX,BY) = A\textrm{Cov}(X,Y)B^T$$

#### Exercise

Suppose that we want to estimate a quantity modeled by a random variable $Z_0$ as a linear combination of known quanties
$Z_1,\dots, Z_n$ stored in a vector $$Z=\left[\begin{array}{c}Z_1\\ \vdots\\ Z_n\end{array}\right]$$

We will denote $$Z_0^\star = \sum_{i=1}^n \lambda_i Z_i = \lambda^T Z$$ this (random) estimator.

We know the covariance matrix of the full vector $(Z_0,Z_1,\dots,Z_n)$ that we write with blocks for convenience:

$$\left[\begin{array}{cc}\sigma_0^2 & c_0^T \\
c_0 & C\end{array}\right]$$


where 

* $\sigma^2_0 = \textrm{Var}(Z_0)$
* $c_0 = \textrm{Cov}(Z,Z_0)$
* $C$ is the covariance matrix of $Z$.

Compute the variance of the error $$Z_0^\star-Z_0$$

##### Solution

$\textrm{Var}(Z_0^\star-Z_0) = \textrm{Cov}(Z_0^\star-Z_0,Z_0^\star-Z_0)$

$\phantom{\textrm{Var}(Z_0^\star-Z_0)} = \textrm{Var}(Z_0) -2 \textrm{Cov}(Z_0^\star,Z_0) + \textrm{Var}(Z_0)$

$\phantom{\textrm{Var}(Z_0^\star-Z_0)} = \textrm{Var}(\lambda^TZ) -2 \textrm{Cov}(\lambda^T Z,Z_0) + \sigma_0^2$

$\phantom{\textrm{Var}(Z_0^\star-Z_0)} = \lambda^T\textrm{Var}(Z)\lambda -2 \lambda^T\textrm{Cov}( Z,Z_0) + \sigma_0^2$

$\phantom{\textrm{Var}(Z_0^\star-Z_0)} = \lambda^TC\lambda -2 \lambda^Tc_0 + \sigma_0^2$

### Correlation coefficient

The covariance is a measure of the link between two variables. However it depends on the scale of each variable. To have a similar measure which is invariant by rescaling, we can use the correlation coefficient:

$$\rho(X,Y)=\frac{\textrm{Cov}(X,Y)}{\sqrt{\textrm{Var}(X)\textrm{Var}(Y)}}$$

When the correlation coefficient is equal to $1$ or $-1$, we have

$$Y=aX+b$$ 

with 

* $a>0$ if $\rho(X,Y)=1$ 
* $a<0$ if $\rho(X,Y)=-1$

Note that $\rho(X,Y)$ can be equal to $0$ even if the variables are strongly linked.

The usual example is a variable $X$ with a pair density ($f(-x)=f(x)$) and $Y=X^2$:

$$\textrm{Cov}(X,Y)=\textrm{Cov}(X,X^2)=E[X^3]-E[X]E[X^2]=E[X^3]=\int_{\mathbb{R}} x^3f(x)dx =0$$

