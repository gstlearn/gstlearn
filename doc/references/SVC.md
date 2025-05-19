The variable of interest $\{z(s), s \in \mathcal{D} \}$ is interpreted
as the realization of some non stationary random function, $Z(s)$ with $Z(s) \in \mathbb{R}$ and $\mathcal{D} \subset \mathbb{R}^d$ for $d \leq 3$.

### The model of kriging with external drifts

We consider first the case where the coefficients are constant but unknown: $$
Z(s) = \sum_{l = 0}^{L} a_l \, f^l(s) + Y(s).
$$ The coefficients of the spatial regression $a_l$ are unknown and the residual $Y$ is a centred stationary random function whose covariance function is $C$.

The two first moments of $Z$ are therefore:

-   $m(s) = \mathbb{E}\{Z(s)\} = \sum_{l = 0}^{L} a_l \, f^l(s)$

-   $\text{Cov}(Z(s), Z(s')) = \mathbb{E}\{(Z(s)- m(s))(Z(s')- m(s'))\} = C(s - s')$

The Best Linear Unbiased Predictor (BLUP) of $Z$ at the target point $\mathbf{o}$ is kriging with external drift usually considered in geostatistics

$$\hat{Z}(\mathbf{o}) = \sum_{\alpha \in A} \lambda^{\alpha}Z(\alpha)$$

with the kriging weights obtained minimizing the variance of the estimation error with the unbiasedness conditions. 

The observations used to predict the variable $Z$ at the target $\mathbf{o}$ may be altered by some Gaussian additive
noise ($\epsilon \sim \mathcal{N}(0,1)$): 

$$\tilde{Z}(\alpha)= Z(\alpha) + \tau \, \varepsilon_{\alpha}$$ 

In this case, the kriging system can be written using matrix notation as 

$$\left[
\begin{array}{cc}
\mathbf{C}+ \tau^2 \mathbf{I} & \mathbf{F}  \\
\mathbf{F}^T  & \mathbf{0}
\end{array}
\right] \times
\left[
\begin{array}{c}
\Lambda \\
\mathbf{m}
\end{array}
\right] = \left[
\begin{array}{c}
\mathbf{c} \\
\mathbf{d}
\end{array}
\right]$$ 

where the matrices and vectors used in the kriging system are:

-   the data covariance matrix,
$\mathbf{C}=[C(\alpha - \beta)]_{\alpha, \beta \in A}$

-   the drift matrix at data locations,
$\mathbf{F}=[f^{l}(\alpha)]_{\alpha\in A, l \in 0:L}$ with  $(f^0 = \mathbf{1})$,

-   the vector of the covariance between data and target,
$\mathbf{c}=[C(\alpha - \mathbf{o})]_{\alpha \in A}$,

-   the vector of the drifts at the target location,
$\mathbf{d}=[f^{l}(\mathbf{o})]_{l \in 0:L}$,

-   the vector of the kriging weights,
$\Lambda = [\lambda^{\alpha}]_{\alpha \in A}$,

-   the vector of the Lagrange coefficients ,
$\mathbf{m}=[\mu_l]_{l \in 0:L}$.

The variance of the kriging error is given by: 
$$
\text{Var}(Z(\mathbf{o}) - \hat{Z}(\mathbf{o})) = C(0) -  \Lambda^T \, \mathbf{c}  - \mathbf{m}^T \, \mathbf{d}
$$

With this model, different Best Linear Unbiased Estimators (BLUPs) can be considered to estimate different target variables. The first/left part of the kriging system is unchanged; only its second/right part must be modified to take into account the actual target variable. They are respectively,

i)  the estimation of the unobserved **variable** at target $\mathbf{o}$, $\hat{Z}(\mathbf{o})$, 
$$
\begin{array}{ccccc}
\mathbf{c} & = & \left[\text{Cov}(\tilde{Z}(\alpha), Z(\mathbf{o}))\right]_{\alpha \in A} & = & [C(\alpha - \mathbf{o})]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [f^l(\mathbf{o})]_{l \in 0:L}
\end{array}.
$$

The variance of the variable to be estimated is
$\text{Var}(Z(\mathbf{o})) = C(0)$.

ii) the estimation of the **residual** at target $\mathbf{o}$, $\hat{Y}(\mathbf{o})$, 
$$
\begin{array}{ccccc}
\mathbf{c} & = & \left[\text{Cov}(\tilde{Z}(\alpha), Y(\mathbf{o}))\right]_{\alpha \in A} & = & [C(\alpha - \mathbf{o})]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [0]_{l \in 0:L}
\end{array}.
$$

The variance of the residual to be estimated is
$\text{Var}(Z(\mathbf{o}) - m(\mathbf{o})) = C(0)$.

iii) the estimation of the **drift** at target $\mathbf{o}$, $\hat{m}(\mathbf{o}) = \hat{(\sum_{l\in 0:L} a_l\, f^l(0))}$, 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), \sum_{l\in 0:L} a_l\, f^l(\mathbf{o})\right]_{\alpha \in A} & = & [0]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [f^l(\mathbf{o})]_{l \in 0:L}
\end{array}
$$

The variance of the drift to be estimated is
$\text{Var}(m(\mathbf{o})) = 0$.

iv) the estimation of the unknown **coefficients** at target
$\mathbf{o}$ for ${l_0} \in 0:L$, $\hat{a}_{l_0}]$, 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), a_{l_0}\right]_{\alpha \in A} & = & [0]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [\delta^l_{l_0}]_{l \in 0:L}
\end{array}
$$

The variance of the coefficient to be estimated is
$\text{Var}(a_{l_0}) = 0$.

In all cases, the variance of the estimation error of $X$ is given by:
$$
\text{Var}(X - \hat{X}) = Var(X)  -  \Lambda^T \, \mathbf{c}  - \mathbf{m}^T \, \mathbf{d}
$$

### The linear model with spatially varying effects

A more complex model can be considered: some coefficients of the residual are varying (or not) across the domain of interest, the others are fixed. Latter on, the spatially varying coefficients are called the *spatial effects*, and the constant coefficients the *fixed effects*.
The spatially effects are interpreted as the realization of some stationary random functions. The model for the variable of interest is therefore: 
$$
{Z}(\alpha) = \sum_{l = 0}^{P} A_l(\alpha) \, f^l(\alpha) + \sum_{l = P+1}^{L} a_l \, f^l(\alpha).
$$ 
The error $\varepsilon \sim \mathcal{N}(0,\tau)$ may be added to the observations. This model includes the previous case:

-   we obtain the estimation of $Z$ by kriging with external drift (**KED**) if $P = 0$ and $\tau = 0$ with $L$ fixed effects,

-   we obtain the estimation of $Z$ by ordinary kriging (**OK**) if $L = 0$ and $\tau = 0$, considering only the unknown mean.

More complex cases are obtained:

-   The estimation with spatially varying effects only if $P = L$.

-   The mixed case with a residual, $P$ spatially varying effects, and $L-P$ fixed effect ($0 < P < L$).

If the varying coefficient $(A_l)_{l \in 0:P}$ are stationary random
functions with unknown means $(a_l)_{l \in 0:P}$ and cross covariance functions
$$\text{Cov}(A_l(s), A_p(s + \mathbf{h})) = C_{lp} (\mathbf{h})$$
the model of the variable of interest is second order non stationary random function whose first two moments are:

-   $m(s) = \mathbb{E}\{Z(s)\} = \sum_{l = 0}^{L} a_l \, f^l(s)$, and

-   $C(s,s') = \text{Cov}(Z(s), Z(s')) = \sum_{l,l'\in 0:P} f^{l}(s) C_{l,l'}(s - s') f^{l'}(s')$.

Using this model, different BLUPs can be proposed,
$$\hat{X} = \sum_{\alpha \in A} \lambda^{\alpha} \tilde{Z}(\alpha)$$ 
The weights are defined by the kriging system but the matrix of the covariance between observations $\mathbf{C}$ and the vector of the covariance between observations and the target $\mathbf{c}$ are computed using non stationary covariance functions.

With this model, different Best Linear Unbiased Estimators (BLUPs) can be considered. The first part of the kriging system is unchanged and only its second part must be modified to take into account the actual target variable to be estimated. They are respectively:

i)  the estimation of the unobserved **variable** at target $\mathbf{o}$, $\hat{Z}(\mathbf{o})$ 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), Z(\mathbf{o}))\right]_{\alpha \in A} & = & [\sum_{l,l' \in 0:P} f^{l}(\alpha) C_{l,l'}(\alpha - \mathbf{o}) f^{l'}(0)]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [f^l(\mathbf{o})]_{l \in 0:L}
\end{array}
$$

The variance of the variable to be estimated is
$\text{Var}(Z(\mathbf{o})) = C(\mathbf{o},\mathbf{o})$.

ii) the estimation of the **residual** at target $\mathbf{o}$, $\hat{(Z(\mathbf{o}) - m(\mathbf{o}))}$
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), Z(\mathbf{o}) - m(\mathbf{o}))\right]_{\alpha \in A} & = & [\sum_{l,l' \in 0:P} f^{l}(\alpha) C_{l,l'}(\alpha - \mathbf{o}) f^{l'}(0)]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [0]_{l \in 0:L}
\end{array}
$$

The variance of the variable to be estimated is
$\text{Var}(Z(\mathbf{o}) - m(\mathbf{o})) = C(\mathbf{o},\mathbf{o})$.

iii) the estimation of the **drift** at target $\mathbf{o}$ $\hat{m}(\mathbf{o}) = \hat{(\sum_{l\in 0:L} a_l\, f^l(0))}$ 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), \sum_{l\in 0:L} a_l\, f^l(\mathbf{o})\right]_{\alpha \in A} & = & [0]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [f^l(\mathbf{o})]_{l \in 0:L}
\end{array}
$$

The variance of the variable to be estimated is $\text{Var}(m(\mathbf{o})) = 0$.

iv) the estimation of the spatial **effects** at target $\mathbf{o}$ for
${l_0} \in 0:P$, $\hat{A}_{l_0}(\mathbf{o})]$ 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), A_{l_0}(\mathbf{o})\right]_{\alpha \in A} & = & [\sum_{l \in 0:P} f^{l}(\alpha) C_{l,l'}(\alpha - \mathbf{o}) ]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [\delta^l_{l_0}]_{l \in 0:L}
\end{array}
$$

The variance of the variable to be estimated is
$\text{Var}(A_{l_0}(\mathbf{o})) = C_{l_0, l_0}(0)$.

v)  the estimation of the fixed* *effects** at target $\mathbf{o}$ for
${l_0} \in (P+1):L$, $\hat{a}_{l_0}]$, 
$$
\begin{array}{ccccc}
\mathbf{c} & = &\left[\text{Cov}(\tilde{Z}(\alpha), a_{l_0}\right]_{\alpha \in A} & = & [0]_{\alpha \in A}\\
\mathbf{d} & = &\left[\sum_{\beta \in A} \lambda^{\beta} f^l(\beta) \right]_{l \in 0:L} & = & [\delta^l_{l_0}]_{l \in 0:L}
\end{array}
$$

The variance of the variable to be estimated is 
$\text{Var}(a_{l_0}) = 0$.

As in previous section, the variance of the estimation error of the variable $X$ is given by: 
$$
\text{Var}(X - \hat{X}) = Var(X)  -  \Lambda^T \, \mathbf{c}  - \mathbf{m}^T \, \mathbf{d}
$$

### The equations of the kriging systems

In this section, we detail the estimation of the spatially varying effect $A_{l_0}(\mathbf{o})$, for target $\mathbf{o}$ and for $l_0 \in 0:P$, using the Best Linear Unbiased Predicator (BLUP):

-   The linear predictor (LP) gives: 
$$
\hat{A}_{l_0}(\mathbf{o}) = \sum_{\alpha \in A} \lambda_{l_0}^{\alpha} \, \tilde{Z}(\alpha)
$$

-   The unbiased condition (U)
$\mathbb{E}\{A_{l_0}(\mathbf{o}) - \hat{A}_{l_0}(\mathbf{o})\} = 0$
defines the $L+1$ conditions ($l \in 0:L$): 
$$
\sum_{\alpha \in A} \lambda_{l_0}^{\alpha} \, f^l(\alpha) = \delta_{l_0}^l.
$$

-   The optimal condition (B) implies that the variance of the estimation error be minimal under the unbiased conditions 
$$
\text{Var} \{A_{l_0}(\mathbf{o}) - \hat{A}_{l_0}(\mathbf{o})\} = C_{{l_0}{l_0}}(\mathbf{o}) - 
2 \, \sum_{\alpha \in A} \lambda_{l_0}^{\alpha} c_{l_0}(\mathbf{o}, \alpha) +
\sum_{\alpha, \beta \in A} \lambda_{l_0}^{\alpha} [C(\alpha, \beta) + \tau^2 \delta(\alpha - \beta)] \lambda_{l_0}^{\beta},
$$ 
with the non stationary covariance functions 
$$
C(s, s') = \text{Cov} (Z(s), Z(s')) = \sum_{l,l' \in 0:P} f^{l}(s) \, C_{l,l'}(s - s') \,  f^{l'}(s')
$$
and 
$$
c_l(s, s') = \text{Cov} (A_l(s), Z(s')) = \sum_{l' \in 0:P} C_{l,l'}(s - s') f^{l'}(s').
$$

The kriging system to solve is therefore: 
$$
\begin{array}{ccc}
\sum_{\beta \in A} [C(\alpha, \beta) + \tau^2\delta(\alpha-\beta)] \, \lambda_{l_0}^{\beta} +
\sum_{l \in 1:L} \mu_l \, f^l(\alpha) & = & c_{l_0}(\mathbf{o},\alpha) \\
\sum_{\beta \in A} f^l(\beta) \lambda_{l_0}^{\beta} & = & \delta_{l_0}^l
\end{array}
$$

In matrix form, the kriging system boils down to

$$\left[
\begin{array}{cc}
\mathbf{C} + \tau^2 \mathbf{I} & \mathbf{F} \\
\mathbf{F}'                    & \mathbf{0}
\end{array}
\right] 
\times
\left[
\begin{array}{c}
\mathbf{\Lambda}  \\
\mathbf{m}
\end{array}
\right]=\left[
\begin{array}{c}
\mathbf{c} \\
\mathbf{d}
\end{array}
\right]$$ 
where

-   the $[n, n]$ matrices defined by $\mathbf{M} = [C(\alpha, \beta)]$ and $\mathbf{I} = [\delta(\alpha - \beta)]$,

-   the $[n, L+1]$ matrix defined by $\mathbf{F} = [f^l(\alpha)]$,

-   the $[n]$ vectors $\mathbf{\Lambda} = [\lambda_{l_0}^{\alpha}]$ and $\mathbf{c} = [c_{l_0}(\mathbf{o},\alpha)]$,

-   the $[L+1]$ vectors $\mathbf{m} = [\mu_l]$ and $\mathbf{d} = [\delta_{l_0}^l]$.

The number of observations is noted $n = \#\{A\}$.

The number of linear predictors is noted $L$.

The variance of the estimation error is: 
$$
\mathbf{Var} (A_{l_0}(\mathbf{o}) - \hat{A}_{l_0}(\mathbf{o})) = 
C_{l_0, l_0}(\mathbf{o},\mathbf{o}) - \sum_{\alpha \in A} \lambda_{l_0}^{\alpha} c_{l_0}(\mathbf{o},\alpha) -
\mu_{l_0}  =
C_{l_0, l_0}(\mathbf{o},\mathbf{o}) - \Lambda^T \mathbf{c} - \mathbf{m}^T\mathbf{d}
$$

The BLU predictor of the observed field $Z$ requires to solve a similar kriging system where only the $[L+1]$ vector $\mathbf{d} = [\delta_{l_0}^l \, f^{l}(\mathbf{o})]$ of the right part of the kriging system is modified.

The variance of the estimation error is: 
$$
\mathbf{Var} (Z(\mathbf{o}) - \hat{Z}(\mathbf{o})) = 
C(\mathbf{o},\mathbf{o}) - \sum_{\alpha \in A} \lambda_{l_0}^{\alpha} c_{l_0}(\mathbf{o},\alpha) -\sum_{l \in 0:L} \mu_l f^l(\mathbf{o}) =
C(\mathbf{o},\mathbf{o}) - \Lambda^T \mathbf{c} - \mathbf{m}^T\mathbf{d}
$$

Finally, the estimation of a fixed effect $a_{l_0}$ ($l_0 \in 0:L$) requires to modify the right part of the kriging system with the $[n]$ vector $\mathbf{c} = \mathbf{0}$ and the $[L+1]$ vector $\mathbf{d} = [\delta_{l_0}^l]$.

The variance of this estimation error is therefore: 
$$
\mathbf{Var} (\hat{a}_{l_0}) = - \mu_{l_0}
$$
