**Non-Linear Geostatistics**

We consider three supports:

* The support of the samples, considered as **points** and noted $x$,

* The support of the selection, the **bloc**, noted $v$,

* The support of the reporting, the **panels**, noted $V$.

The variable of interest $z$ is defined over the domain $S \subset \mathbb{R}^d$ with $d = 2$. The domain is uniformly divided into panels and each panel is regularly subdivided into blocs, hence defining two regular grids, the grid of panels and the grid of blocs. 

> **Hence, we will have a data set and two grids.**

The regionalised variable $z(x)$ for $x \in D \subset \mathbb{R}^d$ is modeled as a transform of the stationary Gaussian Random Function, $Z(x) = \phi(Y(x))$ where

* $Y(x)$ is a stationary Gaussian function, centered and normalized.

> $E\{Y(x)\} = 0$ and $Cov\{Y(x), Y(x+h)\} = \rho_Y(h) = 1 -\gamma_Y (h)$

* $\phi$ is a continous and one to one and mapping, called the Gaussian anamorphosis.

> As $Var(Z) < +\infty$, $\phi \in L^2(g)$ and 
> it can be expressed as a linear combination of Hermite polynomials 
> $\phi(y) = \sum_{n = 1}^{+\infty} \phi_n H_n(y)$.

Non linear Geostatistics implements non linear estimators of non linear transforms of the variable of interest $Z$. Two issues are addressed, the selection and the change of support, with the following tasks common in geosciences (i.e., in a mining context or for environmental studies), 

* predicting recovered quantities with selection over the actual value. 

> For example, recovered mineral resources will be defined by the selected ore at a given cutoff, $T(z_c) = 1_{Z \geq z_c}$, and the metal contained in the selected ore, $Q(z_c) =  Z \times 1_{Z \geq z_c}$.

* taking into account the support of the selection. 

> The average value of the volume $v$ is noted $Z(v) = \frac{1}{|v|} \int_v Z(u) du$. 

A first task is to predict marginal distribution of $Z(v)$, or the histogram, knowing the histogram of $Z$ and its spatial structure. A second task is to predict the conditional distribution of $Z(v)$ given the spatial the prior spatial model and some observations. The first question is referred to as the global recoverable resources, and the second one as the local recoverable resources. 

Three estimators will be illustrated:

* the conditional expectation (**EC**)

* the disjunctive kriging (**DK**)

> EC and DK can be used to evaluate **point** recovery functions at a non observed point $x_0$, $1_{Z(x_0) \geq z_c}$ and $Z(x_0) \times 1_{Z(x_0) \geq z_c}$. 
> They can also be used to evaluate **bloc** recovery functions for a block $v$, $1_{Z(v) \geq z_c}$ and $Z(v) \times 1_{Z(v) \geq z_c}$.
> DK can also evaluate recovery functions average on a bigger support, e.g. $\frac{1}{N} \sum_{i=1}^{N} 1_{Z(v_i) \geq z_c}$ is the recovered ore on the panels $V = \cup_{i=1}^{N}v_i$.

* the uniform conditioning (**UC**)

> UC computes block recovery functions averaged on a panel.

Two change of support models are available for a Gaussian model, they are noted respectively $DGM-1$ and $DGM-2$.
