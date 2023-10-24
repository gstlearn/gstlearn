**SPDE on the Sphere**

The aim of this tutorial is to show how to use gstlearn to simulate the solution of 

$$(\kappa^2-\Delta_{\mathcal{S}_R})^{\alpha/2}Z = \sigma\mathcal{W}$$

on the sphere $\mathcal{S}_R$ of radius $R$.

- $\Delta_{{\mathcal{S}_R}}$ is the Laplace-Beltrami operator, i.e, it acts on each point of the sphere as the usual Laplacian on the tangent plane at this point. 

- $\kappa$ is the inverse of the scale parameter

- $\alpha \geq 2$ is an integer describing the smoothness of the solution.

- $\mathcal{W}$ is a Gaussian white-noise suitably normalized such as $\sigma^2$ is the variance of the solution.

In this notebook, we will define the covariance of Matérn on the sphere, as the covariance of the solution of this SPDE (other extensions of the Matérn function are possible). By analogy with the Euclidian case, its smoothness parameter will be defined by $\nu = \alpha -1$. To compute the covariance function with respect on the geodetic distance, one have to use a decomposition on the Legendre polynomial (see below).

We also treat the more general case
$$P^{1/2}(-\Delta_{\mathcal{S}_R})Z = \sigma\mathcal{W}$$

where $P$ is a polynom positive on $\mathbb{R}^+$