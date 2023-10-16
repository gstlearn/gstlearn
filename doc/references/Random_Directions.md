**n-sphere and (quasi) random directions**

The $n$-sphere is defined as
$$
S^n = \{s \in \mathbb{R}^{n+1} : ||x|| = 1\}
$$
and a direction in $\mathbb{R}^{n+1}$ is a point on the half $n$-sphere. 

A simple approach to generating a uniform point on $S^n$ uses the fact that the multivariate normal distribution with independent standardized components is radially symmetric, i.e., it is invariant under orthogonal rotations. Therefore, if 
$Y \sim \mathcal{N}({\bf 0}_{n+1}, {\bf I}_{n+1})$, then $S_n = Y /||Y||$ has the uniform distribution on the unit $n$-sphere. 

For the simulation of a direction, i.e. a point on $S^n_{+}$, the last axis is selected as a reference and points with a negative coordinate along this axis are replaced by their symmetric points relatively to the origin.

Finally, in order to build a quasi random values on the half $n$-sphere, the normal variables are simulated using the inverse method and the pseudo random generator on $[0,1]^{n+1}$ is replaced by the Van der Corput sequence which is a quasi random generator on $[0,1]^{n+1}$.
