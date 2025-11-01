**Spatio-Temporal Model**

This Spatio_temporal model corresponds to the following family of SPDEs :

$$\left(\frac{\partial}{\partial t} + c\left[(1-\nabla H. \nabla)^{\alpha} + v.\nabla\right]\right)Z(s,t)=\sqrt{c}W_T(t)\otimes X_S(s)$$

where 

- $c>0$ is a time scale parameter
- $H$ is an anisotropic matrix 
- $v$ is a velocity vector
- $W_T$ is a temporal white-noise
- $X_S$ is a spatially structured  noise, solution of

    $$(1-\nabla H_S \nabla)X_S = \mathcal{W}_S$$
    