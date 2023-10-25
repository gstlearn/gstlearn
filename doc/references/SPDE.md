**SPDE**

Lindgren et al. (2011) defined an explicit link between Gaussian fields and Gaussian markov random fields. This  approach is implemented in the SPDE module of $gstlearn$. The random field is a weak solution of some stochastic partial differential equation solved using the finite element method. It follows that the standard tools of Geostatistics, kriging and stochastic simulation, can be rewritten using sparse linear algebra.

The SPDE model represents the domain by a mesh and random field by its values for each mesh cell 
$\mathbf{Z} \sim \mathcal{N}(\mathbf{0,Q^{-1}})$ where the precision matrix $\mathbf{Q = \Sigma^{-1}}$ is sparse. 

The latent vector $\mathbf{Z}$ can be linearly interpolated to:

* the target locations $\mathbf{Y_T = A_g \, Z}$, or 

* the observation locations $\mathbf{Y = A_d \, Z + \tau W}$. 

In the latter case, $\tau$ is the standard deviation of the noise modelling the error of observation and the modelling error. The number of observations is $p$.

In this case, the precision matrix of the vector $\mathbf{(Z, Y_D)}$ is:
$$
  \mathbf{
    \tilde{Q} = \tilde{\Sigma}^{-1}=
      \begin{bmatrix}\mathbf{Q+\tau^{-2}A_d^T A_d} & \mathbf{-\tau^{-2}A_d^T} \\ \mathbf{-\tau^{-2}A_d} & \mathbf{\tau^{-2}I_p} \end{bmatrix}
  } 
$$
From this expression the kriging and the conditional simulation can be derived:
    
* the Kriging of $\mathbf{Z}$ is $E\{\mathbf{Z|Y = y}\} = \tau^{-2}\mathbf{(Q + \tau^{-2}A_d^TA_d)^{-1}A_d^T y}$

* the conditional variance is $Cov\{\mathbf{Z|Y = y}\} = \mathbf{(Q + \tau^{-2}A_d^TA_d)^{-1}}$

* the non conditional simulation is $\mathbf{(S, S_D)\sim \mathcal{N}(0, \tilde{Q}^{-1})}$
        
* the conditional simulation of the latent vector is $\mathbf{S_{|Y=y} = S + \tau^{-2} (Q+ \tau^{-2}A_d^TA_d)^{-1} A_d^T(y - S_D)}$
        
The estimated or simulated latent vector can be linearly interpolated to any target grid.
