**SPDE: Non-conditional Simulations**

The latent vector (of length $n$) $\mathbf{Z \sim \mathcal{N}(0, Q^{-1})}$ is defined on the meshing by its sparse
precision matrix $\mathbf{Q = \Sigma^{-1}}$. $\mathbf{Q}$ is factorized by the Cholesky method: $\mathbf{Q = L\, L^{T}}$.

Thus, the Gaussian vector can be rewritten $\mathbf{Z = (L^T)^{-1} \, U}$ with $\mathbf{U \sim \mathcal{N}(0, I_n)}$.

Finally the Gaussian vector collecting the values of the random field at the grid nodes $Y$ 
is achieved by the interpolation of $\mathbf{Z}$ on the mesh $\mathbf{Y_g = A_{g} \, Z}$.

To compute a non conditional simulation on the grid:

0) Compute the projection matrix of the latent vector $\mathbf{Z}$ to the grid/
1) Compute $\mathbf{L}$ the Cholesky decomposition of the sparse precision matrix $\mathbf{Q}$
2) Compute the normal Gaussian vector $\mathbf{U \sim N(0, I_n)}$,
3) Compute the latent vector $\mathbf{Z}$ solving the sparse linear system $\mathbf{L^{T} \, Z = U}$,
4) Compute the interpolation of the latent vector to the target grid $\mathbf{Y_g = A_{g} \, Z}$
