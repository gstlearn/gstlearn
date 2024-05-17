### Covariance matrix

When we have several variables $z^{(1)},\dots,z^{(p)}$, we can compute their covariance matrix $\Sigma$ which stores the covariances between each pair of variable.

$$\Sigma = \left[
\begin{array}{cccc}
\textrm{var}(z^{(1)})         & \textrm{cov}(z^{(1)},z^{(2)}) &\dots  & \textrm{cov}(z^{(1)},z^{(p)})\\
\textrm{cov}(z^{(2)},z^{(1)}) & \textrm{var}(z^{(2)})         & \dots & \textrm{cov}(z^{(2)},z^{(p)})\\
\vdots & \vdots & \ddots & \vdots \\
\textrm{cov}(z^{(p)},z^{(1)}) &  \textrm{cov}(z^{(p)},z^{(2)})&\dots  & \textrm{var}(z^{(p)})\\
\end{array}\right]$$

Note that this matrix is symmetric.

If the variables (centered by their means) are stored in a matrix $Z_c$ (one column per variable), then 

$$\Sigma = \frac{1}{n} Z_c^TZ_c$$ where $^T$ designates the transposition.

In other words, $Z_c^T$ is the matrix where each line is a variable.

