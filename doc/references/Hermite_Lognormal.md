**Hermite expansion in Lognormal case**

We use the lognormal model to test some computations done with Hermite polynomials.

We consider the second order stationary model $Z_{\lambda}(x) = e^{\lambda Y(x) - \frac{1}{2} \lambda^2}$ where $Y(x)$ is a centered Gaussian model with autocorrelation $\rho$:

* $E\{Y\} = 0$

* $Cov(Y(x), Y(x+h)) = E\{Y(x)Y(x+h)\} = \rho(h)$

* $E\{Z_{\lambda}\} = 1$

* $Cov(Z_{\lambda}(x), Z_{\lambda}(x+h)) = E\{Z_{\lambda}(x)Z_{\lambda}(x+h)\} - 1 = e^{\lambda^2\rho(h)} - 1$


The anamorphosis $\phi_{\lambda}$ maps the Gaussian field into the lognormal SOS model 
$Z_{\lambda} = \phi_{\lambda}(Y) = \sum_{n=0}^{+\infty}\phi_n{(\lambda)}H_n(Y)$.

The Hermite coefficients are $\phi_n{(\lambda)} = \frac{(-\lambda)^n}{\sqrt{n!}}$ and we have
$$
C_{\lambda}(h) = e^{\lambda^2 \rho } - 1 = \sum_{n = 1}^{+\infty} \phi_n^2(\lambda) \rho^{n}.
$$

The recursion formula to compute the values of the Hermite polynomials are: $H_0(y) = 1$, $H_1(y)= -y$, and
$$
H_{n+1}(y) = -\frac{1}{\sqrt{n+1}} y H_n(y) - \sqrt{\frac{n}{n+1}} H_{n-1}(y)
$$
