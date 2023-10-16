**Hermite Polynomials**

The objective is to compute

$$
\psi_1(\alpha,\beta)=\int \phi(\alpha + \beta \, u) \, g(u) \, du
$$

and

$$
\psi_2(\alpha,\beta) = \int \phi^2(\alpha + \beta \, u) \, g(u) \, du
$$

for a function $\phi$ in order to compute the conditional expectation and the conditional variance, usually defined by its expansion in terms of Hermite polynomials:

$$
\phi(y) = \sum_{n=0}^{N} a_n \eta_n(y)
$$

In particular, we may be interested in the case where $\alpha=r y$ and $\beta = \sqrt { { 1-r^2 } }$. Note that $|\beta|< 1$.

Many methods can be used to compute these draw.matrix:

* the Monte-Carlo integration (easy)

* the computation of the Hermite coefficient for $\phi$ and $\phi^2$ 

$$
\int_{-\infty}^{+\infty} \eta_n(\alpha + \beta \, u)\, g(u) \, du = (1-\beta^2)^n \, \eta_n(\frac{\alpha}{\sqrt{1-\beta^2}})
$$
