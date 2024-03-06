**Covariance on the Sphere**

The covariance between two points with great-circle distance $d$  on the sphere of radius $R$ is given by
$$C(d) = \frac{\sigma^2}{\sum_{i=0}^\infty f(i)}\sum_{i=0}^\infty f(i) P_i\left(\cos \frac{d}{R}\right)$$

where the $P_i$'s  are the Legendre polynomials computed with the following reccurence formula

$$P_0(x) = 1.$$

$$P_1(x) = x$$

$$P_{n+1}(x)=\frac{(2n+1)xP_n(x) - n P_{n-1}(x)}{n+1}$$

For $n\geq 0$, $$f(n) = \frac{2 n}{ (R^2\kappa^2 + n ( n + 1))^\alpha}$$

For numerical computations, the sums are truncated at **N**.

For more details on the covariances on sphere, see 
[Lantuejoul, Freulon and Renard (2019)](https://link.springer.com/content/pdf/10.1007/s11004-019-09799-4.pdf)
