**Hermite expansion for Truncated variable**

The observed variable is 

$$
Y_{y_c} = \phi_{y_c}(Y) = Y\times 1_{Y \geq y_c} + y_c \times 1_{Y < y_c} =
\sum_{n\geq 0} \phi_n(y_c) \times H_n(Y)
$$

where normalized Hermite polynomials are $H_n(n) = \frac{1}{\sqrt{n!}} \frac{g^{(n)}(y)}{g(y)}$.

The computation of the coefficients $\phi_n(y_c)$ gives

* $\phi_0(y_c) = g(y_c) + y_c \times G(y_c)$,

* $\phi_1(y_c) = G(y_c) - 1$,

* $\phi_n(y_c) = g(y_c)  \frac{H_{n-2}(y_c)}{\sqrt{n\times(n-1)}}$ for $n > 1$.
