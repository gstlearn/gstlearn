**Van der Corput Directions**

The van der Corput sequence is the simplest one-dimensional low-discrepancy sequence over the unit interval; it was first described in 1935 by the Dutch mathematician J. G. van der Corput. It is constructed by reversing the base-$p$ representation of the sequence of natural numbers (1, 2, 3, …).

The $p$-ary representation of the positive integer $n (\geq 1)$ is

$$
n = \sum_{k = 0}^{L-1} d_k(n) \, p^{k}
$$

where $p$ is the base in which the number $n$ is represented, and $0 \leq d_p(n) < p$, i.e. the k-th digit in the $p$-ary expansion of $n$. The $n$-th number in the van der Corput sequence is

$$
g_p(n) = \sum_{k = 0}^{L-1} d_k(n) , \frac{1}{p^{k+1}}
$$
