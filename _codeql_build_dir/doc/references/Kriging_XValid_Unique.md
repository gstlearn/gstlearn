# Cross-Validation Option (Unique Neighborhood)

This is an interesting case for:

- cross-validating the value of the variable(s) at a Data Location (called the Target)
- in a multivariate case (say with N variables)
- based on the information in the input Db (note that all 'N' variable(s) do not have to be known in the 'heterotopic' case)

When working in **Unique** Neighborhood, all the  information (from the input data base) do not change for cross-validated target. But the Cross-validated information changes at each target.

Hence the interest of benefiting of the inversion of the Kriging System (Covariance and Drift parts). 

In this note, we tackle the problem when a Drift part is present (the Simple Kriging case is much simpler). We also consider the multivariate case, although this is not explicitely mentioned in the equations (we even talk of the *Kriging* system, although it should become the *CoKriging* one in multivariate case).

Let us consider the generic Kriging System

$$
    \begin{bmatrix}
        \Sigma & X \\
        X^t  & 0
    \end{bmatrix}^{-1} 
    \begin{bmatrix}
        \lambda \\
        -\mu
    \end{bmatrix}
    =
    \begin{bmatrix}
        \Sigma_0 \\
        X_0^t
    \end{bmatrix}
$$

where:

- $\Sigma$ designates the data-to-data covariance matrix
- $X$ designates the drift matrix at data location
- $\lambda$ designates the matrix of Kriging weights
- $\mu$ designates the matrix of Lagrange weights
- $\Sigma_0$ is the data-to-target covariance matrix
- $X_0$ is the drift matrix at target location

At this stage, we consider the cross-validation option:
- the data (left-hand side of the Kriging Syetem) are the actual Data (target excluded)
- the target (right-hand side) concerns the sample to be cross-validated.

The interesting feature of processing the cross-validation option in the scope of the Unique Neighborhood, comes from the fact that a global Kriging Matrix (denoted $S$) can be established (the target is *added* to the initial left-hand side and located in first position for better legibility)

$$
    S=\begin{bmatrix}
   \sigma_{00} & \Sigma_0 \\
   \Sigma_0^t  & \Sigma
    \end{bmatrix}
$$

where $\sigma_{00}$ denotes the target-to-target covariance matrix.

This matrix is never modified in Unique Neighborhood (even if the row and column corresponding to the Target sample varies). Therefore the matrix $S$ can be inverted once for all. We identify the same subdivision in the inverse matrix as in the $S$ one, introducing $\alpha$, $\beta$ and $\delta$ matrices. The inverse being calculated, these sub-matrices are considered as known in the rest of this note.

$$
  S^{-1}=\begin{bmatrix}
   \sigma_{00} & \Sigma_0 \\
   \Sigma_0^t  & \Sigma
  \end{bmatrix}^{-1}=
  \begin{bmatrix}
   \alpha & \beta \\
   \beta^t  & \delta
  \end{bmatrix}
$$

In Ordinary Kriging, the full Ordinary Kriging system can be written as follows

$$
    \begin{bmatrix}
        \sigma_{00} & \Sigma_0 & X_0\\
         \Sigma_0^t  & \Sigma & X \\
         X_0^t & X^t & 0
    \end{bmatrix}^{-1}
    \begin{bmatrix}
        0 \\
        \lambda \\
        -\mu
    \end{bmatrix}
    =
    \begin{bmatrix}
        \omega \\
         \Sigma_0^t \\
         X_0^t
    \end{bmatrix}
$$

Note the specific writing of the solution matrix: the presence of the $0$ part (at the top of this matrix) ensures that solving the systems leads to the cross-validation system described earlier in this note.

Conversely, this writing impose to determine the new unknown matrix $\omega$ present in the new right-hand side. As soon as this new parameter is calculated, the inversion of this complete Kriging System leads to the expected cross-validation results.

**How to guess the value of $\omega$?**

Expanding the first line, we get
$$
    \omega = \Sigma_0  \lambda - X_0 \mu
$$  

From the second line

$$
    \lambda = \Sigma^{-1}( \Sigma_0^t + X  \mu)
$$

Replacing in the third line

$$
    X^t \lambda = X_0^t = X^t \Sigma^{-1} ( \Sigma_0^t + X  \mu)
$$

which leads to

$$
    \mu = (X^t \Sigma^{-1} X)^{-1}  (X_0^t - X^t \Sigma^{-1} \Sigma_0)
$$

And finally the value of $\omega$

$$
    \omega = \Sigma_0 \Sigma^{-1} \Sigma_0^t + (\Sigma_0 \Sigma^{-1} X-X_0) (X^t \Sigma^{-1} X)^{-1} (X_0^t-X^t \Sigma^{-1} \Sigma_0^t)
$$

To calculate this coefficient, we need to calculate the following terms:

- $a_1 =\Sigma_0 \Sigma^{-1} \Sigma_0^t$

- $a_2 = X^t \Sigma^{-1} X$

- $a_3 = X_0^t-X^t \Sigma^{-1} \Sigma_0^t$

**Calculating the three key terms**

We need to recall the Schur decomposition of the generic matrix
$$
  \begin{bmatrix}
   A & B \\
   B^t  & D
  \end{bmatrix}
$$

We introduce the term

$$
    C = (A - BD^{-1}B^t)^{-1}
$$

Then

$$
  \begin{bmatrix}
   A & B \\
   B^t  & D
  \end{bmatrix}^{-1}=
  \begin{bmatrix}
    C & -CBD^{-1} \\
    -D^{-1}B^tC & D^{-1} + D^{-1}B^tCBD^{-1} 
  \end{bmatrix}
$$

We apply this Schur decomposition to the matrix $S$.

From the top-left term of the Schur decomposition

$$
    \alpha = (\sigma_{00} - \Sigma_0  \Sigma^{-1} \Sigma_0^t)^{-1}
$$

which directly leads to the value of $a_1$

$$
 a_1 = \sigma_{00} - \alpha^{-1}
$$

The top-right term of the Schur decomposition leads to 

$$
    \beta = -\alpha \times \Sigma_0 \times \Sigma^{-1}
    \implies
    \Sigma_0 \times \Sigma^{-1} = - \alpha^{-1} \times \beta
$$

Multiplying by $X$ to the right and introducing the term $\epsilon$

$$
    \epsilon = \Sigma_0 \times \Sigma^{-1} \times X = - \alpha^{-1} \times \beta \times X
$$

We can derive the value of $a_3$

$$
    a_3 = X_0^t - \epsilon^t
$$

From the lower-right term of the Schur decomposition

$$
    \delta = \Sigma^{-1} + \Sigma^{-1} \Sigma_0^t \alpha  \Sigma_0  \Sigma^{-1}
$$

Multiplying by $X^t$ on the left and by $X$ on the right, it comes

$$
    X^t \delta  X = 
    (X^t  \Sigma^{-1} X) + 
    (X^t \Sigma^{-1} \Sigma_0^t) \alpha (\Sigma_0 \Sigma^{-1} X)
$$

We finally derive the value of $a_2$

$$
    a_2 = X^t \Sigma^{-1} X = (X^t \delta  X) - (\epsilon^t  \alpha \epsilon)
$$

This prooves that we can derive the cross-validation vector of weight and Lagrange multipliers from the global matrix inverted once for all.
