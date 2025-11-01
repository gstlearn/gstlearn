#### Multiple linear regression

When we have several variables $x^{(1)},\dots,x^{(p)}$ to explain an interest variable $y$ we can also use a linear regression

$$y=\sum_{j=1}^p \beta_j x^{(j)} + \beta_0 + R$$

Note that for convenience, we will rewrite the relation 

$$y=\sum_{j=0}^p \beta_j x^{(j)}+R$$ 

where the variable $x^{(0)}$ is equal to $1$.

Last, we can rewrite more compactly

$$y = \beta^T X +R$$

where $$\beta = \left[\begin{array}{c}\beta_0 \\ \vdots \\ \beta_p\end{array}\right]$$

and $X$ is the table with all the observations. The first column contains $1$'s and then each column is a variable 
$$X  = \left[\begin{array}{cccc} 1 & x^{(1)} & \dots & x^{(p)}\end{array}\right]$$

As in the simple linear regression case, we will try to minimize

$$||R||^2=||y-\beta^TX||^2$$

We can show that 

$$\hat\beta = (X^TX)^{-1}X^Ty$$

