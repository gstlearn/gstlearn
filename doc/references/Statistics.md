# Statistics #

### Location

#### Mean

$$bar{z} = \frac{1}{n}\sum_{i=1}^n z_i$$

#### Median

The median $m$ of the set of values is a value such as half of the total observations are below and half is above.

* If the size $n$ is odd, it is the value of the $z_{\left(\frac{n+1}{2}\right)}$ where $z_{(i)}$ is the value of the $i$ th observation (when ordered in increasing order).

* If $n$ is even, we can take $\frac{z_{\left(\frac{n}{2}\right)}+z_{\left(\frac{n}{2}+1\right)}}{2}$.

The median is less sensitive than the mean to extreme values.

#### Quartiles

The quartiles are the values which divide the samples as follows:

* lower quartile: 25\% of the individuals are below
* upper quartile: 25\% of the individuals are above

####  Quantiles (p)

We can generalize to any proportion $p$.

The $p-$quantile denoted $q_p$ is a value which divides the samples such as a proportion $p$ of the individuals are below the quantile. 

The median is $q_{\frac{1}{2}}$. The lower quartile is $q_{\frac{1}{4}}$ and the upper quartile is $q_{\frac{3}{4}}$.

### Dispersion

The dispersion summaries try to measure how the data are spread.

#### Range

$$\max_{i=1,\dots,n}{z_i}-\min_{i=1,\dots,n}{z_i}$$

#### Inter-quartiles distance

$$q_{\frac{3}{4}}-q_{\frac{1}{4}}$$

#### Variance 

The variance measures the average distance between each individual to the mean:

$$\frac{1}{n}\sum_{i=1}^n(z_i-\bar{z})^2$$

It is sometimes more convenient to use the equivalent formula


$$\frac{1}{n}\sum_{i=1}^nz_i^2-\bar{z}^2 = \bar{z^2}-\bar{z}^2$$

Note that for statistical reasons, one often prefers to use

$$\frac{1}{n-1}\sum_{i=1}^n(z_i-m)^2$$

for the variance. The two formulas give close results when $n$ is large.

To have a measure in the same unit as the variable, one often consider the standard deviation.

$$\sqrt{\frac{1}{n}\sum_{i=1}^n(z_i-\bar{z})^2}$$

## Distribution

### Histogram

To have a good idea of the distribution of a variable, one can compute the histogram.

The idea is 

* divide the range of the variable $[min,Max]$ into small intervals. Here, we only treat the case were all intervals have the same size
* compute the number of samples in each interval.


Normalized histogram rescales the ordinate such as the total surface is equal to 1.

### Cumulated histogram

We can represent the cumulated histogram. It is a function which computes, for each value, the proportion of individuals below this value. 
It can be written as 
    
   $$F(z_c) =\frac{1}{n}\sum_{i=1}^n 1\!\!\!1_{]z_{i},+\infty]}(z_c)$$
   
   where $1\!\!\!1_A$ is the indicator function of the set $A$:
   
   $$1\!\!\!1_A(x)=\left\{\begin{array}{ccc}1 &\textrm{ if } & x\in A\\
   0 & \textrm{ otherwise } & \end{array}
   \right.$$

### Quantile function

If we inverse the two axes, we obtain the quantile function which gives, for each value $p$, the quantile $q_p$.

$$q(p) = F^{-1}(p)$$

### Ore

In mine, we often consider the ore function 
$$T(z_c) = 1-F(z_c)$$

Indeed, it gives the proportion of the data which are above a cut-off.

### Metal

$$Q(z_c) =\frac{1}{n}\sum_{i=1}^n z_i1\!\!\!1_{]z_{i},+\infty]}(z_c)$$


### Grade 

$$m(z_c)=\frac{Q(z_c)}{T(z_c)}$$

#### $Q(T)$ curve

We just represent the **Metal** with respect to the **Ore** for various cut-off values $z_c$.

#### Conventional benefit

$$B(z_c) = Q(z_c)-z_cT(z_c)$$

## Bivariate statistics

Now we consider two variables:

* $z^{(1)}=(z_1^{(1)},\dots,z_n^{(1)})$
* $z^{(2)}=(z_1^{(2)},\dots,z_n^{(2)})$

and we will study their relationship.

### Covariance

We can compute the covariance between the two vectors $z^{(1)}$ and  $z^{(2)}$.

$$\textrm{cov}(z^{(1)},z^{(2)}) = \frac{1}{n}\sum_{i=1}^n (z^{(1)}_i-\bar{z}^{(1)})(z^{(2)}_i-\bar{z}^{(2)})$$

where $\bar{z}^{(j)}$ is the mean of the variable $z^{(j)}$ with $j=1,2$.

### Correlation coefficient

The covariance depends on the scale of $z^{(1)}$ and $z^{(2)}$. In order to have a scale invariant measure, we can use the correlation coefficient 
$$\rho = \frac{\textrm{cov}(z^{(1)},z^{(2)})}{\sqrt{\textrm{var}(z^{(1)})\textrm{var}(z^{(2)})}}$$

The correlation coefficient lies within $[-1,1]$.

When it is equal to $-1$ or $1$, the variables are linked by a linear relationship

$$z^{(2)}=a.z^{(1)}+b$$

where the sign of $a$ corresponds to the sign of $\rho$.

When $\rho=0$, we say that the variables are uncorrelated. But they can still have a link (not linear).

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

### Scatter plot

We can represent the scatter plot between the two variables (only isotopic samples are represented).

Here the relation could be considered as linear. Let's try to find the coefficents of the regression line.

### Linear regression

#### Simple linear regression

We can model the relationship between $z^{(1)}$ and $z^{(2)}$ by using a linear regression.
 model 
$$z^{(2)}=az^{(1)}+b + R$$ where $R$ is a residual.

We try to find $(a,b)$ by minimizing the sum of the squared difference between $z^{(2)}$ and $az^{(1)}+b$: 

$$||R||^2 =\sum_{i=1}^n(z^{(2)}_i - (az^{(1)}_i+b))^2.$$

We can show that the coefficients $a$ and $b$ can be estimated by

$$\hat a = \frac{\textrm{cov}(z^{(1)},z^{(2)})}{\textrm{var}(z^{(1)})}$$

and $b$ by 

$$\hat b = \bar{z}^{(2)}-\hat a\bar{z}^{(1)}$$

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

### Regression

To represent the two variables, we can perform a 2d histogram.


Then we could look at the histogram of $z_2$ for a given class of $z_1$.

For instance, if we consider the 3rd class, $z_1\in[1.14,1.67]$ :

It shows the conditional distribution of $z_2$ knowing that $z_1\in[1.14,1.67]$.

#### Conditional mean (or regression)

In the same spirit, we can consider the conditional mean (mean of $z_2$ for different class of $z_1$). 

It is named conditional mean (or regression).

