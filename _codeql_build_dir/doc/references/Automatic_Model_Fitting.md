**Automatic Model Fitting**

The aim is to find a model $\gamma$ which minimizes a weighted sum of the square difference between the model and the empirical variogram over all the computation lags:

$$\mathbf\gamma = \arg\min_\Gamma \sum_{i=1}^N w_i (\Gamma(h_i)-\hat\gamma(h_i))^2$$

Here $N$ designates the number of lags and $\hat\gamma(h_i)$ is the empirical variogram computed for lag $h_i$.

The weights $w_i$ generally give more importance to the points of the empirical variogram which have been computed with a large number of pairs $n_i$. They are also more important for the small distances. For instance, we can take  $$w_i =\frac{n_i}{h_i}$$

This problem is difficult. So we will choose a parametric form for $\mathbf\gamma$ and we will try to find the "best" parameters. More details are given [here](https://link.springer.com/article/10.1007/s11004-012-9434-1).

Note that this tool can help you to fit a model but the quality of the result is not always good. In this notebook, we also show some ways to help the algorithm to find a good result.