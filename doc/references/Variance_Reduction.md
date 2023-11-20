# Variance reduction with support

Here we compare the empirical variance reduction with the one given by the model :

The empirical punctual variance is obtained by 

$$\hat\sigma^2 = \frac{1}{N}\sum_{i=1}^N z^2(x_i)- \left(\frac{1}{N}\sum_{i=1}^N z(x_i))\right)^2$$

The empirical block variance 

$$\hat\sigma_v^2 = \frac{1}{N_v}\sum_{i=1}^{N_v} z^2(v_i)- \left(\frac{1}{N_v}\sum_{i=1}^{N_v} z(v_i))\right)^2$$

The true (empirical) variance reduction is $$Empirical = \hat\sigma^2-\hat\sigma_v^2$$

The variance reduction computed by the model is 

$$Theoretical = \bar{\gamma}(v,v) = C(0)-\bar{C}(v,v)$$
