# Covariance over block

To compute $$\bar{C}(v,v)=\frac{1}{|v|^2}\int_v\int_v C(x,y)dxdy$$, you need to define $v$ and the model.

The support $v$ must be defined by its extension spanned over the space dimension.

We also need to define the discretization $[N_1,N_2]$ since $\bar{C}(v,v)$ is approximated by

$$\frac{1}{N_1}\frac{1}{N_2}\sum_{i=1}^{N_1}\sum_{j=1}^{N_2} C(x_i,x_j)$$ where $x_i$ and $x_j$ are some points inside the block $v$.
