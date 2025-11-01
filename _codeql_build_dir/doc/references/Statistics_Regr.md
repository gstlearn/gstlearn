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

