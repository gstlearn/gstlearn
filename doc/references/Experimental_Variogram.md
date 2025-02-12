**Experimental Variogram**

The experimental (isotropic) variogram $\gamma$ is a  function defined as

$$\gamma(h)=\frac{1}{2\vert N(h)\vert}\sum_{(i,j) \in N(h)}\big\vert z(x_i)-z(x_j)\big\vert^2, \quad h\ge 0,$$

where $N(h)$ is set of all pairs of data points separated by a distance $h$ (called *lag*):
$$ N(h) = \bigg\lbrace (i,j) : \Vert x_j-x_i\Vert = h\bigg\rbrace_{1\le i\le j\le n},$$

and $\vert N(h)\vert$ is the cardinal of $N(h)$. In practice, when computing $\gamma(h)$, we look for pairs of data points separated by a distance $h \pm \tau h$ where $\tau > 0$ is a tolerance on the separation distance $h$. In other words, $N(h)$ is replaced by
$$ \widehat N(h) = \bigg\lbrace (i,j) : (1-\tau)h \le \Vert x_j-x_i\Vert \le (1+\tau) h\bigg\rbrace_{1\le i\le j\le n}$$

To compute an experimental variogram, we start by creating a `VarioParam` object containing the parameters of the variogram. This is done using the function `VarioParam_createOmniDirection`. We can specify the number of lags $h$ for which the experimental variogram is computed (argument `nlag`), and the distance between these lags (argument `dpas`), as well as the tolerance $\tau$ on the lags (argument `toldis`).

Then, the experimental variogram is computed in two steps. First, a `Vario` object is initialized from the `VarioParam` object  and the `Db` containing the data points. Then, the values of the experimental variogram at the lags specified by  the `VarioParam` object  are computed using the method `compute` of the `Vario` object (which returns an error code, `0` meaning that no error was detected).

Note : The variable $z$ for which we wish to define the experimental variogram should be the only variable in the `Db` with a `z` locator (i.e. it should have locator `z1` and the other variables should not have a locator starting with `z`). This can be done bu using the method `setLocator` of the `Db` object containing the data. If several variables with `z` locators are present in the `Db`, then cross-variograms between are also computed (this subject will be covered in the course on multivariate analysis). 

In the next example, we compute an experimental variogram with $40$ lags separated by a distance $10$ (meaning that we take $h =10i$ for $i=0, ..., 39$), and consider a tolerance $\tau = 10\%$ for the variogram computations. We use the `Db` `dat`, and select the variable `January_temp` as our variable of interest (by setting its locator to "z").
