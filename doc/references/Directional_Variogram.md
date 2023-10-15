**Directional Variogram**

The experimental directional variogram $\gamma$ is a  function defined as
$$\gamma(\theta,h)=\frac{1}{2\vert N(\theta, h)\vert}\sum_{(i,j) \in N(\theta, h)}\big\vert z(x_i)-z(x_j)\big\vert^2, \quad 0^{\circ}\le \theta <360^{\circ}, \quad h\ge 0$$

where $N(\theta, h)$ is set of all pairs of data points separated by a vector of size $h$ and along the direction $\theta$ (in degrees):
$$ N(\theta, h) = \bigg\lbrace (i,j) : \Vert x_j-x_i\Vert = h \quad\text{and the vector } \vec{u}=(x_j-x_i) \text{ is along the direction } \theta\bigg\rbrace_{1\le i\le j\le n},$$

In practice, when computing $\gamma(\theta, h)$, we once gain consider a tolerance $\tau$ on the separation distance $h$, and also consider a tolerance $\eta>0$ is also considered for the direction angle. In other words, $N(h)$ is replaced by
 $$\widehat N(\theta, h) = \bigg\lbrace (i,j) : (1-\tau)h \le \Vert x_j-x_i\Vert \le (1+\tau) h \quad\text{and the vector } \vec{u}=(x_j-x_i) \text{ is along the direction } \theta \pm \eta \bigg\rbrace_{1\le i\le j\le n},$$
 
 Much like their isotropic counterparts, experimental directional variograms are computed as `Vario` objects, which can be created from he `VarioParam` object (containing the parameters of the variogram) and a `Db` containing the data points. 

This time, the `VarioParam` object is created using the function `VarioParam_createMultiple`. There, we specify the number $K$ of directions $\theta$ for which we wish to compute the an experimental variogram (argument `ndir`), as well as the reference angle $\theta_0$ of the first direction (argument `angref`, default = $0$) so that the directions $\theta$ = $\theta_0 + i(180/K)$ for $i=0,..., K-1$ are considered. We can also specify the number of lags $h$ for which the experimental variogram is computed (argument `npas`), and the distance between these lags (argument `npas`), as well as the tolerance $\tau$ on the lags (argument `toldis`). Then, the experimental variogram is computed just as in the isotropic case.

Note: When initializing the `VarioParam` object as described above, the angle tolerance $\eta$ is automatically set to $\eta = (90/K)$, meaning that we span the set of possible directions.

In the following example, we create an experimental variogram in the $4$ directions $\theta = 0^{\circ}, 45^{\circ}, 90^{\circ}, 135^{\circ}$.