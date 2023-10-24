**Variogram Map**

The experimental variogram map is a map centered at the origin, which represents the value of experimental directional variogram across all directions $0^{\circ} \le \theta< 360^{\circ}$.

To compute an experimental variogram map, we use the function `db_vmap_compute` which we supply with the `Db` containing the data. The output is a `Db` containing a grid representing the variogram map values.