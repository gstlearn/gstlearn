**Variogram Cloud**

The data is modeled as *samples of a regionalized* variable $z$, i.e. as evaluations at locations $x_1,..,x_n$ of a variable $z$ defined across a spatial domain: 
$$\lbrace z_i = z(x_i) : i = 1, ..., n\rbrace.$$

The variogram cloud is the set of pair of points defined as
$$ \big\lbrace \big( \Vert x_i - x_j\Vert,  \big\vert z(x_i)-z(x_j)\big\vert^2 \big) \quad\text{where}\quad 1\le i\le j\le n \big\rbrace $$

In **gstlearn**, variogram clouds are computed as grids.