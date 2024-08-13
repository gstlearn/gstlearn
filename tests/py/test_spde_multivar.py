import gstlearn as gl

A = gl.ProjMatrix()
v=gl.VectorProjMatrix()
v.push_back(A)
pmulti = gl.ProjMatrixMulti(v)

