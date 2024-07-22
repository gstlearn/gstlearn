import gstlearn as gl

n = 5

b = [2, 4, 6, 8, 10]
B = gl.VectorEigen(b)

o = gl.VectorDouble(n, 0)
O = gl.VectorEigen(n)

I = gl.ScaleOp(n, 2.)
s = gl.LinearOpCGSolver(I)

print(b)

s.solve(b, o)
print(o)

s.solve(B, O)
print(O.getValues())
