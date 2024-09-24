import gstlearn as gl
n = 5
b = [2, 4, 6, 8, 10]
o = gl.VectorDouble(n, 0.)
c = gl.VectorDouble(n, 0)
I = gl.ScaleOp(n, 2.)
s = gl.LinearOpCGSolver(I)
print(b)
s.solve(b, o)
o.display()
I.evalDirect(o, c)
c.display()