import gstlearn as gl

n = 5
I = gl.IdentityEigenCG(n)
x = gl.VectorDouble(n, 0)
X = gl.VectorEigen(x)

b = [2, 4, 6, 8, 10]
B = gl.VectorEigen(b)

print(b)

I.evalInverse(b, x)
print(x)

X.fill(0)
I.evalInverse(B, X)
print(X.getValues())

X.fill(0)
s = gl.LinearOpEigenCGSolver(I)
s.solve(B, X)
print(X.getValues())
