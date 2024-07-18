import gstlearn as gl

n = 5
I = gl.IdentityEigenCG(n)
x = gl.VectorDouble(n, 0)
b = [2, 4, 6, 8, 10]

print(b)

I.evalInverse(b, x)
print(x)