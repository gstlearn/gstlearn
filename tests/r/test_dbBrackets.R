# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))

set.seed(32421)
a = DbGrid_create(nx=c(2,2),dx=c(1.,1.))
nech = a$getSampleNumber()
a$display()

a["var1"] = rnorm(nech)
a$display()

print(a["var1"])
a["var1"] = rnorm(nech)
print(a["var1"])

a["var2"] = rnorm(nech)
a$display()
print(a["var*"])

a["var*"] = a["var*"]>0 
print(a["var*"])

a["newvar"] = c(runif(nech), rnorm(nech))
a$display()
print(a["newvar*"])

v = a["newvar*"]
v[1,1] = NA
a["newvar*"] = v
print(a["newvar*"])

