suppressWarnings(suppressMessages(library(gstlearn)))

##################
# Vector of Double
##################
err = mestitle(0, "Double values")

# - created using VH functions
tab = VectorHelper_sequenceVD(0, 9., 1.)
tab[3] = getTEST()
print("Type of VD created using VH function")
print(class(tab))

# - created using the VD constructor only
tab = VectorDouble(c(2, 3, 4))
print("Type of VD created using VD constructor")
print(class(tab))
print("Identifying the vector type (double)")
err = tab$identify()

##############
# Vector of Id
##############
err = mestitle(0, "Integer values")

# - created using VH function
itab = VectorHelper_sequence(10)
itab[3]        = getITEST()
print("Type of VD created using VH function")
print(class(itab))

# - created using the VD constructor only
itab = VectorInt(c(1, 2, 4, 6))
print("Type of VD created using VD constructor")
print(class(itab))
print("Identifying a vector of Id")
err = itab$identify()
