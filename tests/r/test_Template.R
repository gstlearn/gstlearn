suppressWarnings(suppressMessages(library(gstlearn)))

# - created using VH functions
#tab = VectorHelper_sequenceVD(0, 9., 1.)
#tab[3] = getTEST()
#print("Type of VD created using VH function")
#print(class(tab))

# - created using the VD constructor only
tab = VectorDouble(c(2, 3, 4, 6, 7))
print("Type of VD created using VD constructor")
print(class(tab))
print("Identifying the vector type (double)")
err = tab$identify()

# - created using VH function
#itab = VectorHelper_sequence(10)
#itab[3]        = getITEST()
#print("Type of VD created using VH function")
#print(class(itab))

# - created using the VD constructor only
itab = VectorInt(c(1, 2, 4, 6))
print("Type of VD created using VD constructor")
print(class(itab))
print("Identifying a vector of Id")
err = itab$identify()

# - manipulation of VD
err = tab$dump("Initial Vector of Double")
tab$add(tab)
err = tab$dump("Add this to itself")
tab$addCst(12)
err = tab$dump("Add a constant to this")
tabaux = tab$addVec(tab)
tabaux$dump("Add 'this' and return a new VD")

# - manipulation of VI
err = itab$dump("Initial Vector of Id")
itab$add(itab)
err = itab$dump("Add this to itself")
itab$addCst(12)
err = itab$dump("Add a constant to this")
itabaux = itab$addVec(itab)
itabaux$dump("Add 'this' and return a new VI")
