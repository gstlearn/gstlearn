suppressWarnings(suppressMessages(library(gstlearn)))



# - created using the VD constructor only
rtab1 = VectorDouble(c(2, 3, 4, 6, 7))
rtab1$dump("rtab1")
err = rtab1$identify()
print("Moyenne = ", rtab1$mean())
print("Mediane = ", rtab1$median())
print("Minimum = ", rtab1$minimum())
print("Maximum = ", rtab1$maximum())
print("Variance = ", rtab1$variance())
print("St. dev. = ", rtab1$stdv())

rtab2 = VectorDouble(c(-1.1, np.nan, 6.1, 4.1, -0.2))
rtab2$dump("rtab2")
rtab1$add(rtab2)
rtab1$dump("'rtab1 = 'rtab1 + 'rtab2'")
rtab1$addCst(12.)
rtab1$dump("'rtab1' = 'rtab1' + 12")
tabaux = rtab1$addVec(rtab2)
tabaux$dump("'out' = 'rtab1' + 'rtab2'")

# - created using the VD constructor only
itab1 = VectorInt(c(1, 2, 4, 6, 3, 5))
itab1$dump("itab1");
itab1$identify()
print("Moyenne = ", itab1$mean())
print("Mediane = ", itab1$median())
print("Variance = ", itab1$variance())
print("St. dev. = ", itab1$stdv())

itab2 = VectorInt(c(-5, 4., 2, 0, 1, 8))
itab2$dump("itab2");
itab1$add(itab2)
itab1$dump("'itab1 = 'itab1 + 'itab2'")
itab1$addCst(12)
itab1$dump("'itab1' = 'itab1' + 12")
itabaux = itab1$addVec(itab2)
itabaux$dump("'out' = 'itab1' + 'itab2'")