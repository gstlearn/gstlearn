suppressWarnings(suppressMessages(library(gstlearn)))

# - Based on VectorDouble
rtab1 = VectorDouble(c(2, 3, 4, 6, 7))
err = rtab1$dump("'rtab1'", FALSE)
err = rtab1$identify()
print(paste0("Moyenne = ", round(rtab1$mean(), 4)))
print(paste0("Mediane = ", round(rtab1$median(), 4)))
print(paste0("Minimum = ", round(rtab1$minimum(), 4)))
print(paste0("Maximum = ", round(rtab1$maximum(), 4)))
print(paste0("Variance = ", round(rtab1$variance(), 4)))
print(paste0("St. dev. = ", round(rtab1$stdv(), 4)))

rtab2 = VectorDouble(c(-1.1, NA, 6.1, 4.1, -0.2))
err = rtab2$dump("'rtab2'", FALSE)
tabaux = VectorDouble(VectorHelper_addVD(rtab1, rtab2))
err = tabaux$dump("'out' = 'rtab1' + 'rtab2'", FALSE)
tabaux = VectorDouble(VectorHelper_addVDCst(rtab1, 3.2))
err = tabaux$dump("'out' = 'rtab1' + 3.2", FALSE)

# - Based on VectorInt
itab1 = VectorInt(c(1, 2, 4, 6, 3, 5))
err = itab1$dump("'itab1'", FALSE);
err = itab1$identify()
print(paste0("Moyenne = ", round(itab1$mean(), 4)))
print(paste0("Mediane = ", round(itab1$median(), 4)))
print(paste0("Variance = ", round(itab1$variance(), 4)))
print(paste0("St. dev. = ", round(itab1$stdv(), 4)))

itab2 = VectorInt(c(-5, 4., 2, 0, 1, 8))
err = itab2$dump("'itab2'", FALSE)
itabaux = VectorInt(VectorHelper_addVI(itab1, itab2))
err = itabaux$dump("'out' = 'itab1' + 'itab2'", FALSE)
itabaux = VectorInt(VectorHelper_addVICst(itab1, 2))
err = itabaux$dump("'out' = 'itab1' + 2", FALSE)
