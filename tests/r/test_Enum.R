suppressWarnings(suppressMessages(library(gstlearn)))

a = ECov_CUBIC()

key = a$getKey()
key

descr = a$getDescr()
descr
toLower(descr)
toUpper(descr)

val = a$getValue()
val

b = ECov_fromKey(key)
b$getDescr()
