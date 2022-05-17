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

# La ligne suivante ne fonctionne pas encore car
# a["var*"] est considere comme une liste (de listes)

#a["var*"] = a["var*"]>0 
#print(a["var*"])

# La ligne suivante ne fonctionne pas car
# le mode d'entree n'est pas reconnu

#a["newvar"] = list(runif(nech), rnorm(nech))
#a$display()
#print(a["newvar*"])

# Les lignes suivantes n'ont pas de sens puisque newvar n'a pas ete cree

#v = a["newvar*"]
#v[0,0] = NA
#a["newvar*"] = v
#print(a["newvar*"])

