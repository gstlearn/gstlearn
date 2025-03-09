suppressWarnings(suppressMessages(library(gstlearn)))

#Tensor for the first variable
ranges = c(1,3)
angles = c(30, 0)

#Rescaling for the 2 other variables
rescale = c(3., 2.)  

#Smoothness parameters
params = c(1.,2.,3.)

cov = CorMatern(ranges,angles,rescale,params)

#Print the maximum possible correlation (h = 0)
for (i in 0:2)
  {for (j in 0:2)
    cat(round(cov$getCorMax(i,j),3)," ")
   cat("\n")
}

#Cross-covariances between 2 points

p1 = SpacePoint(c(1.4,2.5))
p2 = SpacePoint(c(0.3,3.4))

for (i in 0:2)
{for (j in 0:2)
  cat(round(cov$eval(p1,p2,i,j),3)," ")
  cat("\n")
}

