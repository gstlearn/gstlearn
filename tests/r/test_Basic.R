suppressWarnings(suppressMessages(library(gstlearn)))

# We create a grid of 150 by 100 square cells of 1m edge.

mygrid = DbGrid_create(nx=c(100,150), dx=c(1.,1.))

#We create a geostatistical Model constituted of a single Spherical anisotropic structure with a sill of 1 and a shortest range of 30m and a longest one of 50m. The orientation of the long range is in the direction 30 degrees from East counted counter-clockwise.

mymodel = Model_createFromParam(ECov_SPHERICAL(), sill=1, ranges=c(50,30), angles=c(30,0))

#We perform one non-conditional simulation using the Turning Band method (the number of bands is set to 1000). This results in adding a field (commonly referred to as a variable) in mygrid which is called Data. Typing the name of the grid data base is an easy way to get a summary of its contents.

err = simtub(NULL, mygrid, mymodel, nbtuba=1000, namconv=NamingConvention("Data"))
mygrid

# We sample the grid on a set of 100 samples randomly located within the area covered by the grid. Typing the name of the newly created data base (mypoints) give the list of variables already stored in this data base.

mypoints = Db_createSamplingDb(mygrid, number=100)
mypoints

#We calculate a variogram along the two main directions of the Model, i.e. 30 and 120 degrees.

varioparam = VarioParam_createSeveral2D(angles=c(30,120), nlag=25, dlag=2)
myvario = Vario_computeFromDb(varioparam, mypoints)

#In the next step, we consider that the initial model (mymodel) has been correctly rendered through the simulated information at the 100 samples, as demonstrated in the experimental variogram. This is the reason why we will use this initial model for the further steps of estimation and conditional simulations.

#For this subsequent steps, we need to define the neighborhood which will serve in designating the set of samples, close to the target, which will be used for the estimation or simulation at the target location (called myneigh). Due to the small number of points, we decide to use the whole set of available samples, wherever the target is located: this is known as a Unique neighborhood.

myneigh = NeighUnique_create()

#We have all the ingredients to perform the estimation using Kriging. When typing the name of the (output) grid data base, we can check that 2 variables have been added:

#- the one corresponding to the estimation result (called Kriging.Simu.estim)
#- the one corresponding to the standard deviation of the estimation error (called kriging.Simu.stdev)

err = kriging(dbin=mypoints, dbout=mygrid, model=mymodel, neigh=myneigh)
mygrid

#We can now construct 2 conditional simulations, using the Turning Bands algorithm again, but conditioned by the variable informed at the samples of mypoints. The newly created variables in mygrid are called *Simu.Data.1" and "Simu.Data.2". Each simulation outcome reproduces the spatial characteristics provided by the Model and honors the information provided at sample points.

err = simtub(dbin=mypoints, dbout=mygrid, model=mymodel, neigh=myneigh, nbtuba=1000, nbsimu=2)
