suppressWarnings(suppressMessages(library(gstlearn)))

# We create a grid of 150 by 100 square cells of 1m edge.

mygrid = DbGrid_create(nx=c(100,150), dx=c(1.,1.))

# We create a geostatistical Model constituted of a single Spherical anisotropic structure with a sill of 1 and a shortest range of 30m and a longest one of 50m. The orientation of the long range is in the direction 30 degrees from East counted counter-clockwise.

mymodel = Model_createFromParam(ECov_SPHERICAL(), sill = 1, ranges = c(50, 30), angles = c(30, 0))
mymodel$display()

#We perform one non-conditional simulation using the Turning Bands method (with a large number of bands set to 10000). This results in adding a field (commonly referred to as a variable) in mygrid which is called Data. Typing the name of the grid data base is an easy way to get a summary of its contents.

nbtuba = 10000
seed = 14341
err = simtub(NULL, mygrid, mymodel, nbtuba=nbtuba, namconv=NamingConvention("Data"), seed=seed)
mygrid$display()

########################
p = plot.init(asp=1)
p = p + plot.raster(mygrid, name="Data", flagLegend=TRUE, legendName="Simu", limits=c(-3.5, 3.5))
p = p + plot.decoration(title="Initial Simulated field", xlab="Easting", ylab="Northing")
ggsave("Figure1_r.png", p, dpi=600)
plot.end(p)
########################

# We sample the grid on a set of 100 samples randomly located within the area covered by the grid. Typing the name of the newly created data base (mypoints) give the list of variables already stored in this data base.

mypoints = Db_createSamplingDb(mygrid, number=100)
mypoints$display()

#We calculate a variogram along the two main directions of the Model, i.e. 30 and 120 degrees.

varioparam = VarioParam_createSeveral2D(angles=c(30,120), nlag=25, dlag=2)
myvario = Vario_computeFromDb(varioparam, mypoints)
myvario$display()

########################
p = plot.init()
p = p + plot.vario(myvario, model=mymodel, flagLegend=TRUE)
p = p + plot.decoration(title="Variogram in 2 main directions")
ggsave("Figure2_r.png", p, dpi=600)
plot.end(p)
########################

#In the next step, we consider that the initial model (mymodel) has been correctly rendered through the simulated information at the 100 samples, as demonstrated in the experimental variogram. This is the reason why we will use this initial model for the further steps of estimation and conditional simulations.

#For this subsequent steps, we need to define the neighborhood which will serve in designating the set of samples, close to the target, which will be used for the estimation or simulation at the target location (called myneigh). Due to the small number of points, we decide to use the whole set of available samples, wherever the target is located: this is known as a Unique neighborhood.

myneigh = NeighUnique_create()

#We have all the ingredients to perform the estimation using Kriging. When typing the name of the (output) grid data base, we can check that 2 variables have been added:

#- the one corresponding to the estimation result (called Kriging.Simu.estim)
#- the one corresponding to the standard deviation of the estimation error (called kriging.Simu.stdev)

err = kriging(dbin=mypoints, dbout=mygrid, model=mymodel, neigh=myneigh)
mygrid$display()

########################
p1 = plot.init(asp=1)
p1 = p1 + plot.raster(mygrid, name="Kriging.Data.estim", flagLegend=TRUE, legendName="Estim", limits=c(-3.5, 3.5))
p1 = p1 + plot.symbol(mypoints, nameSize="Data", color="yellow", sizeRange=c(0.1, 2))
p1 = p1 + plot.decoration(title="Estimation", xlab="Easting", ylab="Northing")
p2 = plot.init(asp=1)
p2 = p2 + plot.raster(mygrid, name="Kriging.Data.stdev", palette="Spectral", flagLegend=TRUE, legendName="St.Dev.")
p2 = p2 + plot.symbol(mypoints, nameSize="Data", flagCst=TRUE, color="black")
p2 = p2 + plot.decoration(title="St. Dev. of Estimation Error", xlab="Easting", ylab="")
p = ggarrange(p1, p2, ncol = 2, nrow = 1)
plot.end(p)
ggsave("Figure3_r.png", p, dpi = 600)
########################

# We can now construct 2 conditional simulations, using the Turning Bands algorithm again, but conditioned by the variable informed at the samples of mypoints. The newly created variables in mygrid are called *Simu.Data.1" and "Simu.Data.2". Each simulation outcome reproduces the spatial characteristics provided by the Model and honors the information provided at sample points.

err = simtub(dbin = mypoints, dbout = mygrid, model = mymodel, neigh = myneigh, nbtuba = 1000, nbsimu = 2)
mygrid$display()

########################
p1 = plot.init(asp=1)
p1 = p1 + plot.raster(mygrid, name="Simu.Data.1", limits=c(-3.5, 3.5))
p1 = p1 + plot.symbol(mypoints, nameSize="Data", color="yellow", sizeRange=c(0.1, 2))
p1 = p1 + plot.decoration(title="Simulation 1", xlab="Easting", ylab="Northing")
p2 = plot.init(asp=1)
p2 = p2 + plot.raster(mygrid, name="Simu.Data.2", flagLegend=TRUE, legendName="Simu2", limits=c(-3.5, 3.5))
p2 = p2 + plot.symbol(mypoints, nameSize="Data", color="yellow", sizeRange=c(0.1, 2))
p2 = p2 + plot.decoration(title="Simulation 2", xlab="Easting", ylab="")
p = ggarrange(p1, p2, ncol = 2, nrow = 1)
plot.end(p)
ggsave("Figure4_r.png", p, dpi=600)
########################
