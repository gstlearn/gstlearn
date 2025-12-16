# # Script in python

# Next we load the packages used by the subsequent script, such as *gstlearn* (note that any function belonging to this library will have a prefix *gl*. Similarly, we load the graphic extension (called *gslearn.plot*).

import gstlearn as gl


# We create a grid of 150 by 100 square cells of 1m edge.

mygrid = gl.DbGrid.create(nx=[100, 150], dx=[1.0, 1.0])

# We create a geostatistical Model constituted of a single Spherical anisotropic structure with a sill of 1 and a shortest range of 30m and a longest one of 50m. The orientation of the long range is in the direction 30 degrees from East counted counter-clockwise.

mymodel = gl.Model.createFromParam(
    gl.ECov.SPHERICAL, sill=1, ranges=[50, 30], angles=[30, 0]
)
mymodel.display()

# We perform one non-conditional simulation using the Turning Bands method (with a large number of bands set to 10000). This results in adding a field (commonly referred to as a *variable*) in *mygrid* which is called *Data*. Typing the name of the grid data base is an easy way to get a summary of its contents.

nbtuba = 10000
seed = 14341
err = gl.simtub(
    None, mygrid, mymodel, nbtuba=nbtuba, namconv=gl.NamingConvention("Data"), seed=seed
)
mygrid.display()

########################
# fig, ax = plt.subplots()
# plt.gca().set_aspect('equal')
# ax.raster(mygrid, name="Data", flagLegend=True, vmin=-3.5, vmax=3.5)
# ax.decoration(title="Initial Simulated Field", xlabel="Easting", ylabel="Northing")
# fig.savefig('Figure1.png', dpi=600)
# plt.close(fig)
########################

# We sample the grid on a set of 100 samples randomly located within the area covered by the grid. Typing the name of the newly created data base (*mypoints*) give the list of variables already stored in this data base.

mypoints = gl.DbGrid.createSamplingDb(mygrid, number=100)
mypoints.display()

# We calculate a variogram along the two main directions of the Model, i.e. 30 and 120 degrees.

varioparam = gl.VarioParam.createSeveral2D(angles=[30, 120], nlag=25, dlag=2)
myvario = gl.Vario.computeFromDb(varioparam, mypoints)
myvario.display()

########################
# fig, ax = plt.subplots()
# ax.varmod(myvario, model=mymodel, flagLegend=True)
# ax.decoration(title="Variogram in 2 Main Directions")
# fig.savefig('Figure2.png', dpi=600)
# plt.close(fig)
########################

# In the next step, we consider that the initial model (*mymodel*) has been correctly rendered through the simulated information at the 100 samples, as demonstrated in the experimental variogram. This is the reason why we will use this initial model for the further steps of estimation and conditional simulations.
#
# For this subsequent steps, we need to define the neighborhood which will serve in designating the set of samples, close to the target, which will be used for the estimation or simulation at the target location (called *myneigh*). Due to the small number of points, we decide to use the whole set of available samples, wherever the target is located: this is known as a *Unique* neighborhood.

myneigh = gl.NeighUnique.create()

# We have all the ingredients to perform the estimation using *Kriging*.
# When typing the name of the (output) grid data base, we can check that 2 variables have been added:
#
# - the one corresponding to the estimation result (called *Kriging.Simu.estim*)
# - the one corresponding to the standard deviation of the estimation error (called *kriging.Simu.stdev*)

err = gl.kriging(dbin=mypoints, dbout=mygrid, model=mymodel, neigh=myneigh)
mygrid.display()

########################
# cmap = gp.getColorMap(100,'YlOrBr')
# fig, ax = plt.subplots(1,2)
# ax[0].raster(mygrid, name="Kriging.Data.estim", flagLegend=True, vmin=-3.5, vmax=3.5)
# ax[0].symbol(mypoints, nameSize="Data", c="yellow", sizmin=1, sizmax=30)
# ax[0].decoration(title="Estimation", xlabel="Easting", ylabel="Northing")
# ax[1].raster(mygrid, name="Kriging.Data.stdev", cmap=cmap, flagLegend=True)
# ax[1].decoration(title="St. Dev. of Estimation Error", xlabel="Easting")
# ax[1].symbol(mypoints, nameSize="Data", flagCst=True, c="black")
# fig.tight_layout()
# fig.savefig('Figure3.png', dpi=600)
# plt.close(fig)
########################

# We can now construct 2 conditional simulations, using the Turning Bands algorithm again, but conditioned by the variable informed at the samples of *mypoints*. The newly created variables in *mygrid* are called *Simu.Data.1" and "Simu.Data.2". Each simulation outcome reproduces the spatial characteristics provided by the Model and honors the information provided at sample points.

err = gl.simtub(
    dbin=mypoints, dbout=mygrid, model=mymodel, neigh=myneigh, nbtuba=1000, nbsimu=2
)
mygrid.display()

########################
# fig, ax = plt.subplots(1,2)
# ax[0].raster(mygrid, name="Simu.Data.1", vmin=-3.5, vmax=3.5)
# ax[0].symbol(mypoints, nameSize="Data", flagCst=False, c="yellow", sizmin=1, sizmax=30)
# ax[0].decoration(title="Simulation 1", xlabel="Easting", ylabel="Northing")
# ax[1].raster(mygrid, name="Simu.Data.2", flagLegend=True, vmin=-3.5, vmax=3.5)
# ax[1].decoration(title="Simulation 2", xlabel="Easting")
# ax[1].symbol(mypoints, nameSize="Data", flagCst=False, c="yellow", sizmin=1, sizmax=30)
# fig.tight_layout()
# fig.savefig('Figure4.png', dpi=600)
# plt.close(fig)
########################
