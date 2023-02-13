Preamble
--------

Statistics
----------

<table>
<caption>Statistics on observations</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Mean</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Elevation</td>
<td align="right">87.97351</td>
</tr>
<tr class="even">
<td align="left">Temperature</td>
<td align="right">2.81457</td>
</tr>
</tbody>
</table>

<table>
<caption>Statistics on the grid</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Mean</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Elevation</td>
<td align="right">241.152</td>
</tr>
</tbody>
</table>

Observations
------------


Map of the elevations on Grid
-----------------------------


Temperature vs. Elevation
-------------------------


Temperature vs. Elevation (statistiques)
----------------------------------------

    tab <- dbStatisticsMultiT(dat, c("Temperature", "Elevation"),
                              oper=EStatOption_MEAN(), flagMono=FALSE)
    knitr::kable(tab$toTL(), caption = "Statistics on the grid", digits=3)

<table>
<caption>Statistics on the grid</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Temperature</th>
<th align="right">Elevation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Temperature</td>
<td align="right">2.815</td>
<td align="right">2.815</td>
</tr>
<tr class="even">
<td align="left">Elevation</td>
<td align="right">87.974</td>
<td align="right">146.441</td>
</tr>
</tbody>
</table>

Variograms and Non-stationarity
-------------------------------

Directional variograms (N-S and E-W): no drift or with first order
global trend

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


Global Trend
------------

Find the coefficients of the global regression (first order polynomial)

    err = regression(dat, name0="Temperature", mode=2, model=model, verbose=TRUE)

    ## 
    ## Linear Regression
    ## -----------------
    ## - Calculated on 151 active values
    ## - Explanatory Variable #1 = 3.521360
    ## - Explanatory Variable #2 = -0.007466
    ## - Explanatory Variable #3 = 0.001978
    ## - Initial variance        = 1.019788
    ## - Variance of residuals   = 0.735557

Modelling Variogram of Raw Data with Anisotropy
-----------------------------------------------

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


Ordinary Kriging - Cross-Validation
-----------------------------------

Perform Cross-validation (using Ordinary Kriging) and calculate the Mean
Squared Error.

    err = xvalid(dat, fitmod, unique_neigh,
                 namconv=NamingConvention("OK",TRUE,TRUE,FALSE))

<table>
<caption>Cross-validation with Ordinary Kriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Mean</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.esterr</td>
<td align="right">-0.0406019</td>
</tr>
<tr class="even">
<td align="left">OK.Temperature.stderr</td>
<td align="right">-7.9338923</td>
</tr>
</tbody>
</table>

Ordinary Kriging - Estimation
-----------------------------

    err = kriging(dat, target, fitmod, unique_neigh,
                  namconv=NamingConvention("OK"))


Ordinary Kriging - Statistics
-----------------------------

<table>
<caption>Statistics on the Ordinary Kriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">0.599</td>
<td align="right">5.084</td>
<td align="right">2.806</td>
<td align="right">0.910</td>
</tr>
<tr class="even">
<td align="left">OK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.064</td>
<td align="right">0.883</td>
<td align="right">0.460</td>
<td align="right">0.127</td>
</tr>
</tbody>
</table>

Fitting Variogram of Residuals
------------------------------

    fitmodUK = Model()
    err = fitmodUK$fit(vario_res2dir,
                       types=ECov_fromKeys(c("SPHERICAL")),
                       optvar=Option_VarioFit(FALSE,FALSE))

    p = plot.varmod(vario_res2dir, fitmodUK)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    plot.decoration(p, title="Temperature (°C)")


Note that the residuals seem isotropic, hence use isotropic option for
Fitting.

Universal Kriging - Cross-Validation
------------------------------------

Perform Cross-validation (using Universal Kriging) and calculate the
Mean Squared Error.

    err = xvalid(dat, fitmodUK, unique_neigh,
                 namconv=NamingConvention("UK",TRUE,TRUE,FALSE))

<table>
<caption>Cross-validation with Ordinary Kriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">UK.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-4.336</td>
<td align="right">1.397</td>
<td align="right">-0.352</td>
<td align="right">0.939</td>
</tr>
<tr class="even">
<td align="left">UK.Temperature.stderr</td>
<td align="right">151</td>
<td align="right">-47.724</td>
<td align="right">1.075</td>
<td align="right">-8.056</td>
<td align="right">6.170</td>
</tr>
</tbody>
</table>

Universal Kriging - Estimation
------------------------------

    err = kriging(dat, target, fitmodUK, unique_neigh,
                  namconv=NamingConvention("UK"))


Universal Kriging - Statistics
------------------------------

<table>
<caption>Statistics on the Universal Kriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">UK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">0.471</td>
<td align="right">4.902</td>
<td align="right">2.514</td>
<td align="right">0.821</td>
</tr>
<tr class="even">
<td align="left">UK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.069</td>
<td align="right">0.891</td>
<td align="right">0.482</td>
<td align="right">0.140</td>
</tr>
</tbody>
</table>

Comparing Ordinary and Universal Krigings
-----------------------------------------

    p = plot.correlation(target, "OK*estim", "UK*estim", flagDiag=TRUE, bins=100)
    plot.decoration(p, xlab="Ordinary Kriging",ylab="Universal Kriging")


Comparing Temperature and Elevation
-----------------------------------

    plot.correlation(dat,"Elevation","Temperature",flagRegr=TRUE,asPoint=TRUE)


Comparing Temperature and Elevation
-----------------------------------

    dat$setLocators(c("Temperature", "Elevation"), ELoc_Z())

    ## NULL

Bivariate Modelling
-------------------

    varioexp2var = Vario_create(varioparam, dat)
    err = varioexp2var$compute()
    fitmod2var = Model()
    err = fitmod2var$fit(varioexp2var,
                         types=ECov_fromKeys(c("NUGGET","EXPONENTIAL","CUBIC")))

    p = plot.varmod(varioexp2var, fitmod2var)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    plot.decoration(p, title="Temperature (°C) and Elevation")


Cokriging with elevation - Cross-Validation
-------------------------------------------

Most of the processes are more time-consuming in Unique Neighborhood. We
create a small neighborhood for demonstration.

    moving_neigh = NeighMoving_create(radius = 1000, nmaxi = 10)

Perform Cross-validation (Bivariate Model) and calculate the Mean
Squared Error.

    err = xvalid(dat, fitmod2var, moving_neigh,
                 namconv=NamingConvention("COK",TRUE,TRUE,FALSE))

<table>
<caption>Cross-validation with Cokriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">COK.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-4.582</td>
<td align="right">1.435</td>
<td align="right">-0.407</td>
<td align="right">0.926</td>
</tr>
<tr class="even">
<td align="left">COK.Temperature.stderr</td>
<td align="right">151</td>
<td align="right">-55.040</td>
<td align="right">1.337</td>
<td align="right">-8.735</td>
<td align="right">7.099</td>
</tr>
</tbody>
</table>

Cokriging with elevation - Estimate
-----------------------------------

    err = kriging(dat, target, fitmod2var, unique_neigh,
                  namconv=NamingConvention("COK"))


Cokriging with elevation - Statistics
-------------------------------------

<table>
<caption>Statistics on Cokriging</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">COK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-0.042</td>
<td align="right">4.920</td>
<td align="right">2.553</td>
<td align="right">0.911</td>
</tr>
<tr class="even">
<td align="left">COK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.063</td>
<td align="right">1.422</td>
<td align="right">0.388</td>
<td align="right">0.136</td>
</tr>
</tbody>
</table>

Comparing Kriging and CoKriging
-------------------------------


Note that CoKriging produces estimates which are mostly larger than
Kriging estimates.

Kriging the residuals
---------------------

*Z*<sub>2</sub>(*s*)=*b* + *a**Z*<sub>1</sub>(*s*)+*R*(*s*)

    regr = regression(dat, "Temperature", "Elevation", flagCste=TRUE, verbose=TRUE)

    ## 
    ## Linear Regression
    ## -----------------
    ## - Calculated on 151 active values
    ## - Constant term           = 3.611970
    ## - Explanatory Variable #1 = -0.009064
    ## - Initial variance        = 1.019788
    ## - Variance of residuals   = 0.363298

    b = regr$coeffs[1]
    a = regr$coeffs[2]

    err = dbRegression(dat, "Temperature", "Elevation",
                         namconv = NamingConvention("Regr",TRUE,TRUE,FALSE))

<table>
<caption>Statistics on the residual</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Regr.Temperature</td>
<td align="right">151</td>
<td align="right">-1.359</td>
<td align="right">1.795</td>
<td align="right">0</td>
<td align="right">0.603</td>
</tr>
</tbody>
</table>

Kriging the residuals - Correlation
-----------------------------------

    plot.correlation(dat,"Elevation","Regr*",flagRegr=TRUE,asPoint=TRUE)


Kriging the residuals - Variogram of the residual
-------------------------------------------------

    dat$setLocator("Regr*",ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    varioexpR = Vario(varioparam, dat)
    err = varioexpR$compute()
    fitmodR = Model()
    err = fitmodR$fit(varioexpR,
                      types=ECov_fromKeys(c("NUGGET","SPHERICAL","ORDER1_GC")))

    p = plot.varmod(varioexpR, fitmodR)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    plot.decoration(p, title="Temperature Residual")


Kriging the residuals - Variogram of the residual
-------------------------------------------------

    err = kriging(dat, target, fitmodR, unique_neigh,
                  namconv=NamingConvention("ROK"))


Kriging the residuals - Computing the estimate
----------------------------------------------

*Z*<sub>2</sub><sup>⋆</sup> = *b* + *a**Z*<sub>1</sub>(*s*)+*R*(*s*)<sup>*O**K*</sup>

    ROK_estim = target["ROK.Regr*estim"] + b + a * target["Elevation"]
    uid = target$addColumns(ROK_estim,"KR.Temperature.estim")


Kriging the residuals - Correlation
-----------------------------------

Correlation between Ordinary Kriging and CoKriging


Kriging the residuals - Correlation
-----------------------------------

Correlation between Ordinary Kriging and Kriging with Residuals


Kriging the residuals - Correlation
-----------------------------------

<table>
<caption>Statistics on OK and Kriging of the resisual</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">0.599</td>
<td align="right">5.084</td>
<td align="right">2.806</td>
<td align="right">0.910</td>
</tr>
<tr class="even">
<td align="left">COK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-0.042</td>
<td align="right">4.920</td>
<td align="right">2.553</td>
<td align="right">0.911</td>
</tr>
<tr class="odd">
<td align="left">KR.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-8.097</td>
<td align="right">5.108</td>
<td align="right">1.445</td>
<td align="right">1.906</td>
</tr>
</tbody>
</table>

Using Elevation Map as External Drift
-------------------------------------

Preparing the data bases

    err = dat$setLocator("Temperature",ELoc_Z(),cleanSameLocator=TRUE)
    err = dat$setLocator("Elevation",ELoc_F(),cleanSameLocator=TRUE)

Linear Regression of Temperature knowing Elevation on Data
----------------------------------------------------------

We perform the Linear Regression of Temperature knowing Elevation
(specified as External Drift)

    regr = regression(dat, "Temperature", mode=1, flagCste=TRUE, verbose=TRUE)

    ## 
    ## Linear Regression
    ## -----------------
    ## - Calculated on 151 active values
    ## - Constant term           = 3.611970
    ## - Explanatory Variable #1 = -0.009064
    ## - Initial variance        = 1.019788
    ## - Variance of residuals   = 0.363298

Experimental Variogram of Residuals (External Drift)
----------------------------------------------------

    varioKED = Vario(varioparam, dat)

    model = Model_create()
    err = model$setDriftIRF(order=0, nfex=1)
    err = varioKED$compute(model=model)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


Model of Residuals (External Drift)
-----------------------------------

    modelKED = Model()
    err = modelKED$fit(varioKED,
                       types=ECov_fromKeys(c("NUGGET","CUBIC","GAUSSIAN")))
    modelKED$setDriftIRF(order = 0, nfex = 1)

    ## NULL

    plot.varmod(varioKED, modelKED)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


Kriging with External Drift - Cross-Validation
----------------------------------------------

Perform Cross-validation (External Drift) and calculate the Mean Squared
Error.

    err = xvalid(dat, modelKED, unique_neigh,
                 namconv=NamingConvention("KED",TRUE,TRUE,FALSE))

<table>
<caption>Cross-validation with Kriging with External Drift</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">KED.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-1.577</td>
<td align="right">1.001</td>
<td align="right">-0.009</td>
<td align="right">0.414</td>
</tr>
<tr class="even">
<td align="left">KED.Temperature.stderr</td>
<td align="right">151</td>
<td align="right">-15.061</td>
<td align="right">0.441</td>
<td align="right">-7.569</td>
<td align="right">3.143</td>
</tr>
</tbody>
</table>

Kriging with External Drift - Estimate
--------------------------------------

    err = kriging(dat, target, modelKED, unique_neigh,
                  namconv=NamingConvention("KED"))


Kriging with External Drift - Statistics
----------------------------------------

<table>
<caption>Statistics on the Kriging with External Drift</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">KED.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-6.0035163</td>
<td align="right">4.7725035</td>
<td align="right">1.7779984</td>
<td align="right">1.5399727</td>
</tr>
<tr class="even">
<td align="left">KED.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.3117728</td>
<td align="right">0.6147864</td>
<td align="right">0.3958358</td>
<td align="right">0.0508596</td>
</tr>
</tbody>
</table>

Comparing OK with KED
---------------------


Note that negative Estimates are present when using External Drift.

Summary of Cross-validation scores
----------------------------------

Comparing the cross-validation Mean Squared Errors

    tab <- dbStatisticsMonoT(dat, opers=opers, names = (c("*.Temperature.esterr")))
    knitr::kable(tab$toTL(),
      caption = "Statistics of the cross-validation error for the temperature"
      )

<table>
<caption>Statistics of the cross-validation error for the temperature</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-2.126215</td>
<td align="right">1.851849</td>
<td align="right">-0.0406019</td>
<td align="right">0.5605902</td>
</tr>
<tr class="even">
<td align="left">UK.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-4.336191</td>
<td align="right">1.397189</td>
<td align="right">-0.3520194</td>
<td align="right">0.9386221</td>
</tr>
<tr class="odd">
<td align="left">COK.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-4.582072</td>
<td align="right">1.434923</td>
<td align="right">-0.4074638</td>
<td align="right">0.9255625</td>
</tr>
<tr class="even">
<td align="left">KED.Temperature.esterr</td>
<td align="right">151</td>
<td align="right">-1.576803</td>
<td align="right">1.001177</td>
<td align="right">-0.0092625</td>
<td align="right">0.4144602</td>
</tr>
</tbody>
</table>

Statistics on the estimates
---------------------------

Computing the statistics on the various estimates.

    tab <- dbStatisticsMonoT(target, opers=opers, 
                             names = (c("*.Temperature.estim")))
    knitr::kable(tab$toTL(), digits=3,
      caption = "Statistics of the estimates for the temperature"
      )

<table>
<caption>Statistics of the estimates for the temperature</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">0.599</td>
<td align="right">5.084</td>
<td align="right">2.806</td>
<td align="right">0.910</td>
</tr>
<tr class="even">
<td align="left">UK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">0.471</td>
<td align="right">4.902</td>
<td align="right">2.514</td>
<td align="right">0.821</td>
</tr>
<tr class="odd">
<td align="left">COK.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-0.042</td>
<td align="right">4.920</td>
<td align="right">2.553</td>
<td align="right">0.911</td>
</tr>
<tr class="even">
<td align="left">ROK.Regr.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-0.771</td>
<td align="right">1.586</td>
<td align="right">0.019</td>
<td align="right">0.455</td>
</tr>
<tr class="odd">
<td align="left">KR.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-8.097</td>
<td align="right">5.108</td>
<td align="right">1.445</td>
<td align="right">1.906</td>
</tr>
<tr class="even">
<td align="left">KED.Temperature.estim</td>
<td align="right">3092</td>
<td align="right">-6.004</td>
<td align="right">4.773</td>
<td align="right">1.778</td>
<td align="right">1.540</td>
</tr>
</tbody>
</table>

Statistics on the estimates
---------------------------

Computing the statistics on the various estimates.

    tab <- dbStatisticsMonoT(target, opers=opers, 
                             names = (c("*.Temperature.stdev")))
    knitr::kable(tab$toTL(), digits=3,
      caption = "Statistics of the kriging Std. for the temperature"
      )

<table>
<caption>Statistics of the kriging Std. for the temperature</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Count</th>
<th align="right">Minimum</th>
<th align="right">Maximum</th>
<th align="right">Mean</th>
<th align="right">St. Dev.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">OK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.064</td>
<td align="right">0.883</td>
<td align="right">0.460</td>
<td align="right">0.127</td>
</tr>
<tr class="even">
<td align="left">UK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.069</td>
<td align="right">0.891</td>
<td align="right">0.482</td>
<td align="right">0.140</td>
</tr>
<tr class="odd">
<td align="left">COK.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.063</td>
<td align="right">1.422</td>
<td align="right">0.388</td>
<td align="right">0.136</td>
</tr>
<tr class="even">
<td align="left">ROK.Regr.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.304</td>
<td align="right">0.504</td>
<td align="right">0.362</td>
<td align="right">0.031</td>
</tr>
<tr class="odd">
<td align="left">KED.Temperature.stdev</td>
<td align="right">3092</td>
<td align="right">0.312</td>
<td align="right">0.615</td>
<td align="right">0.396</td>
<td align="right">0.051</td>
</tr>
</tbody>
</table>


