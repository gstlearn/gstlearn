Introduction
============

This case study is meant to demonstrate how to use gstlearn for
non-linear geostatistics with the Gaussian model.

We consider three supports:

-   The support of the samples, considered as **points** and noted *x*,

-   The support of the selection, the **bloc**, noted *v*,

-   The support of the reporting, the **panels**, noted *V*.

The variable of interest *z* is defined over the domain
*S* ⊂ ℝ<sup>*d*</sup> with *d* = 2. The domain is uniformly divided into
panels and each panel is regularly subdivided into blocs, hence defining
two regular grids, the grid of panels and the grid of blocs.

> **Hence, we will have a data set and two grids.**

The regionalised variable *z*(*x*) for *x* ∈ *D* ⊂ ℝ<sup>*d*</sup> is
modeled as a transform of the stationary Gaussian Random Function,
*Z*(*x*)=*ϕ*(*Y*(*x*)) where

-   *Y*(*x*) is a stationary Gaussian function, centered and normalized.

> *E*{*Y*(*x*)} = 0 and
> *C**o**v*{*Y*(*x*),*Y*(*x* + *h*)} = *ρ*<sub>*Y*</sub>(*h*)=1 − *γ*<sub>*Y*</sub>(*h*)

-   *ϕ* is a continous and one to one and mapping, called the Gaussian
    anamorphosis.

> As *V**a**r*(*Z*)&lt; + ∞, *ϕ* ∈ *L*<sup>2</sup>(*g*) and it can be
> expressed as a linear combination of Hermite polynomials
> $\\phi(y) = \\sum\_{n = 1}^{+\\infty} \\phi\_n H\_n(y)$.

Non linear Geostatistics implements non linear estimators of non linear
transforms of the variable of interest *Z*. Two issues are addressed,
the selection and the change of support, with the following tasks common
in geosciences (i.e., in a mining context or for environmental studies),

-   predicting recovered quantities with selection over the actual
    value.

> For example, recovered mineral resources will be defined by the
> selected ore at a given cutoff,
> *T*(*z*<sub>*c*</sub>)=1<sub>*Z* ≥ *z*<sub>*c*</sub></sub>, and the
> metal contained in the selected ore,
> *Q*(*z*<sub>*c*</sub>)=*Z* × 1<sub>*Z* ≥ *z*<sub>*c*</sub></sub>.

-   taking into account the support of the selection.

> The average value of the volume *v* is noted
> $Z(v) = \\frac{1}{|v|} \\int\_v Z(u) du$.

A first task is to predict marginal distribution of *Z*(*v*), or the
histogram, knowing the histogram of *Z* and its spatial structure. A
second task is to predict the conditional distribution of *Z*(*v*) given
the spatial the prior spatial model and some observations. The first
question is referred to as the global recoverable resources, and the
second one as the local recoverable resources.

Three estimators will be illustrated:

-   the conditional expectation (**EC**)

-   the disjunctive kriging (**DK**)

> EC and DK can be used to evaluate **point** recovery functions at a
> non observed point *x*<sub>0</sub>,
> 1<sub>*Z*(*x*<sub>0</sub>)≥*z*<sub>*c*</sub></sub> and
> *Z*(*x*<sub>0</sub>)×1<sub>*Z*(*x*<sub>0</sub>)≥*z*<sub>*c*</sub></sub>.
> They can also be used to evaluate **bloc** recovery functions for a
> block *v*, 1<sub>*Z*(*v*)≥*z*<sub>*c*</sub></sub> and
> *Z*(*v*)×1<sub>*Z*(*v*)≥*z*<sub>*c*</sub></sub>. DK can also evaluate
> recovery functions average on a bigger support, e.g.
> $\\frac{1}{N} \\sum\_{i=1}^{N} 1\_{Z(v\_i) \\geq z\_c}$ is the
> recovered ore on the panels $V = \\cup\_{i=1}^{N}v\_i$.

-   the uniform conditioning (**UC**)

> UC computes block recovery functions averaged on a panel.

Two change of support models are available for a Gaussian model, they
are noted respectively *D**G**M* − 1 and *D**G**M* − 2.

Initially, we generate a realization of the variable on a fine grid,
modeling the regionalized variable defined on the point support. This
realization, or simulation, is sampled to defined (*n**p* points
uniformly sampled on the domain) the input data set.

Initialisation
--------------

We will use the Geostatistics library **gstlearn**.

    # Set the Seed for the Random Number generator
    law_set_random_seed(32131)

    ## NULL

    # Defining the format for dumping statistics of variables in a Db
    dbfmt = DbStringFormat()
    dbfmt$setFlags(flag_resume=TRUE,flag_vars=FALSE,flag_locator=TRUE)

    ## NULL

    # Defining global options
    flag.debug = FALSE

Setting the trace: this option allows dumping all the elements to check
calculations when the rank of the target is equal to 1 (first block,
first panel, ...). This option is important for debugging but it creates
lots of printout.

    OptDbg_reset()

    ## NULL

    if (flag.debug) OptDbg_setReference(1)

Defining color scales common to all future representations of estimated
quantities:

    ZEstMin = -3.
    ZEstMax = +3.
    TEstMin = 0.
    TEstMax = 1.
    TStdMin = 0.
    TStdMax = 1.
    QEstMin = 0.
    QEstMax = 3.
    QStdMin = 0.
    QStdMax = 2.

Generate initial grids
----------------------

Three grids are defined:

-   The grid of the samples. It representes a reference realization of
    the regionalized variable. This realization is sampled to define the
    data set of the observations. This data set is then used to evaluate
    the ability of the non linear technique to reproduce the selectivity
    curves computed on the reference simulation.

-   The grid of the panels
-   The grid of the blocs

<!-- -->

    # grid of samples
    nx_S = c(100,100)
    dx_S = c(0.01, 0.01)
    # grid of panels
    dx_P = 0.25
    # grid of blocs
    nx_B = 5
    dx_B = dx_P / nx_B

    # Generate initial grid
    grid = DbGrid_create(nx = nx_S, dx = dx_S)
    grid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 3
    ## Maximum Number of UIDs       = 3
    ## Total number of samples      = 10000
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.010     0.010
    ## Number :        100       100
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2

    ## NULL

    # Create grid of panels covering the simulated area
    panels = DbGrid_createCoveringDb(grid, dcell=c(dx_P,dx_P))
    panels$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 2
    ## Maximum Number of UIDs       = 2
    ## Total number of samples      = 25
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.250     0.250
    ## Number :          5         5
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2

    ## NULL

    # Discretization with a grid of blocks which covers the simulated area
    blocs = DbGrid_createCoveringDb(grid, dcell=c(dx_B,dx_B))
    blocs$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 2
    ## Maximum Number of UIDs       = 2
    ## Total number of samples      = 441
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.050     0.050
    ## Number :         21        21
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2

    ## NULL

Simulation of the data set
--------------------------

A lognormal model is defined and a simulation is performed on the grid
of samples.

The regionalized variable *z*(*x*) is modeled using a lognormal model
*Z*(*x*)=*m* × *e*<sup>*σ* *Y*(*x*)−1/2*σ*<sup>2</sup></sup> with

-   *Y* a stationary Gaussian Random Function with an exponential
    variogram with a range equal to 0.1 and a sill equal to 1.0. The
    mean of *Y* is null.

-   The lognormal transform is parametrized by its mean *m* and the
    dispersion coefficient *σ*. The first two moments are:
-   Mean *E*{*Z*}=*m*
-   Variance
    *V**a**r*{*Z*}=*m*<sup>2</sup>(*e*<sup>*σ*<sup>2</sup></sup> − 1)

<!-- -->

    # Simulation of the Gaussian variable
    model_init = Model_createFromParam(ECov_EXPONENTIAL(), range=0.1, sill=1.)
    err = simtub(NULL, grid, model_init, namconv=NamingConvention("Y"))
    # Nonlinear transform (lognormal)
    m_Z = 1.5
    s_Z = 0.5
    grid["Z"] = m_Z * exp(s_Z * grid["Y"] - 0.5*s_Z^2)

    plot.grid(grid, "Z")


    p = plot.hist(grid, name = "Y", nbins=100, col='gray', fill='skyblue')
    plot.decoration(p, xlab = "Y", title = "Simulated samples (Y variable)")


    p = plot.hist(grid, name = "Z", nbins = 100, col = "gray", fill = "orange")
    plot.decoration(p, xlab = "Z", title = "Simulated samples (Z variable)")


    opers = EStatOption_fromKeys(c("NUM","MINI","MAXI","MEAN","STDV"))
    stats = dbStatisticsMonoT(grid, c("Y","Z"), opers=opers)
    knitr::kable(stats$toTL(), digits=2,
        caption = "Statistics of the simulated variables on the point support")

<table>
<caption>Statistics of the simulated variables on the point support</caption>
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
<td align="left">Y</td>
<td align="right">10000</td>
<td align="right">-4.30</td>
<td align="right">3.65</td>
<td align="right">0.07</td>
<td align="right">1.01</td>
</tr>
<tr class="even">
<td align="left">Z</td>
<td align="right">10000</td>
<td align="right">0.15</td>
<td align="right">8.20</td>
<td align="right">1.56</td>
<td align="right">0.85</td>
</tr>
</tbody>
</table>

Data extraction
---------------

We create a new data set by extracting few samples from the previous
grid.

    np = 500
    data = Db_createSamplingDb(grid, number=np, names=c("x1","x2","Y","Z"))
    data$setLocator("Z", ELoc_Z())

    ## NULL

    data$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 4
    ## Maximum Number of UIDs       = 4
    ## Total number of samples      = 500
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Y - Locator = NA
    ## Column = 3 - Name = Z - Locator = z1

    ## NULL

    plot.point(data, name_color = "Z", name_size = "Z", show.legend.symbol = TRUE)


    p = plot.hist(data, name = "Y", nbins = 100, col = "gray", fill = "skyblue")
    plot.decoration(p, xlab = "Y", title = "Sampled data (Y variable)")


    p = plot.hist(data, name = "Z", nbins = 100, col = "gray", fill = "orange")
    plot.decoration(p, xlab = "Z", title = "Sampled data (Z variable)")


    stats = dbStatisticsMonoT(data, c("Y","Z"), opers=opers)
    knitr::kable(stats$toTL(), digits=2, 
                 caption = "Statistics of the sampled data on the point support")

<table>
<caption>Statistics of the sampled data on the point support</caption>
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
<td align="left">Y</td>
<td align="right">500</td>
<td align="right">-3.42</td>
<td align="right">2.99</td>
<td align="right">0.03</td>
<td align="right">1.04</td>
</tr>
<tr class="even">
<td align="left">Z</td>
<td align="right">500</td>
<td align="right">0.24</td>
<td align="right">5.92</td>
<td align="right">1.54</td>
<td align="right">0.84</td>
</tr>
</tbody>
</table>

    varZ <- stats$getValue(1,4)^2

Gaussian Anamorphosis with 20 coefficients
------------------------------------------

    anam = AnamHermite_create(nbpoly=20)
    err = anam$fit(data, "Z")
    anam

    ## 
    ## Hermitian Anamorphosis
    ## ----------------------
    ## Minimum absolute value for Y  = -2.8
    ## Maximum absolute value for Y  = 3.9
    ## Minimum absolute value for Z  = 0.274228
    ## Maximum absolute value for Z  = 5.99041
    ## Minimum practical value for Y = -2.8
    ## Maximum practical value for Y = 3.1
    ## Minimum practical value for Z = 0.274228
    ## Maximum practical value for Z = 5.86055
    ## Mean                          = 1.53861
    ## Variance                      = 0.709613
    ## Number of Hermite polynomials = 20
    ## Normalized coefficients for Hermite polynomials (punctual variable)
    ##                [,  0]    [,  1]    [,  2]    [,  3]    [,  4]    [,  5]    [,  6]
    ##      [  0,]     1.539    -0.792     0.273    -0.061    -0.009     0.032    -0.031
    ##      [  7,]     0.010     0.018    -0.021    -0.001     0.016    -0.009    -0.007
    ##      [ 14,]     0.012     0.000    -0.012     0.004     0.010    -0.006

Selectivity curves
------------------

We focus on Tonnage (T) and Metal Quantity (Q) for few cutoffs (Zcuts).

    Zcuts <- c(0.0, 0.5, 0.75, 1.0)
    selectivity = Selectivity_createByKeys(
      keys = c("T", "Q"), zcuts=Zcuts,
      flag_est=TRUE, flag_std=TRUE)

    # Global experimental selectivity, calculated form the experimental Data Set
    knitr::kable(selectivity$eval(data)$toTL(), digits=3,
                 caption = "Selectivity curves computed on data set")

<table>
<caption>Selectivity curves computed on data set</caption>
<thead>
<tr class="header">
<th align="right">Z-Cut</th>
<th align="right">T-estim</th>
<th align="right">Q-estim</th>
<th align="right">B-estim</th>
<th align="right">M-estim</th>
<th align="right">T-stdev</th>
<th align="right">Q-stdev</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">0.00</td>
<td align="right">1.000</td>
<td align="right">1.539</td>
<td align="right">1.539</td>
<td align="right">1.539</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="right">0.50</td>
<td align="right">0.972</td>
<td align="right">1.527</td>
<td align="right">1.041</td>
<td align="right">1.571</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="odd">
<td align="right">0.75</td>
<td align="right">0.858</td>
<td align="right">1.453</td>
<td align="right">0.810</td>
<td align="right">1.694</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="right">1.00</td>
<td align="right">0.718</td>
<td align="right">1.330</td>
<td align="right">0.612</td>
<td align="right">1.853</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
</tbody>
</table>

    # Selectivity in the model, derived from the parameters contained in the Anamorphosis
    knitr::kable(selectivity$evalFromAnamorphosis(anam)$toTL(), digits=3,
                 caption = "Selectivity curves computed on anamorphosis")

<table>
<caption>Selectivity curves computed on anamorphosis</caption>
<thead>
<tr class="header">
<th align="right">Z-Cut</th>
<th align="right">T-estim</th>
<th align="right">Q-estim</th>
<th align="right">B-estim</th>
<th align="right">M-estim</th>
<th align="right">T-stdev</th>
<th align="right">Q-stdev</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">0.00</td>
<td align="right">1.000</td>
<td align="right">1.539</td>
<td align="right">1.539</td>
<td align="right">1.539</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="right">0.50</td>
<td align="right">0.972</td>
<td align="right">1.527</td>
<td align="right">1.041</td>
<td align="right">1.570</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="odd">
<td align="right">0.75</td>
<td align="right">0.867</td>
<td align="right">1.461</td>
<td align="right">0.810</td>
<td align="right">1.684</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="right">1.00</td>
<td align="right">0.716</td>
<td align="right">1.327</td>
<td align="right">0.611</td>
<td align="right">1.853</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
</tbody>
</table>

Transform Data into Gaussian variable
=====================================

    data["Gaussian.Z"] <- VectorHelper_normalScore(data["Z"])
    p = plot.hist(data, name = "Gaussian.Z", nbins = 25, 
                  col = "gray", fill = "yellow") 
    plot.decoration(p, xlab = "Normal score")


    data$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 5
    ## Maximum Number of UIDs       = 5
    ## Total number of samples      = 500
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Y - Locator = NA
    ## Column = 3 - Name = Z - Locator = z1
    ## Column = 4 - Name = Gaussian.Z - Locator = NA

    ## NULL

Variography
===========

Define the variogram calculation parameters:

-   omni-directional variogram,

-   10 lags of 0.025.

<!-- -->

    varioparam = VarioParam_createOmniDirection(npas=10, dpas=0.025)

Variography of the raw variable
-------------------------------

Define the variogram calculation parameters:

    # Computing the experimental variogram
    err = data$setLocator("Z", ELoc_Z())
    vario_raw = Vario_computeFromDb(varioparam, db=data)

    # Fitting the variogram model on the experimental variogram
    model_raw = Model_create()
    err = model_raw$fit(vario_raw, 
                        types = ECov_fromKeys(c("NUGGET", "EXPONENTIAL", "EXPONENTIAL")), 
                        constraints=Constraints(varZ))
    model_raw$display()

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 2
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Nugget Effect
    ## - Sill         =      0.066
    ## Exponential
    ## - Sill         =      0.645
    ## - Range        =      0.117
    ## - Theo. Range  =      0.039
    ## Total Sill     =      0.711

    ## NULL

    p = plot.varmod(vario_raw, model_raw)
    plot.decoration(p, title = "Raw variable")


Variography of the Gaussian variable
====================================

    # Computing of the experimental variogram
    err = data$setLocator("Gaussian.Z", ELoc_Z(), cleanSameLocator=TRUE)
    vario = Vario_computeFromDb(varioparam, data)

    # Fitting the Model on the experimental variogram with a sill equal to one.
    model = Model()
    err = model$fit(vario, 
                    types = ECov_fromKeys(c("EXPONENTIAL", "EXPONENTIAL")), 
                    constraints=Constraints(1.))
    model$display()

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 2
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Exponential
    ## - Sill         =      0.220
    ## - Range        =      0.013
    ## - Theo. Range  =      0.004
    ## Exponential
    ## - Sill         =      0.780
    ## - Range        =      0.137
    ## - Theo. Range  =      0.046
    ## Total Sill     =      1.000

    ## NULL

    p = plot.varmod(vario, model)
    plot.decoration(p, title = "Gaussian variable")


    model_Y = model$clone()

    err = model$setAnam(anam)

Checking the experimental variogram of the Raw variable against its
Model derived from the Model of the Gaussian transform.

    model$setActiveFactor(-1)

    ## NULL

    p = plot.varmod(vario_raw, model)
    plot.decoration(p, title = "Raw variable (model of the Gaussian)")


    model$setActiveFactor(0)

    ## NULL

Creating a Moving Neighborhood
------------------------------

    nmini = 5
    nmaxi = 5
    radius = 1.
    neigh = NeighMoving_create(nmaxi=nmaxi, radius=radius, nmini=nmini)
    neigh

    ## 
    ## Moving Neighborhood
    ## ===================
    ## Space dimension = 2
    ## Minimum number of samples           = 5
    ## Maximum number of samples           = 5
    ## Maximum horizontal distance         = 1

Nonlinear estimates with the Gaussian model
===========================================

A unified workflow for the nonlinear techniques available with the
Gaussian model is proposed below. The estimators are: the conditional
expectation (CE), the disjunctive kriging (DK), and the uniform
conditioning (UC).

The workflow is:

1.  Pre-processing of the input dataset using an anamorphosis. It is
    required for CE and DK.

-   CE: computation of the Gaussian values from the raw values

-   DK: computation of the values of the Hermite polynomials from the
    raw values

-   UC: No pre-processing is required.

1.  Computing the kriging

-   CE: computation of simple kriging of the Gaussian variable
    (estimation value and standard deviation). The standard function
    **kriging** is used. Ordinary kriging may be used.

-   DK: computation of the kriging of the Hermite polynomials from the
    raw values. The function **krigingFactors** is used to loop over the
    factors and define the proper covariance function.

-   UC: computation of ordinary/simple kriging of the raw variable
    (estimation value and variance of the estimator). The standard
    function **kriging** is used.

Here, we have

**kriging(data, target, model, ....)**, or **krigingFactors(data,
target, model, ...)**

1.  Post-processing to calculate the estimator of the selectivity curves
    or any nonlinear transform from the kriging of the factors or
    Gaussian variable.

-   CE: dedicated function **ConditionExpectation(db, anam, selectivity,
    ...)**

-   DK: dedicated function **DisjunctiveKriging(db, anam, selectivity,
    ...)**

-   UC: dedicated function **UniformConditioning(db, anam, selectivity,
    ...)**

Two versions are now considered: the point estimate and the block
estimate.

Point estimate
==============

Two methods are presented:

-   The conditional expectation

-   The disjunctive kriging

Conditional Expectation
-----------------------

The Conditional Expectation is used to estimate **selectivity curves**
for a point support at the grid nodes (i.e. the center of the blocs).

### Preprocessing

The Gaussian variable has been computed earlier.

### Kriging of the Gaussian variable

Simple Kriging is used to estimate the point value of the Gaussian
variable at the grid nodes (point support). The known mean is equal to
zero.

    err = model$delDrift(0)
    err = model$setMean(0, 0.)
    data$setLocator("Gaussian.Z", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    err = kriging(data, blocs, model, neigh, 
                  calcul = EKrigOpt_PONCTUAL(),
                  namconv= NamingConvention("G_Pts"))

    p = plot(blocs, "G_Pts.Gaussian.Z.estim", 
             show.legend.raster=TRUE, legend.name.raster = "Estimation")
    plot.decoration(p, title = "Punctual Simple Kriging of Y")


    p = plot(blocs, "G_Pts.Gaussian.Z.stdev", 
             show.legend.raster=TRUE, legend.name.raster = "Std. Dev.")
    plot.decoration(p, title = "Punctual Simple Kriging of Y")


### Calculating the conditional expectation on blocs

Kriging has been computed, the back transform is computed for a punctual
suport:

*ψ*(*Y*(*x*<sub>0</sub>))<sup>*C**E*</sup> = ∫*ψ*(*y*(*x*<sub>0</sub>)<sup>*S**K*</sup> + *σ*<sub>*S**K*</sub> *u*)*g*(*u*)*d**u*
 It is possible to check the computed values for ore:

$$
(1\_{Y(x\_0)\\geq y\_c})^{CE} = 1 - G(\\frac{y\_c - y(x\_0)^{SK}}{\\sigma\_{SK}})
$$

For change of support calculation DGM kriging (**krigdgm**) should be
used. **To be checked ?**

    err = ConditionalExpectation(db = blocs, 
                                 anam = anam, 
                                 selectivity = selectivity, 
                                 name_est="G_Pts*estim", 
                                 name_std="G_Pts*stdev",
                                 namconv = NamingConvention("Pts_Recovery", FALSE))

    blocs$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 20
    ## Maximum Number of UIDs       = 20
    ## Total number of samples      = 441
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.050     0.050
    ## Number :         21        21
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = G_Pts.Gaussian.Z.estim - Locator = NA
    ## Column = 3 - Name = G_Pts.Gaussian.Z.stdev - Locator = NA
    ## Column = 4 - Name = Pts_Recovery.T-estim-0 - Locator = NA
    ## Column = 5 - Name = Pts_Recovery.T-estim-0.5 - Locator = NA
    ## Column = 6 - Name = Pts_Recovery.T-estim-0.75 - Locator = NA
    ## Column = 7 - Name = Pts_Recovery.T-estim-1 - Locator = NA
    ## Column = 8 - Name = Pts_Recovery.T-stdev-0 - Locator = NA
    ## Column = 9 - Name = Pts_Recovery.T-stdev-0.5 - Locator = NA
    ## Column = 10 - Name = Pts_Recovery.T-stdev-0.75 - Locator = NA
    ## Column = 11 - Name = Pts_Recovery.T-stdev-1 - Locator = NA
    ## Column = 12 - Name = Pts_Recovery.Q-estim-0 - Locator = NA
    ## Column = 13 - Name = Pts_Recovery.Q-estim-0.5 - Locator = NA
    ## Column = 14 - Name = Pts_Recovery.Q-estim-0.75 - Locator = NA
    ## Column = 15 - Name = Pts_Recovery.Q-estim-1 - Locator = NA
    ## Column = 16 - Name = Pts_Recovery.Q-stdev-0 - Locator = NA
    ## Column = 17 - Name = Pts_Recovery.Q-stdev-0.5 - Locator = NA
    ## Column = 18 - Name = Pts_Recovery.Q-stdev-0.75 - Locator = NA
    ## Column = 19 - Name = Pts_Recovery.Q-stdev-1 - Locator = z1

    ## NULL

    # Graphics
    # Metal (the anamorphosis is used)
    p = plot(blocs, "Pts_Recovery.Q-estim-0",
             show.legend.raster=TRUE, legend.name.raster= "Q(0.0)")
    plot.decoration(p, title = "Point conditional expectation")


    # Tonnage (the anamorphosis is not used)
    p = plot(blocs, "Pts_Recovery.T-estim-0",
             show.legend.raster = TRUE, legend.name.raster = "T(0.0)")
    plot.decoration(p, title = "Point conditional expectation")


    plot(blocs, "Pts_Recovery.T-estim-1",
         show.legend.raster=TRUE, legend.name.raster = "T(1.0)")


    plot.decoration(p, title = "Point conditional expectation")


**Values computed by the anamorphosis are invalid.... Est-ce toujours
vrai (DR) ? Si non, supprimer le paragraphe QC suivant.**

    # QC of Conditional Expectation
    Ycuts <- anam$rawToGaussianVector(Zcuts)
    indsel = blocs["G_Pts.Gaussian.Z.stdev"] <= 0.

    # Tonnage à zc = 1.0
    blocs["T_1"] <- 1 - VectorHelper_pnormVec((Ycuts[4] - blocs["G_Pts.Gaussian.Z.estim"]) /
                                               blocs["G_Pts.Gaussian.Z.stdev"])
    blocs["T_1"][indsel] = NA
    p = plot.correlation(blocs, name1="Pts_Recovery.T-estim-1", name2="T_1", flagDiag=TRUE)
    plot.decoration(p, xlab = "Computed value", ylab = "Theoretical value",
                    title = "Conditional Expectation: T(1.0)")


    print(paste0(">>> Sum of differences for T(1.00): ",
                 sum(blocs["T_1"] != blocs["Pts_Recovery.T-estim-1"], 
                     na.rm=TRUE)))

    ## [1] ">>> Sum of differences for T(1.00): 0"

    # Tonnage à zc = 0.75

    blocs["T_0.75"] <- 1 - VectorHelper_pnormVec((Ycuts[3] -
                       blocs["G_Pts.Gaussian.Z.estim"]) / blocs["G_Pts.Gaussian.Z.stdev"])
    blocs["T_0.75"][indsel] = NA
    p = plot.correlation(blocs, name1="Pts_Recovery.T-estim-0.75", name2="T_0.75", 
                         flagDiag=TRUE)
    plot.decoration(p, xlab = "Computed value", ylab = "Theoretical value",
                    title = "Conditional Expectation: T(0.75)")


    print(paste0(">>> Sum of differences for T(0.75): ",
                 sum(blocs["T_0.75"] != blocs["Pts_Recovery.T-estim-0.75"],
                     na.rm=TRUE)))

    ## [1] ">>> Sum of differences for T(0.75): 0"

Disjunctive Kriging with the Gaussian model
-------------------------------------------

### Preprocessing

    #mm <- model$clone()
    # model of the Gaussian variable
    p = plot.varmod(vario, model)
    plot.decoration(p, title = "Variogram of the Gaussian variable")


    # attach the anamorphosis to the model
    err = model$setAnam(anam)

    # model of the raw variable as the anamorphosis is attached to the model
    model$setActiveFactor(-1)

    ## NULL

    p = plot.varmod(vario_raw, model)
    plot.decoration(p, title = "Variogram of the raw variable")


    model$setActiveFactor(0)

    ## NULL

Computing the point factors (i.e., the values of the Hermite
polynomials). The following plot gives the proportion of the variance
explained by the factors. It indicates the relevant number of factors
that should be considered afterwards.

    plot(1:(anam$getNFactor()-1),anam$cumulateVarianceRatio(1.), 
         ylim = c(0,1),
         xlab = "Number of factors", ylab = "Variance proportion",
         type = "b", pch = 19, col = "gray")


Plot the variogram model of the factors (rho^n in the Gaussian model)

    hmax = 0.2
    nf = anam$getNFactor()
    asCov = FALSE # Pb. if asCov = FALSE (should be 1 - rho^n and not (1-rho)^n)

    # model of the factors
    p = NULL
    for (i in 1:nf) 
    {
      model$setActiveFactor(i)
      p <- plot.model(model, hmax=hmax, asCov = asCov, padd=p)
    }
    p = p + geom_hline(yintercept = 1.0, colour = "red")
    plot.decoration(p, title = "Covariance of the different factors")


    # model of the raw variable
    model$setActiveFactor(-1)

    ## NULL

    p = plot.model(model, hmax=hmax, asCov = asCov)
    p = p + geom_hline(yintercept = varZ, colour = "red")
    plot.decoration(p, title = "Covariance of the raw variable")


Computation of the factors: H1, H2 et H3... They correspond to the
Hermite polynomials. Their presence will dictate the number of terms
used in the Hermite expansion in subsequent calculations.

    nfactor = 3
    err = anam$rawToFactor(db = data, nfactor = nfactor,
                           namconv = NamingConvention("F", FALSE))

### Simple Point Kriging over the blocs

What is performed? Compute the (simple) kriging of the factors. Simple
Kriging is used and point or block estimate can be used as kriging is a
linear operator.

    err = KrigingFactors(dbin   = data,   # Input data set (a Db structure)
                         dbout  = blocs,  # Output data set (a DbGrid structure)
                         model  = model, neigh, calcul= EKrigOpt_PONCTUAL(),
                         namconv = NamingConvention("DK_Pts"))
    blocs$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 28
    ## Maximum Number of UIDs       = 28
    ## Total number of samples      = 441
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.050     0.050
    ## Number :         21        21
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = G_Pts.Gaussian.Z.estim - Locator = NA
    ## Column = 3 - Name = G_Pts.Gaussian.Z.stdev - Locator = NA
    ## Column = 4 - Name = Pts_Recovery.T-estim-0 - Locator = NA
    ## Column = 5 - Name = Pts_Recovery.T-estim-0.5 - Locator = NA
    ## Column = 6 - Name = Pts_Recovery.T-estim-0.75 - Locator = NA
    ## Column = 7 - Name = Pts_Recovery.T-estim-1 - Locator = NA
    ## Column = 8 - Name = Pts_Recovery.T-stdev-0 - Locator = NA
    ## Column = 9 - Name = Pts_Recovery.T-stdev-0.5 - Locator = NA
    ## Column = 10 - Name = Pts_Recovery.T-stdev-0.75 - Locator = NA
    ## Column = 11 - Name = Pts_Recovery.T-stdev-1 - Locator = NA
    ## Column = 12 - Name = Pts_Recovery.Q-estim-0 - Locator = NA
    ## Column = 13 - Name = Pts_Recovery.Q-estim-0.5 - Locator = NA
    ## Column = 14 - Name = Pts_Recovery.Q-estim-0.75 - Locator = NA
    ## Column = 15 - Name = Pts_Recovery.Q-estim-1 - Locator = NA
    ## Column = 16 - Name = Pts_Recovery.Q-stdev-0 - Locator = NA
    ## Column = 17 - Name = Pts_Recovery.Q-stdev-0.5 - Locator = NA
    ## Column = 18 - Name = Pts_Recovery.Q-stdev-0.75 - Locator = NA
    ## Column = 19 - Name = Pts_Recovery.Q-stdev-1 - Locator = NA
    ## Column = 20 - Name = T_1 - Locator = NA
    ## Column = 21 - Name = T_0.75 - Locator = NA
    ## Column = 22 - Name = DK_Pts.F.1.estim - Locator = z1
    ## Column = 23 - Name = DK_Pts.F.2.estim - Locator = z2
    ## Column = 24 - Name = DK_Pts.F.3.estim - Locator = z3
    ## Column = 25 - Name = DK_Pts.F.1.stdev - Locator = NA
    ## Column = 26 - Name = DK_Pts.F.2.stdev - Locator = NA
    ## Column = 27 - Name = DK_Pts.F.3.stdev - Locator = NA

    ## NULL

Consistency between the simple kriging of the Gaussian variable and the
kriging of F.1. As F.1 = -Y, its kriging should follow (F.1)^{SK} =
(-Y)^{SK}

    # H_1 = -Y
    p = plot.correlation(blocs, 
                     name1 = "G_Pts.Gaussian.Z.estim",
                     name2 = "DK_Pts.F.1.estim", 
                     flagDiag = FALSE)
    p = p + geom_abline(intercept = 0, slope = -1, colour = "red")
    plot.decoration(p, xlab = "Gaussian variable",
                    ylab = "H_1 = -Y",
                    title = "Values of the Simple Kriging")


    p = plot.correlation(blocs, 
                     name1 = "G_Pts.Gaussian.Z.stdev",
                     name2 = "DK_Pts.F.1.stdev", 
                     flagDiag = TRUE)
    plot.decoration(p, xlab = "Gaussian variable",
                    ylab = "H_1 = -Y",
                    title = "Standard Deviation of the Simple Kriging error")


### Computation of the disjunctive kriging

The final aim is to estimate a *nonlinear* function of the Gaussian
variable
$$\\psi(Y(x)) = \\sum\_{n=0}^{nf} \\psi\_n \\, H\_n(Y(x))$$
 As kriging is a linear operator and the Hermite polynomials are
orthogonal,
$$(\\psi(Y(x))^{DK} =  \\sum\_{n=0}^{nf} \\psi\_n \\, H\_n^{SK}(Y(x))$$
 and
$$Var\\{\\psi(Y(x)) - \\psi(Y(x))^{DK}\\} =  \\sum\_{n=1}^{nf} \\psi\_n^2 \\, \\sigma^2\_{SK}(n)+ \\sum\_{n=nf+1}^{+\\infty} \\psi\_n^2= Var\\{\\psi(Y(x))\\}-\\sum\_{n=1}^{nf} \\psi\_n^2 \\, (1-\\sigma^2\_{SK}(n))$$

We have first to calculate the Hermite coefficients of the selectivity
curve or the nonlinear transform to be evaluated and then compute the
linear combinations.

-   Raw variable, *ψ*(*y*)=*ϕ*(*y*)

Hence, *ψ*<sub>*n*</sub> = *ϕ*<sub>*n*</sub>

-   Recovered ore, or the indicator above the cutoff,
    *ψ*(*y*)=1<sub>*y* ≥ *y*<sub>*c*</sub></sub>,

Hence, *ψ*<sub>0</sub> = 1 − *G*(*y*<sub>*c*</sub>) and
$\\psi\_n = -\\frac{g(y\_c)}{\\sqrt{n}} H\_{n-1}(y\_c)$ for *n* &gt; 0.

-   Recovered metal above the cutoff,
    *ψ*(*y*)=*ϕ*(*y*)×1<sub>*y* ≥ *y*<sub>*c*</sub></sub>,

Hence,
*ψ*<sub>*n*</sub> = ∑<sub>*p* ≥ 0</sub>*ϕ*<sub>*p*</sub> *α*<sub>*n*, *p*</sub>(*y*<sub>*c*</sub>)
with
*α*<sub>*n*, *p*</sub>(*y*<sub>*c*</sub>)=∫<sub>*y*<sub>*c*</sub></sub><sup>+∞</sup>*H*<sub>*n*</sub>(*u*)*H*<sub>*p*</sub>(*u*)*g*(*u*)*d**u*.

We should also be able to compute non linear transforms such as:

-   The value over a limit, *ψ*(*y*)=max(*y*<sub>*c*</sub>, *y*),

Hence,
*ψ*<sub>0</sub> = *g*(*y*<sub>*c*</sub>)+*y*<sub>*c*</sub> *G*(*y*<sub>*c*</sub>),
*ψ*<sub>1</sub> = *G*(*y*<sub>*c*</sub>)−1, and
$\\psi\_n = \\frac{g(y\_c)}{\\sqrt{n\\times(n-1)}} H\_{n-2}(y\_c)$ for
*n* &gt; 1.

-   The lognormal transform,
    $\\psi(y) = m \\times e^{\\sigma\\, y - \\frac{1}{2} \\sigma^2}$

Hence,
$\\psi\_n = m \\, e^{-\\frac{1}{2}\\sigma^2} \\, \\frac{(-\\sigma)^n}{\\sqrt{n!}}$
for *n* ≥ 0.

-   etc.

<!-- -->

    err = DisjunctiveKriging(db = blocs, 
                             anam = anam, 
                             selectivity = selectivity, 
                             name_est="DK_Pts.F.*.estim", 
                             name_std="DK_Pts.F.*.stdev",
                             namconv = NamingConvention("DK_Pts.Recovery", FALSE))

Simple Block Kriging over the blocs
-----------------------------------

Can be used to estimate the average of the point selectivity curve on a
block, e.g.
$$
(\\frac{1}{|v|}\\int\_v 1\_{Y(u) \\geq y\_c} du)^{DK} = \\sum\_{n=0}^{nf} \\psi\_n \\, \\{\\frac{1}{|v|}\\int\_v H\_n(Y(u)) du\\}^{SK}
$$
 The post-processing is identical using **DisjunctiveKriging**.

    # The number of estimated factors is defined by the number of selected variables 
    # in the input data base
    ndisc_B = c(5, 5)
    err = KrigingFactors(dbin   = data,   # Input data set (a Db structure)
                         dbout  = blocs,  # Output data set (a DbGrid structure)
                         model  = model, neigh, calcul= EKrigOpt_BLOCK(), ndisc = ndisc_B,
                         namconv = NamingConvention("DK_Blk"))
    blocs$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 34
    ## Maximum Number of UIDs       = 34
    ## Total number of samples      = 441
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.050     0.050
    ## Number :         21        21
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = G_Pts.Gaussian.Z.estim - Locator = NA
    ## Column = 3 - Name = G_Pts.Gaussian.Z.stdev - Locator = NA
    ## Column = 4 - Name = Pts_Recovery.T-estim-0 - Locator = NA
    ## Column = 5 - Name = Pts_Recovery.T-estim-0.5 - Locator = NA
    ## Column = 6 - Name = Pts_Recovery.T-estim-0.75 - Locator = NA
    ## Column = 7 - Name = Pts_Recovery.T-estim-1 - Locator = NA
    ## Column = 8 - Name = Pts_Recovery.T-stdev-0 - Locator = NA
    ## Column = 9 - Name = Pts_Recovery.T-stdev-0.5 - Locator = NA
    ## Column = 10 - Name = Pts_Recovery.T-stdev-0.75 - Locator = NA
    ## Column = 11 - Name = Pts_Recovery.T-stdev-1 - Locator = NA
    ## Column = 12 - Name = Pts_Recovery.Q-estim-0 - Locator = NA
    ## Column = 13 - Name = Pts_Recovery.Q-estim-0.5 - Locator = NA
    ## Column = 14 - Name = Pts_Recovery.Q-estim-0.75 - Locator = NA
    ## Column = 15 - Name = Pts_Recovery.Q-estim-1 - Locator = NA
    ## Column = 16 - Name = Pts_Recovery.Q-stdev-0 - Locator = NA
    ## Column = 17 - Name = Pts_Recovery.Q-stdev-0.5 - Locator = NA
    ## Column = 18 - Name = Pts_Recovery.Q-stdev-0.75 - Locator = NA
    ## Column = 19 - Name = Pts_Recovery.Q-stdev-1 - Locator = NA
    ## Column = 20 - Name = T_1 - Locator = NA
    ## Column = 21 - Name = T_0.75 - Locator = NA
    ## Column = 22 - Name = DK_Pts.F.1.estim - Locator = NA
    ## Column = 23 - Name = DK_Pts.F.2.estim - Locator = NA
    ## Column = 24 - Name = DK_Pts.F.3.estim - Locator = NA
    ## Column = 25 - Name = DK_Pts.F.1.stdev - Locator = NA
    ## Column = 26 - Name = DK_Pts.F.2.stdev - Locator = NA
    ## Column = 27 - Name = DK_Pts.F.3.stdev - Locator = NA
    ## Column = 28 - Name = DK_Blk.F.1.estim - Locator = z1
    ## Column = 29 - Name = DK_Blk.F.2.estim - Locator = z2
    ## Column = 30 - Name = DK_Blk.F.3.estim - Locator = z3
    ## Column = 31 - Name = DK_Blk.F.1.stdev - Locator = NA
    ## Column = 32 - Name = DK_Blk.F.2.stdev - Locator = NA
    ## Column = 33 - Name = DK_Blk.F.3.stdev - Locator = NA

    ## NULL

Comparing Point and Block estimations and standard deviation of
Estimation errors (this comparison is performed on the results of the
first factors only).

    p = plot.correlation(blocs, name1 = "DK_Pts*1.estim", name2 = "DK_Blk*1.estim",
                         flagDiag=TRUE)
    plot.decoration(p, xlab = "Point estimation", ylab = "Block estimation", 
                    title = "First factor (Gaussian model)")


    p = plot.correlation(blocs, name1 = "DK_Pts*2.estim", name2 = "DK_Blk*2.estim",
                         flagDiag=TRUE)
    plot.decoration(p, xlab = "Point estimation", ylab = "Block estimation", 
                     title = "Second factor (Gaussian model)")


    p = plot.correlation(blocs, name1 = "DK_Pts*3.estim", name2 = "DK_Blk*3.estim",
                     flagDiag = TRUE)
    plot.decoration(p, xlab = "Point estimation", ylab = "Block estimation", 
                     title = "Third factor (Gaussian model)")


Simple Block Kriging over the panels
------------------------------------

Can be used to estimate the average of the point selectivity curve on a
panel.

XF: factors should be used instead if DK has to be computed.

    # The number of estimated factors is defined by the number of selected variables 
    # in the input data base
    ndisc_P = c(10, 10)
    data$setLocator("Gaussian.Z", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    err = kriging(dbin   = data,   # Input data set (a Db structure).
                  dbout  = panels, # Output data set (a DbGrid structure)
                  model  = model, neigh, calcul= EKrigOpt_BLOCK(), ndisc = ndisc_P,
                  namconv = NamingConvention("DK_Blk"))
    panels$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 4
    ## Maximum Number of UIDs       = 4
    ## Total number of samples      = 25
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      0.250     0.250
    ## Number :          5         5
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = DK_Blk.Gaussian.Z.estim - Locator = z1
    ## Column = 3 - Name = DK_Blk.Gaussian.Z.stdev - Locator = NA

    ## NULL

    p = plot(panels, "DK_Blk.*.estim", show.legend.raster = TRUE) 
    plot.decoration(p, xlab = "Easting", ylab = "Northing", 
                title = "Estimate of the first factor (panels)")


    p = plot(panels, "DK_Blk.*.stdev", show.legend.raster = TRUE)
    plot.decoration(p, xlab = "Easting", ylab = "Northing", 
                  title = "Std. of the first factor estimate (panels)")


Cross validation
----------------

Estimators are tested on the fine grid: the reference simulated values
are estimated using the subset of 500 random points. Cross plot of the
actual value versus the estimated value are plotted and scoring rules
evaluated

-   Mean error $ME = \\frac{1}{N} \\sum\_{i = 1}^N (z\_i - z\_i^\*))$
-   Mean absolute error
    $MAE = \\frac{1}{N} \\sum\_{i = 1}^N |z\_i - z\_i^\*|)$
-   Root mean square error
    $RMSE = \\sqrt{\\frac{1}{N} \\sum\_{i = 1}^N (z\_i - z\_i^\*)^2}$

-   Mean standardized error
    $MSE = \\frac{1}{N} \\sum\_{i = 1}^N \\frac{z\_i - z\_i^\*}{\\sigma\_i^\*}$
-   Mean absolute standardized error
    $MASE = \\frac{1}{N} \\sum\_{i = 1}^N |\\frac{z\_i - z\_i^\*}{\\sigma\_i^\*}|$
-   Root mean square standardized error
    $RMSSE = \\sqrt{\\frac{1}{N} \\sum\_{i = 1}^N (\\frac{z\_i - z\_i^\*}{\\sigma\_i^\*})^2}$

TODO: - results with the disjunctive kriging are dubious (see xplot for
Z and Ind) - k-fold cross validation to be implemented ?

### Estimation on the initial fine grid

    # cleaning the target data base
    err = grid$deleteColumns("SK.*")
    err = grid$deleteColumns("CE.*")
    err = grid$deleteColumns("DK.*")
    err = grid$deleteColumns("sel-.*")
    # Selectivity curve : Z
    zcuts <- round(quantile(data$getColumn("Z"), probs = c(0.5, 0.75)),3)
    ycuts <- round(anam$rawToGaussianVector(zcuts),3)

    Z_est = Selectivity_createByKeys(keys = c("Z"), zcuts=c(0.0),
                                     flag_est=TRUE, flag_std=TRUE)

    # Global experimental selectivity, calculated form the experimental Data Set
    grid$setLocator("Z", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    knitr::kable(Z_est$eval(grid)$toTL(), digits=3,
                 caption = "Selectivity curves computed on fine grid")

<table>
<caption>Selectivity curves computed on fine grid</caption>
<thead>
<tr class="header">
<th align="right">Z-Cut</th>
<th align="right">T-estim</th>
<th align="right">Q-estim</th>
<th align="right">B-estim</th>
<th align="right">M-estim</th>
<th align="right">T-stdev</th>
<th align="right">Q-stdev</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1.561</td>
<td align="right">1.561</td>
<td align="right">1.561</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
</tbody>
</table>

    # Simple kriging of Z
    err = model_raw$delDrift(0)
    err = model_raw$setMean(0, mean(data["Z"]))
    data$setLocator("Z", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    err = kriging(data, grid, model_raw, neigh, 
                  calcul = EKrigOpt_PONCTUAL(),
                  namconv= NamingConvention("SK"))

    # Conditional expectation
    err = model_Y$delDrift(0)
    err = model_Y$setMean(0, 0.)
    data$setLocator("Gaussian.Z", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    err = kriging(data, grid, model_Y, neigh, 
                  calcul = EKrigOpt_PONCTUAL(),
                  namconv= NamingConvention("SK"))
    err = ConditionalExpectation(db = grid, 
                                 anam = anam, 
                                 selectivity = Z_est, 
                                 name_est="SK.Gaussian.Z.estim", 
                                 name_std="SK.Gaussian.Z.stdev",
                                 namconv = NamingConvention("CE", FALSE))
    grid$setName("CE.Z-estim", "CE.Z.estim")

    ## NULL

    grid$setName("CE.Z-stdev", "CE.Z.stdev")

    ## NULL

    # Disjunctive kriging
    data$setLocator("F.*", ELoc_Z(), cleanSameLocator=TRUE)

    ## NULL

    err = KrigingFactors(dbin   = data,   # Input data set (a Db structure)
                         dbout  = grid,  # Ouput data set (a DbGrid structure)
                         model  = model, neigh, calcul= EKrigOpt_PONCTUAL(),
                         namconv = NamingConvention("SK"))

    # Disjunctive kriging for Z
    err = DisjunctiveKriging(db = grid, 
                             anam = anam, 
                             selectivity = Z_est , 
                             name_est="SK.F.*.estim", 
                             name_std="SK.F.*.stdev",
                             namconv = NamingConvention("DK", FALSE))
    err = grid$setName("DK.Z-estim", "DK.Z.estim")
    err = grid$setName("DK.Z-stdev", "DK.Z.stdev")

    # Disjunctive kriging for T and Q
    Z_est = Selectivity_createByKeys(keys = c("T", "Q"), zcuts=zcuts,
                                      flag_est=TRUE, flag_std=TRUE)
    err = DisjunctiveKriging(db = grid, 
                             anam = anam, 
                             selectivity = Z_est , 
                             name_est="SK.F.*.estim", 
                             name_std="SK.F.*.stdev",
                             namconv = NamingConvention("DK", FALSE))

### QC of the Disjunctive Kriging functions

    grid$deleteColumns("DK.*.qc.*")

    ## NULL

    N  = 1000
    Y = seq(from = -3.5, to = 3.5, length.out = 1000)
    H = matrix(NaN, nrow = length(Y), ncol = N)
    for (i in seq_along(Y)) {
      H[i,] = hermitePolynomials(Y[i], 1.0, N)
    }

    # Lognormal
    m = 1.0 ; sigma = 1.5
    psi.LN = hermiteLognormal(mean = m, sigma = sigma, nbpoly = N)
    plot(Y, H %*% psi.LN, type = "l", col = "red", 
         xlab = "Gaussian value", ylab = "Raw value", 
         main = paste0("Lognormal m=", m, " and sigma = ", sigma))
    lines(Y, m * exp(sigma*Y - 1/2 * sigma^2), col = "blue", lwd = 2, lty = 2)
    legend("topleft", legend = c("Hermite", "Theoretical"), 
           col = c("red", "blue"), lty = c(1,2))


    # QC of the indicator decomposition (OK)
    # psi.LN.qc = rep(NaN, N)
    # psi.LN.qc = m * (-sigma)^(0:(N-1)) / sqrt(factorial(0:(N-1)))
    # plot(psi.LN, psi.LN.qc)
    # abline(a = 0, b = 1, col = "red")

    # Indicator
    yc = ycuts[1]
    psi.ind = hermiteCoefIndicator(yc = yc, nbpoly = N)
    plot(Y, H %*% psi.ind, type = "l", col = "red", 
         xlab = "Gaussian value", ylab = "Raw value", 
         main = paste0("Indicator of ", yc))
    lines (Y, as.numeric(Y >= yc), col = "blue", lty = 2)
    abline(h = c(0,1), col = "gray", lty = 2)
    legend("topleft", legend = c("Hermite", "Theoretical"), 
           col = c("red", "blue"), lty = c(1,2))


    # QC of the indicator decomposition (OK)
    # psi.ind.qc = rep(NaN, N)
    # psi.ind.qc[1] = 1 - pnorm(yc)
    # psi.ind.qc[-1] = -dnorm(yc)/sqrt(1:(N-1))*hermitePolynomials(yc, 1.0, N-1)
    # plot(psi.ind, psi.ind.qc)
    # abline(a = 0, b = 1, col = "red")

    # Metal
    yc = 2.0
    psi.metal = hermiteCoefMetal(yc = yc, phi = psi.LN)
    plot(Y, H %*% psi.metal, type = "l", col = "red", 
         xlab = "Gaussian value", ylab = "Raw value", 
         main = paste0("Metal for cutoff = ", yc))
    lines (Y, m * exp(sigma * Y - 1/2 * sigma^2) * as.numeric(Y >= yc), 
           col = "blue", lty = 2)
    legend("topleft", legend = c("Hermite", "Theoretical"), 
           col = c("red", "blue"), lty = c(1,2))


    # QC of the metal decomposition (OK)
    # psi.metal.qc = matrix(hermiteIncompleteIntegral(yc, N)$getValues(), N, N) %*% psi.LN
    # plot(psi.metal, psi.metal.qc)
    # abline(a = 0, b = 1, col = "red")


    # generic function to test DK
    compute_DK <- function(psi, name) {
    N   <- length(psi)
    stopifnot(nfactor < N)

    nm_estim    <- paste0(name, ".estim")
    nm_stdev    <- paste0(name, ".stdev")

    nm_qc_estim <- paste0(name, ".qc.estim")
    nm_qc_stdev <- paste0(name, ".qc.stdev")

    grid[nm_qc_estim] <- psi[1] + Hn.est %*% psi[-1][1:nfactor]
    grid[nm_qc_stdev] <- sqrt(Hn.std^2 %*% (psi[-1][1:nfactor]^2)  + sum(psi[-1][(nfactor+1):(N-1)]^2))

    # Raw variable
    p = plot.correlation(grid, name1 = nm_estim, name2 = nm_qc_estim, flagDiag=TRUE)
    plot.decoration(p, xlab = nm_estim, ylab = nm_qc_estim,
                    title = "Disjunctive Kriging - Estimation value")
    p = plot.correlation(grid, name1 = nm_stdev, name2 = nm_qc_stdev, flagDiag=TRUE)
    plot.decoration(p, xlab = nm_stdev, ylab = nm_qc_stdev,
                    title = "Disjunctive Kriging - Standard deviation")
    NULL
    }


    # QC of the DK
    Hn.est <- matrix(grid$getColumns(names = c("SK.F.*.estim"), useSel = FALSE),
                     ncol = nfactor, nrow = grid$getSampleNumber())
    Hn.std <- matrix(grid$getColumns(names = c("SK.F.*.stdev"), useSel = FALSE),
                    ncol = nfactor, nrow = grid$getSampleNumber())

    # estimation of Z using DK
    psi <- anam$getPsiHns()
    compute_DK(psi = psi, name = "DK.Z")

    ## NULL

    # estimation of the tonnage Z > yc
    psi <- hermiteCoefIndicator(ycuts[1], anam$getNFactor())
    print(paste0("Std(I) = ", round(sqrt(psi[1]*(1 - psi[1])),2)))

    ## [1] "Std(I) = 0.5"

    print(paste0("Std(I)*= ", round(sqrt(sum(psi[-1]^2)),2)))

    ## [1] "Std(I)*= 0.47"

    grid$setName(paste0("DK.T-estim-", zcuts[1]), "DK.T1.estim")

    ## NULL

    grid$setName(paste0("DK.T-stdev-", zcuts[1]), "DK.T1.stdev")

    ## NULL

    compute_DK(psi = psi, name = "DK.T1")

    ## NULL

    # estimation of the metal Z > yc
    psi <- hermiteCoefMetal(yc = ycuts[1], phi = anam$getPsiHns())
    grid$setName(paste0("DK.Q-estim-", zcuts[1]), "DK.Q1.estim")

    ## NULL

    grid$setName(paste0("DK.Q-stdev-", zcuts[1]), "DK.Q1.stdev")

    ## NULL

    compute_DK(psi = psi, name = "DK.Q1")

    ## NULL

Scores of the estimators
------------------------

    # Evaluation
    sel <- (grid["SK.Z.stdev"] > 0.2)
    grid["sel"] <- sel
    grid$setLocator("sel", ELoc_SEL(), cleanSameLocator=TRUE)

    ## NULL

    nm <-c("Z",
           paste0(c("SK", "DK", "CE"), ".Z.estim"),
           paste0(c("SK", "DK", "CE"), ".Z.stdev")
    )

    # Correlation actual values vs. estimated values
    for (est in c("SK", "DK", "CE")) {
      p = plot.correlation(grid, name1=paste0(est, ".Z.estim"), name2="Z", flagDiag = TRUE)
      plot.decoration(p, xlab = "Estimated value", ylab = "Actual value",
                      title = paste(est, " of raw variable"))
      }

    # Comparison of the estimators
    p = plot.correlation(grid, name1="CE.Z.estim", name2="SK.Z.estim", flagDiag = TRUE)
    plot.decoration(p, xlab = "Conditional expectation", ylab = "Simple kriging",
                    title = "Point raw variable")


    p = plot.correlation(grid, name1="DK.Z.estim", name2="SK.Z.estim", flagDiag = TRUE)
    plot.decoration(p, xlab = "Disjunctive Kriging", ylab = "Simple kriging",
                    title = "Point raw variable")


    # Maps of the estimators
    for (est in c("SK", "DK", "CE")){
      for (type in c("estim", "stdev")){

      p = plot(grid, paste0(est, ".Z.", type),
               show.legend.raster=TRUE, legend.name.raster=type) 
      plot.decoration(p, xlab = "Easting", ylab = "Northing", 
                      title = paste("Estimator = ", est)) 
      
      if (type == "estim") col_fill = "orange"
      if (type == "stdev") col_fill = "skyblue"
      
      p = plot.hist(grid, name = paste0(est, ".Z.", type), nbins = 50, 
                col = "gray", fill = col_fill)
      plot.decoration(p, xlab = paste0(est, ".Z.", type),
                title = paste("Estimator = ", est))
    }
    }

    # Mono variate statistics of the estimators
    knitr::kable(dbStatisticsMonoT(grid, nm, opers=opers)$toTL(), digits=3,
        caption = "Statistics of the estimated values")

<table>
<caption>Statistics of the estimated values</caption>
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
<td align="left">Z</td>
<td align="right">9500</td>
<td align="right">0.154</td>
<td align="right">8.198</td>
<td align="right">1.562</td>
<td align="right">0.846</td>
</tr>
<tr class="even">
<td align="left">SK.Z.estim</td>
<td align="right">9500</td>
<td align="right">0.523</td>
<td align="right">4.924</td>
<td align="right">1.553</td>
<td align="right">0.496</td>
</tr>
<tr class="odd">
<td align="left">DK.Z.estim</td>
<td align="right">9500</td>
<td align="right">0.430</td>
<td align="right">4.087</td>
<td align="right">1.555</td>
<td align="right">0.448</td>
</tr>
<tr class="even">
<td align="left">CE.Z.estim</td>
<td align="right">9500</td>
<td align="right">0.518</td>
<td align="right">4.305</td>
<td align="right">1.554</td>
<td align="right">0.455</td>
</tr>
<tr class="odd">
<td align="left">SK.Z.stdev</td>
<td align="right">9500</td>
<td align="right">0.473</td>
<td align="right">0.834</td>
<td align="right">0.668</td>
<td align="right">0.071</td>
</tr>
<tr class="even">
<td align="left">DK.Z.stdev</td>
<td align="right">9500</td>
<td align="right">0.552</td>
<td align="right">0.832</td>
<td align="right">0.701</td>
<td align="right">0.052</td>
</tr>
<tr class="odd">
<td align="left">CE.Z.stdev</td>
<td align="right">9500</td>
<td align="right">0.210</td>
<td align="right">1.375</td>
<td align="right">0.681</td>
<td align="right">0.194</td>
</tr>
</tbody>
</table>

    # Scoring rules
    scoringRules <- function(db, nm_val, nm_est, nm_std){
      val  <- matrix(
        db$getColumns(names = c(nm_val, nm_est, nm_std), useSel = TRUE),
        nrow = db$getActiveSampleNumber(),
        ncol = 3, byrow = FALSE)
      err  <- (val[,1] - val[,2])
      nerr <- err / val[,3]
      c(
        mean(err), mean(abs(err)),  sqrt(mean(err^2)),
        mean(nerr),mean(abs(nerr)), sqrt(mean(nerr^2))
      )
    }

    res <- matrix(c(
      scoringRules(grid, "Z", "SK.Z.estim", "SK.Z.stdev"),
      scoringRules(grid, "Z", "DK.Z.estim", "DK.Z.stdev"),
      scoringRules(grid, "Z", "CE.Z.estim", "CE.Z.stdev")),
      nrow = 3, ncol = 6, byrow = TRUE)
    colnames(res) <- c("ME", "MAE", "RMSE", "MSE", "MASE", "RMSSE")
    rownames(res) <- c("SK", "DK", "CE")
    knitr::kable(res, digits=4,
                 caption = "Scores of point estimators of Z using 500 data points")

<table>
<caption>Scores of point estimators of Z using 500 data points</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">ME</th>
<th align="right">MAE</th>
<th align="right">RMSE</th>
<th align="right">MSE</th>
<th align="right">MASE</th>
<th align="right">RMSSE</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">SK</td>
<td align="right">0.0093</td>
<td align="right">0.4839</td>
<td align="right">0.6721</td>
<td align="right">0.0168</td>
<td align="right">0.7225</td>
<td align="right">0.9996</td>
</tr>
<tr class="even">
<td align="left">DK</td>
<td align="right">0.0069</td>
<td align="right">0.4831</td>
<td align="right">0.6732</td>
<td align="right">0.0121</td>
<td align="right">0.6865</td>
<td align="right">0.9540</td>
</tr>
<tr class="odd">
<td align="left">CE</td>
<td align="right">0.0085</td>
<td align="right">0.4807</td>
<td align="right">0.6694</td>
<td align="right">-0.0076</td>
<td align="right">0.6969</td>
<td align="right">0.9108</td>
</tr>
</tbody>
</table>

### k-fold cross validation of the data set

TODO: 2023-01-31, it should be completed (XF)

    K <- 10 # Number of folds
    res_xval <- data$clone() # duplicate the data base
    # cleaning the target data base
    res_xval$deleteColumns("SK.*")
    res_xval$deleteColumns("CE.*")
    res_xval$deleteColumns("DK.*")
    res_xval$deleteColumns("code")

    # definition of the folds
    res_xval["code"] <- sample(x = 1:K, size = res_xval$getSampleNumber(), replace = TRUE)
    res_xval$setLocator("code", ELoc_C(), cleanSameLocator=TRUE)

    # generic function for k-fold cross validation
    # the db should have:
    # - db: the data base with a "code" variable defining the folds
    # - fn_estim: a R function which defines the interpolation process. Its prototype is int f(dbin, dbout)
    # - nameconv: the naming convention 
    # the outputs are
    # - err: code of error
    # - db: computed variables (estimated value and the estimation Std.)

    kfold_compute <- function(db, fn_estim, namconv = NamingConvention("kfold")){
      # TODO
      err = 0
      err
    }
    # the evaluation of the estimator is: 
    # 1) define the fn_estim
    # 2) compute the cross validation using the *kfold_compute* function
    # 3) compute the scoring using the *scoringRules* function
    # 4) display the results *knitr::kable(...))*

    # TODO: not completed

Discussion and further works
============================

Still to be done:

-   mise en oeuvre de la validation croisée et comparaison avec les
    ressources récupérables calculée sur la simulation initiale.

-   Exploiter les résultats du KD, EC et UC pour définir des ressources
    récupérables sur un polygone

References
==========

-   Rivoirard, J. (1994). Introduction to disjunctive kriging and
    non-linear geostatistics. Number 551.021 R626i. Clarendon Press.

