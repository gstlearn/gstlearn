General Introduction
====================

The objective of this document is to provide some guides for
manipulation of Hermite polynomials. In praticular, it will give some
hints for calculating of complex expressions (e.g. metal bove cutoff)
and be double checked using a Monte Carlo approach.

It will also demonstrate the quality of the polynomial expansion for a
given known distribution of the data (i.e. lognormal). Finally it will
serve in demonstrating the ability of representing any function through
a Hermite polynomial expansion.

Part I
======

The objective of the first part is to compute

*ψ*<sub>1</sub>(*α*, *β*)=∫*ϕ*(*α* + *β* *u*) *g*(*u*) *d**u*

and

*ψ*<sub>2</sub>(*α*, *β*)=∫*ϕ*<sup>2</sup>(*α* + *β* *u*) *g*(*u*) *d**u*

for a function *ϕ* in order to compute the conditional expectation and
the conditional variance, usually defined by its expansion in terms of
Hermite polynomials:

$$
\\phi(y) = \\sum\_{n=0}^{N} a\_n \\eta\_n(y)
$$

In particular, we may be interested in the case where *α* = *r**y* and
$\\beta = \\sqrt { { 1-r^2 } }$. Note that |*β*|&lt;1.

Many methods can be used to compute these draw.matrix:

-   the Monte-Carlo integration (easy)

-   the computation of the Hermite coefficient for *ϕ* and
    *ϕ*<sup>2</sup>

$$
\\int\_{-\\infty}^{+\\infty} \\eta\_n(\\alpha + \\beta \\, u)\\, g(u) \\, du = (1-\\beta^2)^n \\, \\eta\_n(\\frac{\\alpha}{\\sqrt{1-\\beta^2}})
$$

Evaluation on the lognormal case
================================

Processing parameters

    nbpoly = 100
    nbsimu = 1000
    params = list(nbpoly=nbpoly, nbsimu=nbsimu, yc=1.5)

We consider the Lognormal case (where the anamorphosis is known
exactly).

    m     <- 2.5
    sigma <- 0.5
    hn     = hermiteLognormal(m, sigma, params$nbpoly)

We calculate the various functions (expectation, ore and metal quantity,
for expectation and standard deviation) for different values of *y*
varying from -3 to 3 and different values of *r* varying from 0 to 1.
The results are stored in a Db (organized as a matrix) whose parameters
are given next.

    y0  = -3
    dy  = 0.1
    dr  = 0.1
    y   = seq(-3,3,by=dy)
    ny  = length(y)
    r   = seq(1,0,by=-dr)
    nr  = length(r)

Comparisons
===========

We evaluate the variable, the ore and metal quantity above a cutoff, in
conditional expectation or in conditional standard deviation.

Calculating all elements using either Hermite expansion or Monte Carlo
Simulations

    draw.matrix(type=1,calcul=1,flag.est=TRUE,params,flag.cut=FALSE)
    draw.matrix(type=2,calcul=1,flag.est=TRUE,params,flag.cut=FALSE)
    draw.correlation(calcul=1,flag.est=TRUE,params)

    draw.matrix(type=1,calcul=2,flag.est=TRUE,params)
    draw.matrix(type=2,calcul=2,flag.est=TRUE,params)
    draw.correlation(calcul=2,flag.est=TRUE,params)

    draw.matrix(type=1,calcul=3,flag.est=TRUE,params)
    draw.matrix(type=2,calcul=3,flag.est=TRUE,params)
    draw.correlation(calcul=3,flag.est=TRUE,params)

    draw.matrix(type=1,calcul=1,flag.est=FALSE,params,flag.cut=FALSE)
    draw.matrix(type=2,calcul=1,flag.est=FALSE,params,flag.cut=FALSE)
    draw.correlation(calcul=1,flag.est=FALSE,params)

    draw.matrix(type=1,calcul=2,flag.est=FALSE,params)
    draw.matrix(type=2,calcul=2,flag.est=FALSE,params)
    draw.correlation(calcul=2,flag.est=FALSE,params)

    draw.matrix(type=1,calcul=3,flag.est=FALSE,params)
    draw.matrix(type=2,calcul=3,flag.est=FALSE,params)
    draw.correlation(calcul=3,flag.est=FALSE,params)

Part II
=======

We use the lognormal model to test some computations done with Hermite
polynomials.

We consider the second order stationary model
$Z\_{\\lambda}(x) = e^{\\lambda Y(x) - \\frac{1}{2} \\lambda^2}$ where
*Y*(*x*) is a centered Gaussian model with autocorrelation *ρ*:

-   *E*{*Y*}=0

-   *C**o**v*(*Y*(*x*),*Y*(*x* + *h*)) = *E*{*Y*(*x*)*Y*(*x* + *h*)} = *ρ*(*h*)

-   *E*{*Z*<sub>*λ*</sub>}=1

-   *C**o**v*(*Z*<sub>*λ*</sub>(*x*),*Z*<sub>*λ*</sub>(*x* + *h*)) = *E*{*Z*<sub>*λ*</sub>(*x*)*Z*<sub>*λ*</sub>(*x* + *h*)} − 1 = *e*<sup>*λ*<sup>2</sup>*ρ*(*h*)</sup> − 1

The anamorphosis *ϕ*<sub>*λ*</sub> maps the Gaussian field into the
lognormal SOS model
$Z\_{\\lambda} = \\phi\_{\\lambda}(Y) = \\sum\_{n=0}^{+\\infty}\\phi\_n{(\\lambda)}H\_n(Y)$.

The Hermite coefficients are
$\\phi\_n{(\\lambda)} = \\frac{(-\\lambda)^n}{\\sqrt{n!}}$ and we have
$$
C\_{\\lambda}(h) = e^{\\lambda^2 \\rho } - 1 = \\sum\_{n = 1}^{+\\infty} \\phi\_n^2(\\lambda) \\rho^{n}.
$$

The recursion formula to compute the values of the Hermite polynomials
are: *H*<sub>0</sub>(*y*)=1, *H*<sub>1</sub>(*y*)= − *y*, and
$$
H\_{n+1}(y) = -\\frac{1}{\\sqrt{n+1}} y H\_n(y) - \\sqrt{\\frac{n}{n+1}} H\_{n-1}(y)
$$

    OptCst_defineByKey("ASP",1)

    ## NULL

The next lines define some functions which can be explitely written in
the lognormal case

    GM_hermite <- function(y, nh){
      h <- matrix(NaN, nrow = length(y), ncol = 1+nh)
      h[,1] <- 1.0
      h[,2] <- -y
      for (i in 2:nh)
      {
        h[,1+i] <- -y*h[,i]/sqrt(i) - sqrt((i-1)/i)*h[,i-1]
      }
      h
    }

    # fonct 
    GM_LN_psi <- function(lambda, nh){
      (-lambda)^{0:nh} / sqrt(factorial(0:nh))
    }
    GM_LN_Cov_theo <- function(lambda){
      function(rho) {exp(lambda^2 * rho) - 1}
    }
    GM_LN_Cov_nh <- function(lambda, nh)
    {
      psi <- GM_LN_psi(lambda = lambda, nh = nh)
      foo <- function(rho) {outer(X = rho, Y = 1:nh, FUN = function(r,n){r^n}) %*% psi[-1]^2}
    }

The next paragraph initiates a data set and runs the previous functions.
The output will serve as reference for comparison with rgstlearn
functions.

    np <- 10000
    nh <- 100
    lambda <- 1.0
    # anamorphosis in the lognormal case (theoretical)
    yy <- qnorm(p = (1:np)/(np+1), mean = 0, sd = 1.0)
    aymin = range(yy)[1]
    aymax = range(yy)[2]
    zz_theo <- exp(lambda*yy - lambda^2/2)

Perfoming the tests with the functions programmed internally
------------------------------------------------------------

Anamorphosis in the lognormal case (Hermite pol.)

    zz_hm   <- GM_hermite(y = yy, nh = nh) %*% GM_LN_psi(lambda = lambda, nh = nh)

    p = plot.XY(yy, zz_theo, join=TRUE, color = "red")
    p = plot.XY(yy, zz_hm, color = "blue", linetype="dashed", padd=p)
    p = plot.decoration(p, xlab = "Y", ylab = "Z", title = "Lognormal transform")
    p


    p = plot.XY(zz_theo, zz_hm, join=TRUE, color = "blue")
    p = plot.decoration(p , xlab = "theoretical", ylab = "Hermite pol.", 
                        title = "Lognormal transform")
    p


Covariance

    r <- seq(from = 0, to =1, by = 0.01)
    C_theo <- GM_LN_Cov_theo(lambda = lambda)
    C_nh   <- GM_LN_Cov_nh(lambda = lambda, nh = nh)
    p = plot.XY(r, C_theo(r), color = "red", join=TRUE)
    p = plot.XY(r,C_nh(r), join=TRUE, color = "blue", linetype="dashed",  padd=p)
    p = plot.decoration(p, xlab = "Covariance of Y", ylab = "Covarariance of Z", 
                        title = "Lognormal transform")
    p


Performing similar tests with gstlearn
--------------------------------------

    db <- Db_createFromSamples(np,ELoadBy_SAMPLE(),zz_theo)
    db$setNameByUID(db$getLastUID(),"zz_theo")

    ## NULL

    ana = AnamHermite_create(nh)
    err = ana$fit(db, name="zz_theo")

    # Patch the Hermite coefficients to have exact numerical values
    localPsi <- GM_LN_psi(lambda = lambda, nh = nh-2)
    err = ana$setPsiHns(localPsi)

    p = plot.anam(ana, aymin=aymin, aymax=aymax)
    p = plot.XY(yy, zz_theo, color="red", linetype="dashed", padd=p)
    p = plot.decoration(p, title="Lognormal Transform")
    p


    #legend("topleft", legend = c("theoretical", "RGeostats"), col = c("red", "black"), lty = c(1,1))

    zz_rg = ana$gaussianToRawVector(yy)
    p = plot.XY(zz_theo, zz_rg, join=TRUE, color="blue")
    p = plot.decoration(p, xlab="theoretical", ylab="Hermite pol.", title="Lognormal Transform")
    p
    #legend("topleft", legend = c("transform", "Y=X"), col = c("blue", "red"), lty = c(1,2))

Covariance

    mod_y <- Model(nvar=1, ndim=2)
    mod_y$addCovFromParam(ECov_SPHERICAL(),range=1,sill=1)

    ## NULL

    h     <- seq(from = 0, to =1, by = 0.01)
    r = mod_y$evalIvarNpas(0,0,h)
    C_theo<- GM_LN_Cov_theo(lambda = lambda)
    cov_theo<- C_theo(rho = r) 

    # Creation of the lognormal model (il faudrait pouvoir créer un modèle gaussien à partir de la covariance de la gaussienne mod_Y et de l'anamorphose gaussienne ana)
    # mod_Z <- model.create.GM(model = mod_Y, anam = ana)
    # cov_rg<- model.eval(model = mod_Z, h = h)

    model.eval.GM <- function(model, anam, h){
      r = model$evalIvarNpas(0,0,h)
      psi <- anam$getPsiHns()
      nh <- length(psi) - 1
      outer(X = r, Y = 1:nh, FUN = function(rho,n){rho^n}) %*% (psi[-1])^2
    }
    cov_rg <- model.eval.GM(model = mod_y, anam = ana, h = h)

    p = plot.XY(r, cov_theo, color = "red", join=TRUE) 
    p = plot.decoration(p, xlab = "Covariance of Y", ylab = "Covariance of Z", 
                        title = "Lognormal transform")
    p = plot.XY(r, cov_rg, color = "blue", linetype = "dashed", padd=p)
    p


    #legend("topleft", legend = c("theoretical", "RGeostats."), col = c("red", "blue"), lty = c(1,2))

The empirical anamorphosis does not reproduce the lognormal transform
and the Hermite coefficient need to be altered.

Part III
========

The observed variable is

*Y*<sub>*y*<sub>*c*</sub></sub> = *ϕ*<sub>*y*<sub>*c*</sub></sub>(*Y*)=*Y* × 1<sub>*Y* ≥ *y*<sub>*c*</sub></sub> + *y*<sub>*c*</sub> × 1<sub>*Y* &lt; *y*<sub>*c*</sub></sub> = ∑<sub>*n* ≥ 0</sub>*ϕ*<sub>*n*</sub>(*y*<sub>*c*</sub>)×*H*<sub>*n*</sub>(*Y*)

where normalized Hermite polynomials are
$H\_n(n) = \\frac{1}{\\sqrt{n!}} \\frac{g^{(n)}(y)}{g(y)}$.

The computation of the coefficients *ϕ*<sub>*n*</sub>(*y*<sub>*c*</sub>)
gives

-   *ϕ*<sub>0</sub>(*y*<sub>*c*</sub>)=*g*(*y*<sub>*c*</sub>)+*y*<sub>*c*</sub> × *G*(*y*<sub>*c*</sub>),

-   *ϕ*<sub>1</sub>(*y*<sub>*c*</sub>)=*G*(*y*<sub>*c*</sub>)−1,

-   $\\phi\_n(y\_c) = g(y\_c) \\frac{H\_{n-2}(y\_c)}{\\sqrt{n\\times(n-1)}}$
    for *n* &gt; 1.

<!-- -->

    np <- 1000
    nh <- 1000
    yc <- qnorm(0.25)
    yy <- qnorm(p = (1:np)/(np+1))
    hm <- GM_hermite(y = yy, nh = nh-1)

Floor variable Z = max(yc, Y)

    psi = hermiteCoefLower(yc, nh)
    zz <- hm %*% psi
    ylim = c(-3.5, 3.5)
    p = plot.XY(yy, zz, join=TRUE, color = "red")
    p = plot.XY(yy, pmax(yc, yy), color = "blue", join=TRUE, linetype = "dashed", padd=p)
    p = plot.decoration(p, xlab = "Y", ylab = "Z", 
                        title = paste0("cut = ",round(yc,2), " nh = ", nh))
    p = plot.geometry(p, ylim=ylim)
    p


    #legend("topleft", legend = c("theoretical", "Hermit. pol."), col = c("blue", "red"), lty = #c(1,2))

QC Indicator

    psi = hermiteIndicatorLower(yc, nh)
    zz <- hm %*% psi
    p = plot.XY(yy, zz, join=TRUE, color = "red")
    p = plot.decoration(p, xlab = "Y", ylab = "Indicator", 
                        title = paste0("cut = ",round(yc,2), " nh = ", nh)) 
    p = plot.geometry(p, ylim = range(-0.5, 1.5))
    p = plot.XY(yy, as.numeric(yy >= yc), color = "blue", join=TRUE, linetype="dashed", padd=p)
    p


    #legend("bottomright", legend = c("theoretical", "Hermit. pol."), col = c("blue", "red"), lty = #c(1,2))
