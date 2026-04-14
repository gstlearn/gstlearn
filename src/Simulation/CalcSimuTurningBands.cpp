/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Model/Model.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/TurningBandDirection.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include <cmath>

namespace gstlrn
{

  CalcSimuTurningBands::CalcSimuTurningBands(
    Id nbsimu,
    Id nbtuba,
    bool flag_check,
    Id seed)
    : ACalcSimulation(nbsimu, seed)
    , _nbtuba(nbtuba)
    , _iattOut(-1)
    , _icase(0)
    , _flagCheck(flag_check)
    , _flagAllocationAlreadyDone(false)
    , _nameCoord()
    , _npointSimulated(0)
    , _field(0.)
    , _theta(0.)
    , _seedBands()
    , _codirs()
    , _modelLocal(nullptr)
  {
  }

  CalcSimuTurningBands::~CalcSimuTurningBands() {}

  bool CalcSimuTurningBands::_resizeTB()
  {
    auto nbsimu = getNbSimu();
    auto nbtuba = getNBtuba();
    if (nbtuba <= 0)
    {
      messerr(" The number of Turning Bands must be positive");
      return false;
    }

    _codirs.clear();
    _seedBands.clear();

    // The class is only partially defined
    if (nbsimu > 0 && nbtuba > 0)
    {
      auto nvar = _getNVar();
      auto ncova = _getNCov();

      /* Allocate the structures for the seeds */

      Id size = nvar * ncova * nbtuba * nbsimu;
      _seedBands.resize(size, 0);

      /* Allocate the structures for the directions */

      Id nbands = nbsimu * nbtuba * ncova;
      _codirs.resize(nbands);
      for (Id i = 0; i < nbands; i++) _codirs[i] = TurningBandDirection();
    }

    return true;
  }

  Id
    CalcSimuTurningBands::_getAddressBand(Id ivar, Id is, Id ib, Id isimu) const
  {
    auto nvar = _getNVar();
    auto ncova = _getNCov();
    auto nbtuba = getNbtuba();
    return ivar + nvar * ((is) + ncova * ((ib) + nbtuba * (isimu)));
  }

  void
    CalcSimuTurningBands::_setSeedBand(Id ivar, Id is, Id ib, Id isimu, Id seed)
  {
    auto iad = _getAddressBand(ivar, is, ib, isimu);
    _seedBands[iad] = seed;
  }

  Id CalcSimuTurningBands::_getSeedBand(Id ivar, Id is, Id ib, Id isimu) const
  {
    auto iad = _getAddressBand(ivar, is, ib, isimu);
    return _seedBands[iad];
  }

  void CalcSimuTurningBands::_dumpSeeds() const
  {
    auto nvar = _getNVar();
    auto ncova = _getNCov();
    auto nbsimu = getNbSimu();
    auto nbtuba = getNBtuba();

    mestitle(1, "Seeds");
    for (Id isimu = 0; isimu < nbsimu; isimu++)
      for (Id ivar = 0; ivar < nvar; ivar++)
        for (Id is = 0; is < ncova; is++)
          for (Id ib = 0; ib < nbtuba; ib++)
          {
            auto iad = _getAddressBand(ivar, is, ib, isimu);
            message(
              "Var=%d Simu=%d Is=%d Ib=%d iad=%d : %d\n", ivar, isimu, is, ib,
              iad, _seedBands[iad]);
          }
  }

  Id CalcSimuTurningBands::_getIBS(Id isimu, Id is, Id ib) const
  {
    auto ncova = _getNCov();
    auto nbtuba = getNBtuba();
    return ib + nbtuba * (is + ncova * isimu);
  }

  /****************************************************************************/
  /*!
   **  Generate directions according to Van Der Corput algorithm.
   **  The count of directions returned is the product of nbtuba by the
   **  number of basic structures
   **
   ** \param[in]  dbout   Db structure
   **
   *****************************************************************************/
  Id CalcSimuTurningBands::_generateDirections(const Db* dbout)
  {
    double x[2];
    auto ndim = _getNDim();
    auto ncova = _getNCov();
    auto nbsimu = getNbSimu();
    auto nbands = getNDirs();
    auto nbtuba = getNbtuba();

    /* Loop on the directions */

    for (Id ibs = 0; ibs < nbands; ibs++)
    {

      /* Decomposition according to basis 2 then 3 */

      for (Id id = 0; id < 2; id++)
      {
        Id n = 1 + ibs;
        Id p = id + 2;
        x[id] = 0;
        double d = id + 2;
        while (n > 0)
        {
          x[id] += (n % p) / d;
          d *= p;
          n /= p;
        }
      }

      /* Direction coefficients */

      double sqr = sqrt(1. - x[1] * x[1]);
      _setCodirAng(ibs, 0, cos(2. * GV_PI * x[0]) * sqr);
      _setCodirAng(ibs, 1, sin(2. * GV_PI * x[0]) * sqr);
      _setCodirAng(ibs, 2, x[1]);
      _setCodirTmin(ibs, MAXIMUM_BIG);
      _setCodirTmax(ibs, MINIMUM_BIG);
      _setCodirScale(ibs, 1.);
    }

    /* Random rotation of the directions */

    double r = 0.;
    double axyz[3];
    for (Id i = 0; i < 3; i++)
    {
      axyz[i] = law_gaussian();
      r += axyz[i] * axyz[i];
    }
    r = sqrt(r);
    for (Id i = 0; i < 3; i++) axyz[i] /= r;
    double theta = 2. * GV_PI * law_uniform(0., 1.);
    _rotateDirections(axyz, theta);

    /* Take the anisotropy into account */

    for (Id isimu = 0; isimu < nbsimu; isimu++)
      for (Id is = 0; is < ncova; is++)
        for (Id ib = 0; ib < nbtuba; ib++)
        {
          Id ibs = _getIBS(isimu, is, ib);
          const CovAniso* cova = _modelLocal->getCovAniso(is);

          // If the covariance has no Range (i.e. Nugget Effect), the rest is non-sense.
          // Nevertheless this code is maintained in order not to disorganize
          // the possible drawing of random numbers.
          if (!cova->hasRange()) continue;

          if (cova->getFlagAniso())
          {
            VectorDouble ranges = cova->getScales();
            double scale = 0.;
            for (Id i = 0; i < 3; i++)
            {
              double val = 0.;
              if (cova->getFlagRotation())
                for (Id j = 0; j < 3; j++)
                {
                  double rot = (i == j) ? 1. : 0.;
                  if (i < ndim && j < ndim)
                    rot = cova->getAnisoRotMat().getValue(i, j);
                  double range = 0.;
                  if (j < ndim) range = ranges[j];
                  if (range > 0.) val += (_getCodirAng(ibs, j) * rot / range);
                }
              else
              {
                double range = 0.;
                if (i < ndim) range = ranges[i];
                if (range > 0.) val += (_getCodirAng(ibs, i) / range);
              }
              scale += val * val;
              axyz[i] = val;
            }
            _setCodirScale(ibs, 1. / sqrt(scale));
            for (Id i = 0; i < 3; i++)
              _setCodirAng(ibs, i, axyz[i] * _getCodirScale(ibs));
          }
          else
          {
            _setCodirScale(ibs, cova->getScaleIso());
          }

          if (dbout->isGrid())
          {
            const auto* dbgrid = dynamic_cast<const DbGrid*>(dbout);
            double t00 = _codirs[ibs].projectGrid(dbgrid, 0, 0, 0);
            _setCodirT00(ibs, t00);
            _setCodirDXP(ibs, _codirs[ibs].projectGrid(dbgrid, 1, 0, 0) - t00);
            _setCodirDYP(ibs, _codirs[ibs].projectGrid(dbgrid, 0, 1, 0) - t00);
            _setCodirDZP(ibs, _codirs[ibs].projectGrid(dbgrid, 0, 0, 1) - t00);

            if (cova->getType() == ECov::SPHERICAL
                || cova->getType() == ECov::CUBIC)
            {
              double scale = _getCodirScale(ibs);
              _setCodirT00(ibs, _getCodirT00(ibs) / scale);
              _setCodirDXP(ibs, _getCodirDXP(ibs) / scale);
              _setCodirDYP(ibs, _getCodirDYP(ibs) / scale);
              _setCodirDZP(ibs, _getCodirDZP(ibs) / scale);
            }
          }
        }

    return 0;
  }

  void CalcSimuTurningBands::_dumpBands() const
  {
    auto nbands = getNDirs();
    bool flagGrid = getDbout()->isGrid();
    for (Id ibs = 0; ibs < nbands; ibs++)
    {
      message("- Band %d/%d\n", ibs + 1, nbands);
      _codirs[ibs].dump(flagGrid);
    }
  }

  /*****************************************************************************/
  /*!
   **  Perform the rotation of a set of normalized direction
   **  coefficients
   **
   ** \param[in]  a      Rotation direction
   ** \param[in]  theta  Rotation angle
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_rotateDirections(double a[3], double theta)
  {
    double dirs[3];
    double ct = cos(theta);
    double st = sin(theta);
    auto nbands = getNDirs();

    /* Loop on the direction coefficients */

    for (Id ibs = 0; ibs < nbands; ibs++)
    {
      for (Id idir = 0; idir < 3; idir++) dirs[idir] = _getCodirAng(ibs, idir);
      GH::rotationGetRandomDirection(ct, st, a, dirs);
      for (Id idir = 0; idir < 3; idir++) _setCodirAng(ibs, idir, dirs[idir]);
    }
  }

  /****************************************************************************/
  /*!
   **  Calculates the data extension for a set of turning bands
   **
   ** \param[in]  db      Db structure
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_minmax(const Db* db)
  {
    double tt;
    if (db == nullptr) return;
    auto nbands = getNDirs();

    if (db->isGrid())
    {
      const auto* dbgrid = dynamic_cast<const DbGrid*>(db);
      Id nx = (dbgrid->getNDim() >= 1) ? dbgrid->getNX(0) : 1;
      Id ny = (dbgrid->getNDim() >= 2) ? dbgrid->getNX(1) : 1;
      Id nz = (dbgrid->getNDim() >= 3) ? dbgrid->getNX(2) : 1;

      /* Case when the data obeys to a grid organization */
      /* This test is programmed for 3-D (maximum) grid  */
      /* as the Turning Bands method is limited to 3-D   */

      for (Id ibs = 0; ibs < nbands; ibs++)
      {
        for (Id iz = 0; iz < 2; iz++)
          for (Id iy = 0; iy < 2; iy++)
            for (Id ix = 0; ix < 2; ix++)
            {
              tt = _codirs[ibs].projectGrid(
                dbgrid, ix * (nx - 1), iy * (ny - 1), iz * (nz - 1));
              if (tt < _getCodirTmin(ibs)) _setCodirTmin(ibs, tt);
              if (tt > _getCodirTmax(ibs)) _setCodirTmax(ibs, tt);
              double delta = _getCodirTmax(ibs) - _getCodirTmin(ibs);
              if (_field < delta) _field = delta;
            }
      }
    }
    else
    {

      /* Case of an isolated set of data */

      for (Id iech = 0; iech < db->getNSample(); iech++)
      {
        if (!db->isActive(iech)) continue;
        for (Id ibs = 0; ibs < nbands; ibs++)
        {
          tt = _codirs[ibs].projectPoint(db, iech);
          if (tt < _getCodirTmin(ibs)) _setCodirTmin(ibs, tt);
          if (tt > _getCodirTmax(ibs)) _setCodirTmax(ibs, tt);
          double delta = _getCodirTmax(ibs) - _getCodirTmin(ibs);
          if (_field < delta) _field = delta;
        }
      }
    }

    _npointSimulated += db->getNSample();
  }

  /****************************************************************************/
  /*!
   **  Calculate the Poisson intensity for the generation
   **  of the Wiener-Levy along the line
   **
   ** \remark  The average number of points per band is calculated:
   ** \remark  - so as too have in average one band between two target
   ** \remark    points (either data or target)
   ** \remark  - to have a number of Poisson points per band lying within
   ** \remark    [nmini; nmaxi]
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_setDensity()

  {
    auto nbtuba = getNbtuba();
    Id nmini = 5;
    Id nmaxi = 5000;

    Id naverage = _npointSimulated / nbtuba;
    if (naverage < nmini) naverage = nmini;
    if (naverage > nmaxi) naverage = nmaxi;
    double scale = _field / naverage;
    _theta = 1. / scale;
  }

  /****************************************************************************/
  /*!
   **  Particular case of the stable model. It must be turned into:
   **  - Exponential : when param is too close to 1
   **  - Gaussian    : when param is too close to 2
   **
   **  Particular case of the Matern model. It must be turned into:
   **  - Exponential : when param is too close to 0.5
   **
   ** \return  The modified type
   **
   *****************************************************************************/
  ECov CalcSimuTurningBands::_particularCase(const CovAniso* cova, double eps)
  {
    double param = cova->getParam();
    const ECov& type = cova->getType();

    switch (type.toEnum())
    {
      case ECov::E_STABLE:
        if (ABS(param - 1.) < eps) return (ECov::EXPONENTIAL);
        if (ABS(param - 2.) < eps) return (ECov::GAUSSIAN);
        return (ECov::STABLE);
        break;

      case ECov::E_MATERN:
        if (ABS(param - 0.5) < eps) return (ECov::EXPONENTIAL);
        break;

      default: break;
    }
    return (type);
  }

  /****************************************************************************/
  /*!
   **  Initialize the array of seeds for the generation of a simulation
   **  using the Turning Bands method
   **
   ** \return  Error return code : 1 for problem; 0 otherwise
   **
   *****************************************************************************/
  Id CalcSimuTurningBands::_initializeSeedBands()
  {
    _setDensity();
    auto ncova = _getNCov();
    auto nvar = _getNVar();
    auto nbsimu = getNbSimu();
    auto nbtuba = getNbtuba();
    TurningBandOperate operTB;
    double correc;

    // Store the initial seed
    Id mem_seed = law_get_random_seed();

    //
    // Important remark: the order of the following loops must not be modified
    // in order to keep the same seeds and not modify the results
    //

    // Loop for fixing the seed for each Simulation / Variable / Covariance / Band
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id isimu = 0; isimu < nbsimu; isimu++)
        for (Id is = 0; is < ncova; is++)
          for (Id ib = 0; ib < nbtuba; ib++)
          {
            Id ibs = _getIBS(isimu, is, ib);
            ECov type = _particularCase(_modelLocal->getCovAniso(is));
            operTB.reset();
            _setSeedBand(ivar, is, ib, isimu, law_get_random_seed());

            Id optionSpectral = _getCorrec(type, is, ibs, operTB, correc);
            if (optionSpectral > 0) continue;
            messerr(
              "The structure (%s) cannot be simulated", type.getDescr().data());
            messerr("using the Turning Bands algorithm");
            return 1;
          }

    // Set the initial seed back
    law_set_random_seed(mem_seed);
    return 0;
  }

  /*****************************************************************************/
  /*!
   **  Generate a migration process
   **
   ** \param[in]  ibs   Rank of the turning band
   ** \param[in]  is    Rank of the structure
   ** \param[in]  scale scale of the exponential
   ** \param[in]  eps   Epsilon value
   **
   ** \param[out] operTB TurningBandOperate structure
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_migrationInit(
    Id ibs,
    Id is,
    double scale,
    TurningBandOperate& operTB,
    double eps)
  {
    DECLARE_UNUSED(is);
    static double vexp1 = 0.1;
    static double vexp2 = 0.1967708298;

    double tmin = _getCodirTmin(ibs);
    double tmax = _getCodirTmax(ibs);

    /* Initializations */

    double delta = tmax - tmin;
    if (scale < delta * eps)
    {
      Id count = static_cast<Id>(ceil(delta / eps));
      for (Id i = 0; i < count; i++) operTB.pushT(law_gaussian());
    }
    else
    {
      double value;
      value = tmin + scale * log(law_uniform(0., 1.));
      operTB.pushT(value);
      value = tmin - scale * log(law_uniform(0., 1.));
      operTB.pushT(value);

      while (value <= tmax)
      {
        value -= scale * log(law_uniform(0., 1.));
        operTB.pushT(value);
      }
    }

    double vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);
    operTB.setVexp(vexp);
  }

  /*****************************************************************************/
  /*!
   **  Generate a dilution process
   **
   ** \param[in]  ibs Rank of the Turning Band
   ** \param[in]  is  Rank of the covariance
   **
   ** \param[out] operTB TurningBandOperate structure
   **
   ** \return Correction factor
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_dilutionInit(
    Id ibs,
    Id is,
    TurningBandOperate& operTB)
  {
    double scale = _getCodirScale(ibs);
    double tmin = _getCodirTmin(ibs);
    double tmax = _getCodirTmax(ibs);
    double tdeb = tmin - scale * law_uniform(0., 1.);

    Id count = 0;
    while (tdeb + count * scale <= tmax)
    {
      double value = (law_uniform(0., 1.) < 0.5) ? -1. : 1.;
      operTB.pushT(value);
      count++;
    }

    operTB.setTdeb(tdeb);

    ECov type = _modelLocal->getCovType(is);
    double correc;
    switch (type.toEnum())
    {
      case ECov::E_SPHERICAL: correc = sqrt(3.); break;

      case ECov::E_CUBIC: correc = sqrt(840.); break;

      default: correc = TEST; break;
    }

    return correc;
  }

  /*****************************************************************************/
  /*!
   **  Initiate the spectral method
   **
   ** \param[in]  ibs    Rank of the Turning Band
   ** \param[in]  is     Rank of the covariance
   ** \param[in]  operTB TurningBandOperate structure
   **
   ** \return Correction factor
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_spectralInit(
    Id ibs,
    Id is,
    TurningBandOperate& operTB)
  {
    double scale = _getCodirScale(ibs);
    double param = _modelLocal->getParam(is);
    ECov type = _modelLocal->getCovType(is);

    double val = 0.;
    double period = 0.;
    switch (type.toEnum())
    {
      case ECov::E_GAUSSIAN:
        for (Id i = 0; i < 3; i++)
        {
          val = law_gaussian();
          period += val * val;
        }
        period = sqrt(2. * period) / scale;
        break;

      case ECov::E_STABLE:
        for (Id i = 0; i < 3; i++)
        {
          val = law_gaussian();
          period += val * val;
        }
        scale = _getScale(param, scale);
        period = sqrt(2. * period) / scale;
        break;

      case ECov::E_SINCARD:
        val = (law_uniform(0., 1.) >= 0.5) ? 1 : -1;
        period = val / scale;
        break;

      case ECov::E_BESSELJ:
        val = law_beta1(1.5, param - 0.5);
        period = sqrt(val) / scale;
        break;

      case ECov::E_MATERN:
        param = sqrt(2. * law_gamma(param));
        for (Id i = 0; i < 3; i++)
        {
          val = law_gaussian();
          period += val * val;
        }
        period = sqrt(period) / (param * scale);
        break;

      default: break;
    }

    operTB.setOmega(period);
    operTB.setPhi(2. * GV_PI * law_uniform(0., 1.));
    operTB.setOffset(0.);
    return sqrt(2.);
  }

  /****************************************************************************/
  /*!
   **  Calculate the scale for 1D process for the stable model
   **
   ** \return  Scale parameter of the 1D process to be simulated
   **
   ** \param[in]  alpha       Third parameter of the stable covariance model
   ** \param[in]  scale       Scale parameter of the model
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_getScale(double alpha, double scale)

  {
    if (alpha < 1) return scale / law_stable_standard_abgd(alpha);
    return scale / sqrt(law_stable_standard_abgd(alpha / 2.));
  }

  /****************************************************************************/
  /*!
   **  Calculate the scale for 1D process for the Matern model (param<0.5)
   **
   ** \return  Scale parameter of the 1D process to be simulated (param<0.5)
   **
   ** \param[in]  param       Third parameter of the Matern covariance model
   ** \param[in]  scale       Scale parameter of the model
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_getScaleKB(double param, double scale)
  {
    return (scale * sqrt(law_beta1(param, 0.5 - param)));
  }

  /*****************************************************************************/
  /*!
   **  Generate the 1D stochastic process
   **
   ** \param[in]  ibs     Rank of the turning band
   ** \param[in]  is      Rank of the covariance
   **
   ** \param[out] operTB  TurningBandOperate structure
   **
   ** \return Value of theta_alpha,3(R)
   **
   **  \remark Y_alpha_1=theta_alpha_1(R)cos(2pi R.x+phi)
   **  \remark used to simulate a GRF with a power semi-variogram h^alpha
   **  \remark according to the method proposed in
   **  \remark Emery, X. and Lantuejoul, C. (2008)
   **  \remark A spectral approach to simulating intrinsec random fields with power
   **  \remark and spline generalized covariance.
   **  \remark In Computational Geosciences 12:121-132
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_power1DInit(
    Id ibs,
    Id is,
    TurningBandOperate& operTB)
  {
    double R, theta_1;
    static double log3s2, log1s2, logap1, logap1s2, logap3s2, as2, coeff,
      coeff3;
    static double alpha_mem = -1.;

    double param = _modelLocal->getParam(is);
    if (ibs == 0 || !isEqual(param, alpha_mem))
    {
      double scale = _getCodirScale(ibs);
      as2 = param / 2.;
      log3s2 = loggamma(1.5);
      log1s2 = loggamma(0.5);
      logap1 = loggamma(1.0 + param);
      logap1s2 = loggamma(0.5 + as2);
      logap3s2 = loggamma(1.5 + as2);
      coeff = 2. * sqrt(exp(logap1) / pow(2. * GV_PI, param)) * pow(scale, as2);
      coeff3 = sqrt(exp(log3s2 + logap1s2 - log1s2 - logap3s2));
      alpha_mem = param;
    }

    R = law_beta2(1 - as2, as2);
    theta_1 = coeff * sqrt((R + 1.) / pow(R, as2 + 1));
    double phi = 2. * GV_PI * law_uniform(0., 1.);

    operTB.setOmega(2. * GV_PI * R);
    operTB.setPhi(phi);
    operTB.setOffset(cos(phi));

    return theta_1 / coeff3;
  }

  /*****************************************************************************/
  /*!
   **  Generate the 1D stochastic process for spline covariance
   **
   ** \param[in]  ibs     Rank of the turning band
   ** \param[in]  k       power of the variogram h^(2k)log(h)
   **
   ** \param[out] operTB  TurningBandOperate structure
   **
   ** \return    value of xi_2k,3(R)
   **
   **  \remark Compute the random elements
   **  \remark used to simulate a GRF with a power semi-variogram GC h^(2k)log(h)
   **  \remark according to the method proposed in
   **  \remark Emery, X. and Lantuejoul, C. (2008)
   **  \remark A spectral approach to simulating intrinsec random fields with power
   **  \remark and spline generalized covariance.
   **  \remark In Computational Geosciences 12:121-132
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_spline1DInit(
    Id ibs,
    Id k,
    TurningBandOperate& operTB)
  {
    double R;
    Id twokm1;
    static double twoPI, twokp1s2, log3s2, logkp3s2, logkp1, coeff;
    static Id twok = 0;
    static Id k_mem = -1;
    double scale = _getCodirScale(ibs);

    if (ibs == 0 || k != k_mem)
    {
      twok = 2 * k;
      twokm1 = twok - 1;
      twokp1s2 = twok + 1. / 2;
      twoPI = 2. * GV_PI;
      log3s2 = loggamma(1.5);
      logkp3s2 = loggamma(k + 1.5);
      logkp1 = loggamma(k + 1.);
      coeff = sqrt(2 * exp(logkp3s2 + logkp1 - log3s2) / pow(GV_PI, twokm1));
      k_mem = k;
    }

    R = law_beta2(1. / 2, 1. / 2);
    double phi = twoPI * law_uniform(0., 1.);

    operTB.setOmega(twoPI * R * pow(scale, twok));
    operTB.setPhi(phi);
    operTB.setOffset(cos(phi));

    return coeff * sqrt((R + 1.) / pow(R, twokp1s2));
  }

  /*****************************************************************************/
  /*!
   **  Generates the process constituted by independent gaussian
   **  variables along a 1D Poisson process. The process consists in
   **  the integration(s) of the previous process Perform the core
   **  allocation
   **
   ** \param[in]  ibs    Rank of the turning band
   ** \param[in]  is     Rank of the covariance
   **
   ** \param[out] operTB TruningBandOperate structure
   **
   ** \return Correction factor
   **
   ** \remark  This procedure allocates memory that should be freed
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_irfProcessInit(
    Id ibs,
    Id is,
    TurningBandOperate& operTB)
  {
    ECov type = _modelLocal->getCovType(is);
    double delta;

    Id level = -1;
    if (type == ECov::LINEAR) level = 0;
    if (type == ECov::ORDER1_GC) level = 0;
    if (type == ECov::ORDER3_GC) level = 1;
    if (type == ECov::ORDER5_GC) level = 2;
    auto nt = operTB.getTsize();
    const auto& t = operTB.getT();

    /* Generation of the Wiener-Levy process and its integrations */

    double val0 = 0.;
    double val1 = 0.;
    double val2 = 0.;
    if (level >= 0) operTB.pushV0(val0);
    if (level >= 1) operTB.pushV1(val1);
    if (level >= 2) operTB.pushV2(val2);

    for (Id i = 1; i < nt; i++)
    {
      val0 += law_gaussian();
      operTB.pushV0(val0);

      if (level == 0) continue;
      delta = t[i] - t[i - 1];
      val1 += val0 * delta;
      operTB.pushV1(val1);

      if (level == 1) continue;
      val2 += val1 * delta + val0 * delta * delta / 2.;
      operTB.pushV2(val2);
    }

    double scale = _getCodirScale(ibs);
    double theta1 = 1. / _theta;
    return _irfCorrec(type, theta1, scale);
  }

  void CalcSimuTurningBands::_spreadRegularOnGrid(
    const DbGrid* dbgrid,
    const CovAniso* cova,
    Id ibs,
    double correc,
    TurningBandOperate& operTB,
    const VectorBool& activeArray,
    VectorDouble& tabvar)
  {
    double t0y, t0z, t0;

    Id ndim = dbgrid->getNDim();
    Id nx = (ndim >= 1) ? dbgrid->getNX(0) : 1;
    Id ny = (ndim >= 2) ? dbgrid->getNX(1) : 1;
    Id nz = (ndim >= 3) ? dbgrid->getNX(2) : 1;

    double t00 = _getCodirT00(ibs);
    double dxp = _getCodirDXP(ibs);
    double dyp = _getCodirDYP(ibs);
    double dzp = _getCodirDZP(ibs);

    t0z = t00;
    Id ind = 0;
    for (Id iz = 0; iz < nz; iz++)
    {
      t0y = t0z;
      t0z += dzp;
      for (Id iy = 0; iy < ny; iy++)
      {
        t0 = t0y;
        t0y += dyp;
        for (Id ix = 0; ix < nx; ix++)
        {
          if (activeArray[ind])
            tabvar[ind] += correc * cova->simulateTurningBand(t0, operTB);
          t0 += dxp;

          ind++;
        }
      }
    }
  }

  void CalcSimuTurningBands::_spreadSpectralOnGrid(
    const DbGrid* dbgrid,
    const CovAniso* cova,
    Id ibs,
    double correc,
    TurningBandOperate& operTB,
    const VectorBool& activeArray,
    VectorDouble& tabvar)
  {
    double c1, s1, c0x, s0x, c0y, s0y, c0z, s0z, cxp, sxp, cyp, syp, czp, szp;
    Id ndim = dbgrid->getNDim();
    Id nx = (ndim >= 1) ? dbgrid->getNX(0) : 1;
    Id ny = (ndim >= 2) ? dbgrid->getNX(1) : 1;
    Id nz = (ndim >= 3) ? dbgrid->getNX(2) : 1;
    _getOmegaPhi(ibs, operTB, &cxp, &sxp, &cyp, &syp, &czp, &szp, &c0z, &s0z);

    Id ind = 0;
    for (Id iz = 0; iz < nz; iz++)
    {
      c0y = c0z;
      s0y = s0z;
      c1 = c0z * czp - s0z * szp;
      s1 = s0z * czp + c0z * szp;
      c0z = c1;
      s0z = s1;
      for (Id iy = 0; iy < ny; iy++)
      {
        c0x = c0y;
        s0x = s0y;
        c1 = c0y * cyp - s0y * syp;
        s1 = s0y * cyp + c0y * syp;
        c0y = c1;
        s0y = s1;
        for (Id ix = 0; ix < nx; ix++)
        {
          if (activeArray[ind])
            tabvar[ind] += correc * cova->simulateTurningBand(c0x, operTB);
          c1 = c0x * cxp - s0x * sxp;
          s1 = s0x * cxp + c0x * sxp;
          c0x = c1;
          s0x = s1;

          ind++;
        }
      }
    }
  }

  void CalcSimuTurningBands::_spreadRegularOnPoint(
    const Db* db,
    const CovAniso* cova,
    Id ibs,
    double correc,
    TurningBandOperate& operTB,
    const VectorBool& activeArray,
    VectorDouble& tabvar)
  {
    double t0;
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!activeArray[iech]) continue;
      t0 = _codirs[ibs].projectPoint(db, iech);
      tabvar[iech] += correc * cova->simulateTurningBand(t0, operTB);
    }
  }

  void CalcSimuTurningBands::_spreadSpectralOnPoint(
    const Db* db,
    const CovAniso* cova,
    Id ibs,
    double correc,
    TurningBandOperate& operTB,
    const VectorBool& activeArray,
    VectorDouble& tabvar)
  {
    double t0;
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!activeArray[iech]) continue;
      t0 = _codirs[ibs].projectPoint(db, iech);
      tabvar[iech] += correc * cova->simulateTurningBand(t0, operTB);
    }
  }

  /**
   * @brief Compute one non conditional simulation on the samples of Db using Turning Bands method
   *
   * @param db The target Db
   * @param isimu The rank of the simulation to compute
   * @param activeArray Array indicating active samples
   * @param tab Array to store the simulation values for all bands
   */
  Id CalcSimuTurningBands::_compute(
    Db* db,
    Id isimu,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    auto* dbgrid = dynamic_cast<DbGrid*>(db);
    bool flagGrid = (dbgrid != nullptr);

    VectorBool activeLoc; // Not used
    VectorVectorDouble tabLoc;
    _allocateForOneSimulation(db, getNVar(), activeLoc, tabLoc, false);

    // Seed must be memorized as it will be modified when generating each band
    // within _computeGrid() or _computePoint() facility.
    Id mem_seed = law_get_random_seed();

    for (Id icov = 0, ncova = _getNCov(); icov < ncova; icov++)
    {
      const CovAniso* cova = _modelLocal->getCovAniso(icov);
      ECov type = _particularCase(cova);
      if (type == ECov::NUGGET) continue;

      // Blank out the array 'tab'
      tabLoc.fillWith(0.);

      // Evaluate the multivariate simulation on the target samples for current structure
      if (flagGrid)
        _computeGrid(dbgrid, cova, type, isimu, icov, activeArray, tabLoc);
      else
        _computePoint(db, cova, type, isimu, icov, activeArray, tabLoc);

      // Cumulate structures for current simulation
      _scaleResults(db, cova, activeArray, tabLoc, tab);
    }

    // Set the initial seed back
    law_set_random_seed(mem_seed);
    return 0;
  }

  /**
   * @brief Save the multivariate simulation result into the Db after:
   * @brief - multiplying by the sill matrix
   * @brief - adding to existing values
   *
   * @param db Db where the result is stored
   * @param cova Covariance where to read the AIC matrix
   * @param activeArray Array indicating active samples
   * @param tabLoc Array containing the non-conditional simulation values for all variables
   * @param tab   Array containing simulation values for all variables
   */
  void CalcSimuTurningBands::_scaleResults(
    Db* db,
    const CovBase* cova,
    const VectorBool& activeArray,
    const VectorVectorDouble& tabLoc,
    VectorVectorDouble& tab) const
  {
    auto nvar = getNVar();
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
      if (activeArray[iech])
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
          for (Id jvar = 0; jvar < nvar; jvar++)
            tab[jvar][iech] += tabLoc[ivar][iech] * cova->getAic(ivar, jvar);
      }
  }

  /**
   * @brief For all simulated values, apply the normation by the number of bands
   *
   * @param db Target Db
   * @param activeArray Array indicating active samples
   * @param tab Array containing simulation values for all bands
   */
  void CalcSimuTurningBands::_normalizeForBands(
    const Db* db,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    auto nech = db->getNSample();
    auto nvar = _getNVar();
    auto nbtuba = getNbtuba();
    double norme = sqrt(1. / nbtuba);

    for (Id iech = 0; iech < nech; iech++)
      if (activeArray[iech])
        for (Id ivar = 0; ivar < nvar; ivar++) tab[ivar][iech] *= norme;
  }

  Id CalcSimuTurningBands::_getCorrec(
    const ECov& type,
    Id is,
    Id ibs,
    TurningBandOperate& operTB,
    double& correc)
  {
    double scale;
    double param;
    double theta1;
    correc = 1.;
    switch (type.toEnum())
    {
      case ECov::E_NUGGET:
        // Next line is simply to let the random number cycle
        (void)law_gaussian();
        return 1;

      case ECov::E_STABLE:
        param = _modelLocal->getParam(is);
        if (param > 1)
        {
          correc = _spectralInit(ibs, is, operTB);
          return 1;
        }
        scale = _getScale(param, 2. * _getCodirScale(ibs));
        _migrationInit(ibs, is, scale, operTB);
        return 2;

      case ECov::E_MATERN:
        param = _modelLocal->getParam(is);
        if (param > 0.5)
        {
          correc = _spectralInit(ibs, is, operTB);
          return 1;
        }
        scale = _getScaleKB(param, _getCodirScale(ibs)) * 2;
        _migrationInit(ibs, is, scale, operTB);
        return 2;

      case ECov::E_EXPONENTIAL:
        scale = _getCodirScale(ibs);
        _migrationInit(ibs, is, 2. * scale, operTB);
        return 2;

      case ECov::E_SPHERICAL:
      case ECov::E_CUBIC: correc = _dilutionInit(ibs, is, operTB); return 2;

      case ECov::E_POWER: correc = _power1DInit(ibs, is, operTB); return 1;

      case ECov::E_SPLINE_GC: correc = _spline1DInit(ibs, 1, operTB); return 1;

      case ECov::E_GAUSSIAN:
      case ECov::E_SINCARD:
      case ECov::E_BESSELJ: correc = _spectralInit(ibs, is, operTB); return 1;

      case ECov::E_LINEAR:
      case ECov::E_ORDER1_GC:
      case ECov::E_ORDER3_GC:
      case ECov::E_ORDER5_GC:
        theta1 = 1. / _theta;
        _migrationInit(ibs, is, theta1, operTB);
        correc = _irfProcessInit(ibs, is, operTB);
        return 2;

      default: return 0;
    }
    return 0;
  }

  /*****************************************************************************/
  /*!
   **  Perform non-conditional simulations on a set of points using
   **  Turning Bands method.
   **
   ** \param[in]  db         Db structure
   ** \param[in]  cova       Covariance Anisotropy structure
   ** \param[in]  type       Covariance type
   ** \param[in]  isimu      Simulation index
   ** \param[in]  is         Covariance index
   ** \param[in]  activeArray  Array indicating active samples
   ** \param[out] tab        Array to store simulation values for one band
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_computePoint(
    Db* db,
    const CovAniso* cova,
    const ECov& type,
    Id isimu,
    Id is,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    TurningBandOperate operTB;
    double correc;

    for (Id ivar = 0, nvar = _getNVar(); ivar < nvar; ivar++)
      for (Id ib = 0, nbtuba = getNbtuba(); ib < nbtuba; ib++)
      {
        Id ibs = _getIBS(isimu, is, ib);
        operTB.reset();
        operTB.setScale(_getCodirScale(ibs));
        operTB.setFlagScaled(false);

        law_set_random_seed(_getSeedBand(ivar, is, ib, isimu));
        Id optionSpectral = _getCorrec(type, is, ibs, operTB, correc);

        // Spreading the values on the points within 'tab'
        if (optionSpectral == 1)
          _spreadSpectralOnPoint(
            db, cova, ibs, correc, operTB, activeArray, tab[ivar]);
        else
          _spreadRegularOnPoint(
            db, cova, ibs, correc, operTB, activeArray, tab[ivar]);
      }

    // Normalize by the count of bands
    _normalizeForBands(db, activeArray, tab);
  }

  /*****************************************************************************/
  /*!
   **  Perform non-conditional simulations on a grid using the
   **  Turning Bands method
   **
   ** \param[in]  dbgrid     Db Grid structure
   ** \param[in]  cova       Covariance Anisotropy structure
   ** \param[in]  type       Covariance type
   ** \param[in]  isimu      Simulation index
   ** \param[in]  is         Covariance index
   ** \param[in]  activeArray  Array indicating active samples
   ** \param[out] tab        Array to store simulation values for one band
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_computeGrid(
    DbGrid* dbgrid,
    const CovAniso* cova,
    const ECov& type,
    Id isimu,
    Id is,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    TurningBandOperate operTB;
    double correc;

    for (Id ivar = 0, nvar = _getNVar(); ivar < nvar; ivar++)
      for (Id ib = 0, nbtuba = getNbtuba(); ib < nbtuba; ib++)
      {
        Id ibs = _getIBS(isimu, is, ib);
        operTB.reset();
        operTB.setScale(_getCodirScale(ibs));
        operTB.setFlagScaled(true);

        law_set_random_seed(_getSeedBand(ivar, is, ib, isimu));
        Id optionSpectral = _getCorrec(type, is, ibs, operTB, correc);

        // Spreading the values on the grid within 'tab'
        if (optionSpectral == 1)
          _spreadSpectralOnGrid(
            dbgrid, cova, ibs, correc, operTB, activeArray, tab[ivar]);
        else
          _spreadRegularOnGrid(
            dbgrid, cova, ibs, correc, operTB, activeArray, tab[ivar]);
      }

    // Normalize by the count of bands
    _normalizeForBands(dbgrid, activeArray, tab);
  }

  /****************************************************************************/
  /*!
   **  Calculate the correction factor for IRF_k models
   **
   ** \return  Correction factor
   **
   ** \param[in]  type    type of polynomial generalized covariance
   ** \param[in]  theta1  Equal to inverse of theta value
   ** \param[in]  scale   Range of the model
   **
   *****************************************************************************/
  double CalcSimuTurningBands::_irfCorrec(
    const ECov& type,
    double theta1,
    double scale)
  {
    switch (type.toEnum())
    {
      case ECov::E_LINEAR:
      case ECov::E_ORDER1_GC: return sqrt((4. * theta1) / scale); break;

      case ECov::E_ORDER3_GC:
        return sqrt((48. * theta1) / scale) / scale;
        break;

      case ECov::E_ORDER5_GC:
        return sqrt((1440. * theta1) / scale) / scale / scale;
        break;

      default: break;
    }

    return TEST;
  }

  void CalcSimuTurningBands::_getOmegaPhi(
    Id ibs,
    TurningBandOperate& operTB,
    double* cxp,
    double* sxp,
    double* cyp,
    double* syp,
    double* czp,
    double* szp,
    double* c0z,
    double* s0z)
  {

    double omega = operTB.getOmega();
    double phi = operTB.getPhi();

    double dxp = _getCodirDXP(ibs);
    double dyp = _getCodirDYP(ibs);
    double dzp = _getCodirDZP(ibs);
    double t00 = _getCodirT00(ibs);

    *cxp = cos(omega * dxp);
    *sxp = sin(omega * dxp);
    *cyp = cos(omega * dyp);
    *syp = sin(omega * dyp);
    *czp = cos(omega * dzp);
    *szp = sin(omega * dzp);

    *c0z = cos(omega * t00 + phi);
    *s0z = sin(omega * t00 + phi);
  }

  bool CalcSimuTurningBands::_simulateTB()
  {
    // Dimension the Turning Bands environment
    if (!_resizeTB()) return false;

    // Generate the Bands directions
    if (_generateDirections(getDbout())) return false;

    // Calculate the extension of the field
    _minmax(getDbout());
    _minmax(getDbin());

    // Populate the bands
    if (_initializeSeedBands()) return false;

    // Factorize the matrix of sills
    _modelLocal->computeAic();

    return true;
  }

  /*****************************************************************************/
  /*!
   **  Add the contribution of the nugget effect to the non-conditional
   **  simulations
   **
   ** \param[in]  db         Db structure
   ** \param[in]  isimu      Index of the simulation
   ** \param[in]  activeArray  Array of active samples
   ** \param[out] tab        Array containing the non-conditional simulation values
   **
   *****************************************************************************/
  void CalcSimuTurningBands::_computeNugget(
    Db* db,
    Id isimu,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    /* Do nothing if there is no nugget effect in the model */
    if (!_modelLocal->hasNugget()) return;

    Id nech = db->getNSample();
    auto ncova = _getNCov();
    auto nvar = _getNVar();

    // Memorize the seed
    Id mem_seed = law_get_random_seed();

    /* Performing the simulation */

    for (Id is = 0; is < ncova; is++)
    {
      ECov type = _modelLocal->getCovType(is);
      if (type != ECov::NUGGET) continue;
      const CovAniso* cova = _modelLocal->getCovAniso(is);

      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        law_set_random_seed(_getSeedBand(ivar, is, 0, isimu));

        for (Id iech = 0; iech < nech; iech++)
        {
          if (!activeArray[iech]) continue;
          double nugget = law_gaussian();
          for (Id jvar = 0; jvar < nvar; jvar++)
            tab[jvar][iech] += nugget * cova->getAic(ivar, jvar);
        }
      }
    }

    // Set the initial seed back
    law_set_random_seed(mem_seed);
  }

  bool CalcSimuTurningBands::_run()
  {
    bool flag_cond = hasDbin(false);
    VectorVectorDouble tabOut;
    VectorVectorDouble tabIn;
    VectorBool activeOut;
    VectorBool activeIn;
    _allocateForOneSimulation(getDbout(), getNVar(), activeOut, tabOut);
    if (flag_cond)
      _allocateForOneSimulation(getDbin(), getNVar(), activeIn, tabIn);

    // Loop on the simulations
    for (Id isimu = 0, nbsimu = getNbSimu(); isimu < nbsimu; isimu++)
    {
      if (getVerbose()) message(">>> computing simulation %d\n", isimu + 1);
      tabOut.fillWith(0);
      if (flag_cond) tabIn.fillWith(0);

      // Non conditional simulations on the target points
      _compute(getDbout(), isimu, activeOut, tabOut);
      _correctMean(activeOut, tabOut);
      _computeNugget(getDbout(), isimu, activeOut, tabOut);
      saveResults(getDbout(), isimu, activeOut, tabOut);

      if (!flag_cond) continue;

      // Non conditional simulations on the data points
      _compute(getDbin(), isimu, activeIn, tabIn);
      _correctMean(activeIn, tabIn);
      _difference(getDbin(), isimu, activeIn, tabIn);
      saveResults(getDbin(), isimu, activeIn, tabIn);
    }

    // Conditional simulations
    if (flag_cond)
    {
      if (_conditionalKriging(getDbin(), getDbout())) return 1;
    }

    // Copy value from data to coinciding grid node
    if (flag_cond)
    {
      _updateDataToTarget(getDbin(), getDbout());
    }

    // Check the simulation at data location
    if (_flagCheck)
    {
      if (_checkGaussianDataToGrid(getDbin(), getDbout())) return false;
    }
    return true;
  }

  /****************************************************************************/
  /*!
   **  Perform the Simulation Process using the Turning Bands Method
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin       Input Db structure
   ** \param[in]  dbout      Output Db structure
   ** \param[in]  model      Model structure
   ** \param[in]  neigh      ANeigh structure
   ** \param[in]  icase      Case for PGS or -1
   ** \param[in]  flag_bayes 1 if the Bayes option is switched ON
   ** \param[in]  flag_pgs   1 if called from PGS
   ** \param[in]  flag_gibbs 1 if called from Gibbs
   ** \param[in]  flag_dgm   1 if the Discrete Gaussian Model is used
   **
   *****************************************************************************/
  Id CalcSimuTurningBands::simulate(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    Id icase,
    Id flag_bayes,
    bool flag_pgs,
    bool flag_gibbs,
    bool flag_dgm)
  {
    setDbin(dbin);
    setDbout(dbout);
    setModelGeneric(model);
    setNeigh(neigh);
    setIcase(icase);
    setFlagBayes(flag_bayes);
    setFlagPGS(flag_pgs);
    setFlagGibbs(flag_gibbs);
    setFlagDGM(flag_dgm);

    if (!run()) return 1;
    return 0;
  }

  /****************************************************************************/
  /*!
   **  Perform the (non-conditional) Simulation(s) using the Turning Bands Method
   **
   ** \return  Error return code
   **
   ** \param[in]  dbiso     Isovalues Db structure
   ** \param[in]  dbgrd     Gradient Db structure
   ** \param[in]  dbtgt     Tangent Db structure
   ** \param[in]  dbout     Output Db structure
   ** \param[in]  model     Model structure
   ** \param[in]  delta     Value of the increment
   **
   *****************************************************************************/
  Id CalcSimuTurningBands::simulatePotential(
    Db* dbiso,
    Db* dbgrd,
    Db* dbtgt,
    Db* dbout,
    ModelGeneric* model,
    double delta)
  {
    setDbout(dbout);
    setModelGeneric(model);
    if (getNbSimu() <= 0 || getNBtuba() <= 0)
    {
      messerr("You must define 'nbsimu', 'nbtuba' and the 'model' beforehand");
      return 1;
    }
    law_set_random_seed(getSeed());

    /* Processing the Turning Bands algorithm */

    if (_generateDirections(dbout)) return 1;
    _minmax(dbout);
    _minmax(dbiso);
    _minmax(dbgrd);
    _minmax(dbtgt);
    if (_initializeSeedBands()) return 1;

    // Factorize the matrix of sills
    _modelLocal->computeAic();

    /* Non conditional simulations on the data points */
    if (dbiso != nullptr)
    {
      VectorVectorDouble tab;
      VectorBool activeArray;
      _allocateForOneSimulation(getDbout(), getNVar(), activeArray, tab);
      for (Id isimu = 0; isimu < getNbSimu(); isimu++)
        _compute(dbiso, isimu, activeArray, tab);
    }

    /* Non conditional simulations on the gradient points */
    if (dbgrd != nullptr)
    {
      VectorVectorDouble tab;
      VectorBool activeArray;
      _allocateForOneSimulation(getDbout(), 1, activeArray, tab);
      for (Id isimu = 0; isimu < getNbSimu(); isimu++)
        _computeGradient(dbgrd, isimu, delta, activeArray, tab);
    }

    /* Non conditional simulations on the tangent points */
    if (dbtgt != nullptr)
    {
      VectorVectorDouble tab;
      VectorBool activeArray;
      _allocateForOneSimulation(getDbout(), 1, activeArray, tab);
      for (Id isimu = 0; isimu < getNbSimu(); isimu++)
      {
        _computeTangent(dbtgt, isimu, delta, activeArray, tab);
      }
    }

    /* Non conditional simulations on the target samples */
    VectorVectorDouble tab;
    VectorBool activeArray;
    _allocateForOneSimulation(getDbout(), getNVar(), activeArray, tab);
    for (Id isimu = 0; isimu < getNbSimu(); isimu++)
    {
      _compute(dbout, isimu, activeArray, tab);

      /* Add the contribution of nugget effect (optional) */
      _computeNugget(dbout, isimu, activeArray, tab);
    }

    return 0;
  }

  bool CalcSimuTurningBands::_check()
  {
    if (!ACalcSimulation::_check()) return false;

    if (!hasDbout()) return false;
    if (!hasModelGeneric()) return false;
    if (hasDbin(false))
    {
      if (!hasNeigh()) return false;
    }

    auto ndim = _getNDim();
    if (ndim > 3)
    {
      messerr("The Turning Band Method is not a relevant simulation model");
      messerr("for this Space Dimension (%d)", ndim);
      return false;
    }

    // Check that the model is of type 'Model'
    _modelLocal = dynamic_cast<Model*>(getModelGeneric());
    if (_modelLocal == nullptr)
    {
      messerr("The model must be of type Model not ModelGeneric)");
      return false;
    }
    Id ncova = _modelLocal->getNCov();
    if (ncova <= 0)
    {
      messerr("The Model should contain some valid covariances");
      return false;
    }
    if (!_modelLocal->isValidForSimulation(ESimuType::TB)) return false;

    if (getNBtuba() <= 0)
    {
      messerr("You must define 'nbsimu' and 'nbtuba'");
      return 1;
    }

    if (_getFlagDGM())
    {
      if (!getDbout()->isGrid())
      {
        messerr("For DGM option, the argument 'dbout'  should be a Grid");
        return false;
      }
      if (!_modelLocal->hasAnam())
      {
        messerr("For DGM option, the Model must have an Anamorphosis attached");
        return false;
      }
      if (!_modelLocal->isChangeSupportDefined())
      {
        messerr("DGM option requires a Change of Support to be defined");
        return false;
      }
    }
    return true;
  }

  bool CalcSimuTurningBands::_preprocess()
  {
    if (!ACalcSimulation::_preprocess()) return false;

    // Prepare the seeds for the Bands
    law_set_random_seed(getSeed());
    if (!_simulateTB()) return false;

    auto nvar = _getNVar();
    auto nbsimu = getNbSimu();

    /* Add the attributes for storing the results */

    if (getDbin() != nullptr)
    {
      if (!_flagAllocationAlreadyDone)
      {
        Id iptr_in = _addVariableDb(1, 2, ELoc::SIMU, 0, nvar * nbsimu);
        if (iptr_in < 0) return false;
      }
    }

    if (!_flagAllocationAlreadyDone)
    {
      _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 0, nvar * nbsimu);
      if (_iattOut < 0) return false;
    }

    // Centering the Data (for DGM)

    if (_getFlagDGM())
    {
      // Centering (only if the output file is a Grid)
      auto* dbgrid = dynamic_cast<DbGrid*>(getDbout());
      if (dbgrid != nullptr)
      {
        // Duplicating the coordinate variable names before centering
        _nameCoord = getDbin()->getNamesByLocator(ELoc::X);
        if (_centerDataToGrid(dbgrid)) return false;
      }
    }

    return true;
  }

  bool CalcSimuTurningBands::_postprocess()
  {
    /* Free the temporary variables */
    _cleanVariableDb(2);

    // Clean variables created for Expansion

    if (_expandInformation(-1, ELoc::F)) return false;
    if (_expandInformation(-1, ELoc::NOSTAT)) return false;

    /* Set the error return flag */

    if (!_flagAllocationAlreadyDone)
      _renameVariable(
        2, VectorString(), ELoc::Z, _getNVar(), _iattOut, String(),
        getNbSimu());

    if (_getFlagDGM())
    {
      if (!_nameCoord.empty()) getDbin()->setLocators(_nameCoord, ELoc::X, 0);
    }

    return true;
  }

  void CalcSimuTurningBands::_rollback()
  {
    _cleanVariableDb(1);
  }

} // namespace gstlrn
