/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/TurningDirection.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Model/Model.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Covariances/CovAniso.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <math.h>

CalcSimuTurningBands::CalcSimuTurningBands(int nbsimu, int nbtuba, bool flag_check, int seed)
    : ACalcSimulation(nbsimu, seed),
      _nbtuba(nbtuba),
      _iattOut(-1),
      _icase(0),
      _flagCheck(flag_check),
      _flagBayes(false),
      _flagPGS(false),
      _flagGibbs(false),
      _flagDGM(false),
      _nameCoord(),
      _bayesMean(),
      _bayesCov(),
      _npointSimulated(0),
      _field(0.),
      _theta(0.),
      _seedBands(),
      _codirs()
{
}

CalcSimuTurningBands::~CalcSimuTurningBands()
{
}

bool CalcSimuTurningBands::_resize()
{
  int nbsimu = getNbSimu();
  int nbtuba = getNBtuba();
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
    int nvar  = _getNVar();
    int ncova = _getNCova();

    /* Allocate the structures for the seeds */

    int size = nvar * ncova * _nbtuba * nbsimu;
    _seedBands.resize(size,0);

    /* Allocate the structures for the directions */

    int nbands = nbsimu * _nbtuba * ncova;
    _codirs.resize(nbands);
    for (int i = 0; i < nbands; i++)
      _codirs[i] = TurningDirection();
  }

  return true;
}

int CalcSimuTurningBands::_getAddressBand(int ivar, int is, int ib, int isimu)
{
  int ncova = _getNCova();
  int nvar = _getNVar();
  return ivar+nvar*((is)+ncova*((ib)+_nbtuba*(isimu)));
}

void CalcSimuTurningBands::_setSeedBand(int ivar, int is, int ib, int isimu, int seed)
{
  int iad = _getAddressBand(ivar, is, ib, isimu);
  _seedBands[iad] = seed;
}

int CalcSimuTurningBands::_getSeedBand(int ivar, int is, int ib, int isimu)
{
  int iad = _getAddressBand(ivar, is, ib, isimu);
  return _seedBands[iad];
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
int CalcSimuTurningBands::_generateDirections(const Db* dbout)
{
  double x[2];
  int ndim = _getNDim();
  int ncova = _getNCova();
  int nbsimu = getNbSimu();
  int nbands = getNDirs();

  /* Loop on the directions */

  for (int ibs = 0; ibs < nbands; ibs++)
  {

    /* Decomposition according to basis 2 then 3 */

    for (int id = 0; id < 2; id++)
    {
      int n = 1 + ibs;
      int p = id + 2;
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
    _setCodirTmin(ibs, 1.e30);
    _setCodirTmax(ibs,-1.e30);
    _setCodirScale(ibs, 1.);
  }

  /* Random rotation of the directions */

  double r = 0.;
  double axyz[3];
  for (int i = 0; i < 3; i++)
  {
    axyz[i] = law_gaussian();
    r += axyz[i] * axyz[i];
  }
  r = sqrt(r);
  for (int i = 0; i < 3; i++)
    axyz[i] /= r;
  double theta = 2. * GV_PI * law_uniform(0., 1.);
  _rotateDirections(axyz, theta);

  /* Take the anisotropy into account */

  int ibs = 0;
  for (int isimu = 0; isimu < nbsimu; isimu++)
    for (int is = 0; is < ncova; is++)
      for (int ib = 0; ib < _nbtuba; ib++, ibs++)
      {
        const CovAniso* cova = getModel()->getCova(is);

        // If the covariance has no Range (i.e. Nugget Effect), the rest is non-sense.
        // Nevertheless this code is maintained to ensure in order not to disorganize
        // the possible drawing of random numbers.
        if (!cova->hasRange()) continue;

        if (cova->getFlagAniso())
        {
          VectorDouble ranges = cova->getScales();
          double scale = 0.;
          for (int i = 0; i < 3; i++)
          {
            double val = 0.;
            if (cova->getFlagRotation())
              for (int j = 0; j < 3; j++)
              {
                double rot = (i == j) ? 1. : 0.;
                if (i < ndim && j < ndim) rot = cova->getAnisoRotMat().getValue(i, j);
                double range = 0.;
                if (j < ndim) range = ranges[j];
                if (range > 0.)
                  val += (_getCodirAng(ibs, j) * rot / range);
              }
            else
            {
              double range = 0.;
              if (i < ndim) range = ranges[i];
              if (range > 0.) val += (_getCodirAng(ibs,i) / range);
            }
            scale += val * val;
            axyz[i] = val;
          }
          _setCodirScale(ibs, 1. / sqrt(scale));
          for (int i = 0; i < 3; i++)
            _setCodirAng(ibs, i, axyz[i] * _getCodirScale(ibs));
        }
        else
        {
          _setCodirScale(ibs, cova->getScale());
        }

        if (dbout->isGrid())
        {
          const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(dbout);
          double t00 = _codirs[ibs].projectGrid(dbgrid, 0, 0, 0);
          _setCodirT00(ibs, t00);
          _setCodirDXP(ibs,_codirs[ibs].projectGrid(dbgrid, 1, 0, 0) - t00);
          _setCodirDYP(ibs,_codirs[ibs].projectGrid(dbgrid, 0, 1, 0) - t00);
          _setCodirDZP(ibs,_codirs[ibs].projectGrid(dbgrid, 0, 0, 1) - t00);

          if (cova->getType() == ECov::SPHERICAL || cova->getType() == ECov::CUBIC)
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
  int nbands = getNDirs();

  /* Loop on the direction coefficients */

  for (int ibs = 0; ibs < nbands; ibs++)
  {
    for (int idir = 0; idir < 3; idir++)
      dirs[idir] = _getCodirAng(ibs, idir);
    GH::rotationGetDirection(ct, st, a, dirs);
    for (int idir = 0; idir < 3; idir++)
      _setCodirAng(ibs, idir, dirs[idir]);
  }
}

/****************************************************************************/
/*!
 **  Calculates the data extension for a set of turning bands
 **
 ** \param[in]  db      Db structure
 **
 *****************************************************************************/
void CalcSimuTurningBands::_minmax(const Db *db)
{
  double tt;
  if (db == nullptr) return;
  int nbands = getNDirs();

  if (db->isGrid())
  {
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(db);
    int nx = (dbgrid->getNDim() >= 1) ? dbgrid->getNX(0) : 1;
    int ny = (dbgrid->getNDim() >= 2) ? dbgrid->getNX(1) : 1;
    int nz = (dbgrid->getNDim() >= 3) ? dbgrid->getNX(2) : 1;

    /* Case when the data obeys to a grid organization */
    /* This test is programmed for 3-D (maximum) grid  */
    /* as the Turning Bands method is limited to 3-D   */

    for (int ibs = 0; ibs < nbands; ibs++)
    {
      for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 2; iy++)
          for (int ix = 0; ix < 2; ix++)
          {
            tt = _codirs[ibs].projectGrid(dbgrid, ix * (nx - 1), iy * (ny - 1),
                                          iz * (nz - 1));
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

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      for (int ibs = 0; ibs < nbands; ibs++)
      {
        if (!db->isActive(iech)) continue;
        tt = _codirs[ibs].projectPoint(db, iech);
        if (tt < _getCodirTmin(ibs)) _setCodirTmin(ibs, tt);
        if (tt > _getCodirTmax(ibs)) _setCodirTmax(ibs, tt);
        double delta = _getCodirTmax(ibs) - _getCodirTmin(ibs);
        if (_field < delta) _field = delta;
      }
    }
  }

  _npointSimulated += db->getSampleNumber();
  return;
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
  int nmini = 5;
  int nmaxi = 5000;

  int naverage = _npointSimulated / _nbtuba;
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
 **  Particular case of the K-Bessel model. It must be turned into:
 **  - Exponential : when param is too close to 0.5
 **
 ** \return  The modified type
 **
 *****************************************************************************/
ECov CalcSimuTurningBands::_particularCase(const ECov &type, double param)
{
  static double eps = 1.e-7;

  switch (type.toEnum())
  {
    case ECov::E_STABLE:
      if (ABS(param - 1.) < eps) return (ECov::EXPONENTIAL);
      if (ABS(param - 2.) < eps) return (ECov::GAUSSIAN);
      return (ECov::STABLE);
      break;

    case ECov::E_BESSEL_K:
      if (ABS(param - 0.5) < eps) return (ECov::EXPONENTIAL);
      break;

    default:
      break;
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
int CalcSimuTurningBands::_initializeSeedBands()
{
  double tdeb, omega, phi, correc, correc0;
  VectorDouble v0;
  VectorDouble v1;
  VectorDouble v2;
  VectorDouble t;

  /* Initializations */

  _setDensity();
  int ncova  = _getNCova();
  int nvar   = _getNVar();
  int nbsimu = getNbSimu();
  double theta1 = 1. / _theta;

  /* Loop on the turning bands */

  int mem_seed = law_get_random_seed();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int ibs = 0;
    for (int isimu = ibs = 0; isimu < nbsimu; isimu++)
      for (int is = 0; is < ncova; is++)
        for (int ib = 0; ib < _nbtuba; ib++, ibs++)
        {
          double tmin  = _getCodirTmin(ibs);
          double tmax  = _getCodirTmax(ibs);
          double scale = _getCodirScale(ibs);
          ECov type    = getModel()->getCovaType(is);
          double param = getModel()->getParam(is);
          type = _particularCase(type, param);
          _setSeedBand(ivar, is, ib, isimu, law_get_random_seed());

          switch (type.toEnum())
          {
            case ECov::E_NUGGET:
              // Next line is simply to let the random number cycle
              (void) law_gaussian();
              break;

            case ECov::E_EXPONENTIAL:
              scale *= 2.;
              t = _migration(tmin, tmax, scale);
              (void) law_uniform(0., 1.);
              break;

            case ECov::E_SPHERICAL:
              t = _dilution(tmin, tmax, scale, &tdeb);
              break;

            case ECov::E_CUBIC:
              t = _dilution(tmin, tmax, scale, &tdeb);
              break;

            case ECov::E_GAUSSIAN:
            case ECov::E_SINCARD:
              _spectral(type, scale, 2., &omega, &phi);
              break;

            case ECov::E_BESSEL_J:
              _spectral(type, scale, param, &omega, &phi);
              break;

            case ECov::E_BESSEL_K:
              if (param > 0.5)
                _spectral(type, scale, param, &omega, &phi);
              else
              {
                scale = _computeScaleKB(param, scale) * 2;
                t = _migration(tmin, tmax, scale);
                (void) law_uniform(0., 1.);
              }
              break;

            case ECov::E_STABLE:
              if (param > 1)
                _spectral(type, scale, param, &omega, &phi);
              else
              {
                scale *= 2;
                scale = _computeScale(param, scale);
                t = _migration(tmin, tmax, scale);
                (void) law_uniform(0., 1.);
              }
              break;

            case ECov::E_POWER:
              _power1D(ib, scale, param, &omega, &phi, &correc, &correc0);
              break;

            case ECov::E_SPLINE_GC:
              _spline1D(ib, scale, 1, &omega, &phi, &correc, &correc0);
              break;

            case ECov::E_LINEAR:
            case ECov::E_ORDER1_GC:
            case ECov::E_ORDER3_GC:
            case ECov::E_ORDER5_GC:
              t = _migration(tmin, tmax, theta1);
              _irfProcess(type, t, v0, v1, v2);
              break;

            default:
              messerr("The structure (%s) cannot be simulated",
                      type.getDescr().c_str());
              messerr("using the Turning Bands algorithm");
              return 1;
          }
        }
  }

  // Set the initial seed back
  law_set_random_seed(mem_seed);
  return 0;
}

/*****************************************************************************/
/*!
 **  Generate a migration process
 **
 ** \param[in]  tmin  minimum value
 ** \param[in]  tmax  maximum value
 ** \param[in]  scale scale of the exponential
 ** \param[in]  eps   Epsilon value
 **
 *****************************************************************************/
VectorDouble CalcSimuTurningBands::_migration(double tmin,
                                          double tmax,
                                          double scale,
                                          double eps)
{
  VectorDouble tab;

  /* Initializations */

  double delta = tmax - tmin;
  if (scale < delta * eps)
  {
    int count = (int) ceil(delta / eps);
    for (int i = 0; i < count; i++)
      tab.push_back(law_gaussian());
  }
  else
  {
    double value;
    value = tmin + scale * log(law_uniform(0., 1.));
    tab.push_back(value);
    value = tmin - scale * log(law_uniform(0., 1.));
    tab.push_back(value);

    while (value <= tmax)
    {
      value -= scale * log(law_uniform(0., 1.));
      tab.push_back(value);
    }
  }
  return tab;
}

/*****************************************************************************/
/*!
 **  Generate a dilution process
 **
 ** \param[in]  tmin minimum value
 ** \param[in]  tmax maximum value
 ** \param[in]  mesh mesh of the random walk
 **
 ** \param[out] start  initial time
 **
 *****************************************************************************/
VectorDouble CalcSimuTurningBands::_dilution(double tmin,
                                         double tmax,
                                         double mesh,
                                         double *start)
{
  VectorDouble tab;
  int count = 0;
  double tdeb = tmin - mesh * law_uniform(0., 1.);

  while (tdeb + count * mesh <= tmax)
  {
    double value = (law_uniform(0., 1.) < 0.5) ? -1. : 1.;
    tab.push_back(value);
    count++;
  }
  *start = tdeb;
  return tab;
}

/*****************************************************************************/
/*!
 **  Prepare a spectral method
 **
 ** \param[in]  type   type of method for generating the period
 ** \param[in]  scale  scale factor
 ** \param[in]  param  factor for the period
 **
 ** \param[out] omega  resulting period
 ** \param[out] phi    uniform phase lying within [0,2 PI]
 **
 *****************************************************************************/
void CalcSimuTurningBands::_spectral(const ECov &type,
                                 double scale,
                                 double param,
                                 double *omega,
                                 double *phi)
{
  double val, period;
  int i;

  period = 0.;
  switch (type.toEnum())
  {
    case ECov::E_GAUSSIAN:
      for (i = 0; i < 3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      period = sqrt(param * period) / scale;
      break;

    case ECov::E_STABLE:
      for (i = 0; i < 3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      scale = _computeScale(param, scale);
      period = sqrt(2. * period) / scale;
      break;

    case ECov::E_SINCARD:
      val = (law_uniform(0., 1.) >= 0.5) ? 1 : -1;
      period = val / scale;
      break;

    case ECov::E_BESSEL_J:
      val = law_beta1(1.5, param - 0.5);
      period = sqrt(val) / scale;
      break;

    case ECov::E_BESSEL_K:
      param = sqrt(2. * law_gamma(param));
      for (i = 0; i < 3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      period = sqrt(period) / (param * scale);
      break;

    default:
      break;
  }

  *omega = period;
  *phi = 2. * GV_PI * law_uniform(0., 1.);

  return;
}

/****************************************************************************/
/*!
 **  Calculate the scale for 1D process for the stable model
 **
 ** \return  Scale parameter of the 1D process to simulate
 **
 ** \param[in]  alpha       Third parameter of the stable covariance model
 ** \param[in]  scale       Scale parameter of the model
 **
 *****************************************************************************/
double CalcSimuTurningBands::_computeScale(double alpha, double scale)

{
  double scale_out;

  if (alpha < 1)
    scale_out = scale / law_stable_standard_abgd(alpha);
  else
    scale_out = scale / sqrt(law_stable_standard_abgd(alpha / 2.));

  return (scale_out);
}

/****************************************************************************/
/*!
 **  Calculate the scale for 1D process for the K-Bessel model (param<0.5)
 **
 ** \return  Scale parameter of the 1D process to simulate (param<0.5)
 **
 ** \param[in]  param       Third parameter of the K-Bessel covariance model
 ** \param[in]  scale       Scale parameter of the model
 **
 *****************************************************************************/
double CalcSimuTurningBands::_computeScaleKB(double param, double scale)

{
  return (scale * sqrt(law_beta1(param, 0.5 - param)));
}

/*****************************************************************************/
/*!
 **  Generate the 1D stochastic process
 **
 ** \param[in]  ib      Rank of the turning band
 ** \param[in]  scale   scale factor
 ** \param[in]  alpha   power of the variogram h^alpha
 **
 ** \param[out] omega   period = 2piR
 **                     where R is a second kind beta variable with parameters
 **                     1-alpha/2 and alpha/2
 ** \param[out] phi     uniform phase lying within [0,2 PI]
 ** \param[out] theta_3 value of theta_alpha,3(R)
 ** \param[out] correc0 value to substract from Y_alpha,2 in order to avoid
 **                     numerical problems when R is very small
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
void CalcSimuTurningBands::_power1D(int ib,
                                double scale,
                                double alpha,
                                double *omega,
                                double *phi,
                                double *theta_3,
                                double *correc0)
{
  double R, theta_1;
  static double twoPI, log3s2, log1s2, logap1, logap1s2, logap3s2, as2, coeff,
      coeff3;
  static double alpha_mem = -1.;

  if (ib == 0 || alpha != alpha_mem)
  {
    as2 = alpha / 2.;
    twoPI = 2. * GV_PI;
    log3s2 = loggamma(1.5);
    log1s2 = loggamma(0.5);
    logap1 = loggamma(1.0 + alpha);
    logap1s2 = loggamma(0.5 + as2);
    logap3s2 = loggamma(1.5 + as2);
    coeff = 2. * sqrt(exp(logap1) / pow(twoPI, alpha)) * pow(scale, as2);
    coeff3 = sqrt(exp(log3s2 + logap1s2 - log1s2 - logap3s2));
    alpha_mem = alpha;
  }

  R = law_beta2(1 - as2, as2);
  theta_1 = coeff * sqrt((R + 1.) / pow(R, as2 + 1));
  *correc0 = cos(*phi);
  *omega = twoPI * R;
  *phi = twoPI * law_uniform(0., 1.);
  *theta_3 = theta_1 / coeff3;
}

/*****************************************************************************/
/*!
 **  Generate the 1D stochastic process for spline covariance
 **
 ** \param[in]  ib      Rank of the turning band
 ** \param[in]  scale   scale factor
 ** \param[in]  k       power of the variogram h^(2k)log(h)
 **
 ** \param[out] omega   period = 2piR
 **                     where R is a second kind beta variable with parameters
 **                     1/2 and 1/2
 ** \param[out] phi     uniform phase lying within [0,2 PI]
 ** \param[out] xi_3    value of xi_2k,3(R)
 ** \param[out] correc0 value to substract from S_2k,3 in order to avoid
 **                     numerical problems when R is very small
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
void CalcSimuTurningBands::_spline1D(int ib,
                                 double scale,
                                 int k,
                                 double *omega,
                                 double *phi,
                                 double *xi_3,
                                 double *correc0)
{
  double R;
  int twokm1;
  static double twoPI, twokp1s2, log3s2, logkp3s2, logkp1, coeff;
  static int twok = 0;
  static int k_mem = -1;

  if (ib == 0 || k != k_mem)
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
  *xi_3 = coeff * sqrt((R + 1.) / pow(R, twokp1s2));
  *correc0 = cos(*phi);
  *omega = twoPI * R * pow(scale, twok);
  *phi = twoPI * law_uniform(0., 1.);
}

/*****************************************************************************/
/*!
 **  Generates the process constituted by independent gaussian
 **  variables along a 1D Poisson process. The process consists in
 **  the integration(s) of the previous process Perform the core
 **  allocation
 **
 ** \param[in]  type   degree of the IRF -> number of integrations
 ** \param[in]  t      Vector of Poisson delimitors
 **
 ** \param[out] v0     Wiener-Levy process
 ** \param[out] v1     First integration of the Wiener-Levy process
 ** \param[out] v2     Second integration of the Wiener-Levy process
 **
 ** \remark  This procedure allocates memory that should be freed
 **
 *****************************************************************************/
void CalcSimuTurningBands::_irfProcess(const ECov &type,
                                   const VectorDouble &t,
                                   VectorDouble &v0,
                                   VectorDouble &v1,
                                   VectorDouble &v2)
{
  double delta;

  int level = -1;
  if (type == ECov::LINEAR || type == ECov::ORDER1_GC) level = 0;
  if (type == ECov::ORDER3_GC) level = 1;
  if (type == ECov::ORDER5_GC) level = 2;

  int nt = (int) t.size();

  v0.clear();
  v1.clear();
  v2.clear();

  /* Generation of the Wiener-Levy process and its integrations */

  double val0 = 0.;
  double val1 = 0.;
  double val2 = 0.;
  v0.push_back(val0);
  v1.push_back(val1);
  v2.push_back(val2);
  for (int i = 1; i < nt; i++)
  {
    val0 += law_gaussian();
    v0.push_back(val0);
    if (level == 0) continue;
    delta = t[i] - t[i - 1];
    val1 += val0 * delta;
    if (level == 1) continue;
    val2 += val1 * delta + val0 * delta * delta / 2.;
    v2.push_back(val2);
  }
}

/*****************************************************************************/
/*!
 **  Calculate the linear model of coregionalization starting from the
 **  coregionalization matrix
 **
 ** \return  The Vector containing the AIC matrix
 **
 ** \remark  In case of error, the message is printed by the routine
 ** \remark  Warning: in the case of linked drift, the test of definite
 ** \remark  positiveness is bypassed as we are not in the scope of the
 ** \remark  linear model of coregionalization anymore.
 ** \remark  As a consequence the array "aic()" is not evaluated
 **
 *****************************************************************************/
VectorDouble CalcSimuTurningBands::_createAIC()
{
  int ncova = _getNCova();
  int nvar  = _getNVar();

  VectorDouble valpro(nvar);
  VectorDouble vecpro(nvar * nvar);
  VectorDouble aic(ncova * nvar * nvar);

  /* Calculate the eigen values and vectors of the coregionalization matrix */

  for (int icov = 0; icov < ncova; icov++)
  {
    if (!is_matrix_definite_positive(
        nvar, getModel()->getCova(icov)->getSill().getValues().data(),
        valpro.data(), vecpro.data(), 0))
    {
      messerr("Warning: the model is not authorized");
      messerr("The coregionalization matrix for the structure %d is not definite positive",
          icov + 1);
      return VectorDouble();
    }

    /* Calculate the factor matrix */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        int ijvar = ivar + nvar * jvar;
        aic[icov*nvar*nvar + ijvar] = vecpro[ijvar] * sqrt(valpro[ivar]);
      }
  }
  return aic;
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a set of points using
 **  Turning Bands method.
 **
 ** \param[in]  db         Db structure
 ** \param[in]  aic        Array 'aic'
 ** \param[in]  icase      Rank of PGS or GRF
 ** \param[in]  shift      Shift before writing the simulation result
 **
 *****************************************************************************/
void CalcSimuTurningBands::_simulatePoint(Db *db,
                                      const VectorDouble &aic,
                                      int icase,
                                      int shift)
{
  double vexp, tdeb, dt0, t0;
  int nt0;
  VectorDouble t;
  VectorDouble v0;
  VectorDouble v1;
  VectorDouble v2;
  VectorDouble tab;
  static double vexp1 = 0.1;
  static double vexp2 = 0.1967708298;
  double phi = 0.;
  double omega = 0.;

  /* Initializations */

  int nech = db->getSampleNumber();
  int ncova = _getNCova();
  int nvar = _getNVar();
  int nbsimu = getNbSimu();
  double theta1 = 1. / _theta;
  double norme = sqrt(1. / _nbtuba);

  /* Core allocation */

  tab.resize(nech,0.);

  /*****************************/
  /* Performing the simulation */
  /*****************************/

  int mem_seed = law_get_random_seed();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int ibs = 0;
    for (int isimu = 0; isimu < nbsimu; isimu++)
      for (int is = 0; is < ncova; is++)
        for (int ib = 0; ib < _nbtuba; ib++, ibs++)
        {
          double tmin  = _getCodirTmin(ibs);
          double tmax  = _getCodirTmax(ibs);
          double scale = _getCodirScale(ibs);
          double param = getModel()->getParam(is);
          ECov type    = getModel()->getCovaType(is);
          double correc = 1.;
          double correc0 = 0.;
          type = _particularCase(type, param);
          law_set_random_seed(_getSeedBand(ivar, is, ib, isimu));

          switch (type.toEnum())
          {
            case ECov::E_NUGGET:
              break;

            case ECov::E_STABLE:
              if (param > 1)
              {
                correc = sqrt(2.);
                _spectral(type, scale, param, &omega, &phi);
                for (int iech = 0; iech < nech; iech++)
                {
                  if (! db->isActive(iech)) continue;
                  t0 = _codirs[ibs].projectPoint(db, iech);
                  tab[iech] = cos(omega * t0 + phi);
                }
              }
              else
              {
                scale *= 2;
                scale = _computeScale(param, scale);
                t = _migration(tmin, tmax, scale);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);
                nt0 = 0;
                for (int iech = 0; iech < nech; iech++)
                {
                  if (!db->isActive(iech)) continue;
                  t0 = _codirs[ibs].projectPoint(db, iech);
                  nt0 = _rankInPoisson(nt0, t0, t);
                  tab[iech] = (2. * t0 > t[nt0 + 1] + t[nt0]) ? -vexp : vexp;
                }
              }
              break;

            case ECov::E_EXPONENTIAL:
              scale *= 2.;
              t = _migration(tmin, tmax, scale);
              vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);
              nt0 = 0;
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                nt0 = _rankInPoisson(nt0, t0, t);
                tab[iech] = (2. * t0 > t[nt0 + 1] + t[nt0]) ? -vexp : vexp;
              }
              break;

            case ECov::E_SPHERICAL:
              correc = sqrt(3.);
              t = _dilution(tmin, tmax, scale, &tdeb);
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                nt0 = _rankRegular(t0, tdeb, scale);
                dt0 = (t0 - tdeb) - scale * nt0;
                double r = dt0 / scale;
                tab[iech] = t[nt0] * (2. * r - 1.);
              }
              break;

            case ECov::E_CUBIC:
              correc = sqrt(840.);
              t = _dilution(tmin, tmax, scale, &tdeb);
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                nt0 = _rankRegular(t0, tdeb, scale);
                dt0 = (t0 - tdeb) - scale * nt0;
                double r = dt0 / scale;
                tab[iech] = t[nt0] * r * (r - 0.5) * (r - 1.);
              }
              break;

            case ECov::E_POWER:
              _power1D(ib, scale, param, &omega, &phi, &correc, &correc0);
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                tab[iech] = cos(omega * t0 + phi) - correc0;
              }
              break;

            case ECov::E_SPLINE_GC:
              _spline1D(ib, scale, 1, &omega, &phi, &correc, &correc0);
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                tab[iech] = cos(omega * t0 + phi) - correc0;
              }
              break;

            case ECov::E_GAUSSIAN:
            case ECov::E_SINCARD:
            case ECov::E_BESSEL_J:
            case ECov::E_BESSEL_K:
              correc = sqrt(2.);
              switch (type.toEnum())
              {
                case ECov::E_GAUSSIAN:
                case ECov::E_SINCARD:
                  _spectral(type, scale, 2., &omega, &phi);
                  break;

                case ECov::E_BESSEL_J:
                case ECov::E_BESSEL_K:
                  _spectral(type, scale, param, &omega, &phi);
                  break;
                default:
                  break;
              }

              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                tab[iech] = cos(omega * t0 + phi);
              }
              break;

            case ECov::E_LINEAR:
            case ECov::E_ORDER1_GC:
            case ECov::E_ORDER3_GC:
            case ECov::E_ORDER5_GC:
              t = _migration(tmin, tmax, theta1);
              _irfProcess(type, t, v0, v1, v2);
              correc = _irfCorrec(type, theta1, scale);
              nt0 = 0;
              for (int iech = 0; iech < nech; iech++)
              {
                if (!db->isActive(iech)) continue;
                t0 = _codirs[ibs].projectPoint(db, iech);
                nt0 = _rankInPoisson(nt0, t0, t);
                tab[iech] = _irfProcessSample(type, nt0, t0, t, v0, v1, v2);
              }
              break;

            default:
              break;
          }

          if (type != ECov::NUGGET)
            for (int iech = 0; iech < nech; iech++)
              if (db->isActive(iech))
                for (int jvar = 0; jvar < nvar; jvar++)
                  db->updSimvar(ELoc::SIMU, iech, shift + isimu, jvar, icase,
                                nbsimu, nvar, 0,
                                tab[iech] * correc * _getAIC(aic, is, jvar, ivar));
        }
  }

  /* Normation */

  for (int isimu = 0; isimu < nbsimu; isimu++)
    for (int iech = 0; iech < nech; iech++)
      for (int jvar = 0; jvar < nvar; jvar++)
        if (db->isActive(iech))
          db->updSimvar(ELoc::SIMU, iech, shift + isimu, jvar, icase, nbsimu,
                        nvar, 1, norme);

  // Set the initial seed back
  law_set_random_seed(mem_seed);
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a set of gradient points using
 **  Turning Bands method.
 **
 ** \param[in]  dbgrd      Gradient Db structure
 ** \param[in]  aic        Array 'aic'
 ** \param[in]  delta      Value of the increment
 **
 ** \remarks The simulated gradients are stored as follows:
 ** \remarks idim * nbsimu + isimu (for simulation at first point)
 ** \remarks idim * nbsimu + isimu + ndim * nbsimu (for simulation at 2nd point)
 ** \remarks At the end, the simulated gradient is stored at first point
 **
 *****************************************************************************/
void CalcSimuTurningBands::_simulateGradient(Db *dbgrd,
                                         const VectorDouble &aic,
                                         double delta)
{
  int jsimu;
  int icase = 0;
  int ndim = dbgrd->getNDim();
  int nbsimu = getNbSimu();

  for (int idim = 0; idim < ndim; idim++)
  {

    /* Simulation at the initial location */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu;
      _simulatePoint(dbgrd, aic, icase, jsimu);
    }

    /* Shift the information */

    for (int iech = 0; iech < dbgrd->getSampleNumber(); iech++)
      if (dbgrd->isActive(iech))
        dbgrd->setCoordinate(iech, idim,
                             dbgrd->getCoordinate(iech, idim) + delta);

    /* Simulation at the shift location */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu + ndim * nbsimu;
      _simulatePoint(dbgrd, aic, icase, jsimu);
    }

    /* Un-Shift the information */

    for (int iech = 0; iech < dbgrd->getSampleNumber(); iech++)
      if (dbgrd->isActive(iech))
        dbgrd->setCoordinate(iech, idim,
                             dbgrd->getCoordinate(iech, idim) - delta);

    /* Scaling */

    for (int isimu = 0; isimu < nbsimu; isimu++)
      for (int iech = 0; iech < dbgrd->getSampleNumber(); iech++)
      {
        if (!dbgrd->isActive(iech)) continue;
        jsimu = isimu + idim * nbsimu + ndim * nbsimu;
        double value2 = dbgrd->getSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                                         2 * ndim * nbsimu, 1);
        jsimu = isimu + idim * nbsimu;
        double value1 = dbgrd->getSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                                         2 * ndim * nbsimu, 1);
        dbgrd->setSimvar(ELoc::SIMU, iech, jsimu, 0, icase,
                         2 * ndim * nbsimu,
                         1, (value2 - value1) / delta);
      }
  }
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a set of tangent points using
 **  Turning Bands method.
 **
 ** \param[in]  dbtgt      Tangent Db structure
 ** \param[in]  aic        Array 'aic'
 ** \param[in]  delta      Value of the increment
 **
 ** \remarks Warning: To perform the simulation of the tangent, we must
 ** \remarks simulated the gradients first. So we need to dimension the
 ** \remarks simulation outcome variables as for the gradients
 **
 *****************************************************************************/
void CalcSimuTurningBands::_simulateTangent(Db *dbtgt,
                                            const VectorDouble &aic,
                                            double delta)
{
  int icase = 0;
  int nvar = _getNVar();
  int nbsimu = getNbSimu();

  /* Perform the simulation of the gradients at tangent points */

  _simulateGradient(dbtgt, aic, delta);

  /* Calculate the simulated tangent */

  for (int isimu = 0; isimu < nbsimu; isimu++)
    for (int iech = 0; iech < dbtgt->getSampleNumber(); iech++)
    {
      if (!dbtgt->isActive(iech)) continue;

      double value = 0.;
      for (int idim = 0; idim < dbtgt->getNDim(); idim++)
        value += dbtgt->getLocVariable(ELoc::TGTE,iech, idim)
            * dbtgt->getSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu,
                               nvar);
      dbtgt->setSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu,
                       nvar, value);
    }
}

/*****************************************************************************/
/*!
 **  Perform non-conditional simulations on a grid using the
 **  Turning Bands method
 **
 ** \param[in]  db         Db structure
 ** \param[in]  aic        Array 'aic'
 ** \param[in]  icase      Rank of PGS or GRF
 ** \param[in]  shift      Shift before writing the simulation result
 **
 *****************************************************************************/
void CalcSimuTurningBands::_simulateGrid(DbGrid *db,
                                     const VectorDouble &aic,
                                     int icase,
                                     int shift)
{
  double vexp, phi, tdeb, omega, dt0, dt, t0, t0y, t0z;
  double cxp, sxp, cyp, syp, czp, szp, c0x, s0x, c0y, s0y, c0z, s0z, c1, s1;
  int nt0, ind;
  VectorDouble t;
  VectorDouble v0;
  VectorDouble v1;
  VectorDouble v2;
  VectorDouble tab;
  static double vexp1 = 0.1;
  static double vexp2 = 0.1967708298;

  /* Initializations */

  int nbsimu = getNbSimu();
  double theta1 = 1. / _theta;
  int nvar   = _getNVar();
  int ncova  = _getNCova();
  int nx     = (db->getNDim() >= 1) ? db->getNX(0) : 1;
  int ny     = (db->getNDim() >= 2) ? db->getNX(1) : 1;
  int nz     = (db->getNDim() >= 3) ? db->getNX(2) : 1;
  int nech   = nx * ny * nz;
  double norme  = sqrt(1. / _nbtuba);

  /* Core allocation */

  tab.resize(nech, 0.);
  omega = 0.;
  phi = 0.;

  /*****************************/
  /* Performing the simulation */
  /*****************************/

  int mem_seed = law_get_random_seed();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int ibs = 0;
    for (int isimu = 0; isimu < nbsimu; isimu++)
      for (int is = 0; is < ncova; is++)
        for (int ib = 0; ib < _nbtuba; ib++, ibs++)
        {
          double tmin = _getCodirTmin(ibs);
          double tmax = _getCodirTmax(ibs);
          double dxp = _getCodirDXP(ibs);
          double dyp = _getCodirDYP(ibs);
          double dzp = _getCodirDZP(ibs);
          double t00 = _getCodirT00(ibs);
          double scale = _getCodirScale(ibs);
          ECov type = getModel()->getCovaType(is);
          double param = getModel()->getParam(is);
          double correc = 1.;
          double correc0 = 0.;
          type = _particularCase(type, param);
          law_set_random_seed(_getSeedBand(ivar, is, ib, isimu));

          switch (type.toEnum())
          {
            case ECov::E_NUGGET:
              break;

            case ECov::E_STABLE:
              if (param > 1)
              {
                correc = sqrt(2.);
                _spectral(type, scale, param, &omega, &phi);
                _getOmegaPhi(ibs, omega, phi,
                             &cxp, &sxp, &cyp, &syp, &czp, &szp, &c0z, &s0z);
                ind = 0;
                for (int iz = 0; iz < nz; iz++)
                {
                  c0y = c0z;
                  s0y = s0z;
                  c1 = c0z * czp - s0z * szp;
                  s1 = s0z * czp + c0z * szp;
                  c0z = c1;
                  s0z = s1;
                  for (int iy = 0; iy < ny; iy++)
                  {
                    c0x = c0y;
                    s0x = s0y;
                    c1 = c0y * cyp - s0y * syp;
                    s1 = s0y * cyp + c0y * syp;
                    c0y = c1;
                    s0y = s1;
                    for (int ix = 0; ix < nx; ix++, ind++)
                    {
                      if (db->isActive(ind)) tab[ind] = c0x - correc0;
                      c1 = c0x * cxp - s0x * sxp;
                      s1 = s0x * cxp + c0x * sxp;
                      c0x = c1;
                      s0x = s1;
                    }
                  }
                }
              }
              else
              {
                scale *= 2;
                scale = _computeScale(param, scale);
                t = _migration(tmin, tmax, scale);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);
                t0z = t00;

                nt0 = 0;
                ind = 0;
                for (int iz = 0; iz < nz; iz++)
                {
                  t0y = t0z;
                  t0z += dzp;
                  for (int iy = 0; iy < ny; iy++)
                  {
                    t0 = t0y;
                    t0y += dyp;
                    for (int ix = 0; ix < nx; ix++, ind++)
                    {
                      if (db->isActive(ind))
                      {
                        nt0 = _rankInPoisson(nt0, t0, t);
                        tab[ind] = (2. * t0 > t[nt0 + 1] + t[nt0]) ? -vexp : vexp;
                      }
                      t0 += dxp;
                    }
                  }
                }
              }
              break;

            case ECov::E_BESSEL_K:
              if (param > 0.5)
              {
                correc = sqrt(2.);
                _spectral(type, scale, param, &omega, &phi);
                _getOmegaPhi(ibs, omega, phi,
                             &cxp, &sxp, &cyp, &syp, &czp, &szp, &c0z, &s0z);
                ind = 0;
                for (int iz = 0; iz < nz; iz++)
                {
                  c0y = c0z;
                  s0y = s0z;
                  c1 = c0z * czp - s0z * szp;
                  s1 = s0z * czp + c0z * szp;
                  c0z = c1;
                  s0z = s1;
                  for (int iy = 0; iy < ny; iy++)
                  {
                    c0x = c0y;
                    s0x = s0y;
                    c1 = c0y * cyp - s0y * syp;
                    s1 = s0y * cyp + c0y * syp;
                    c0y = c1;
                    s0y = s1;
                    for (int ix = 0; ix < nx; ix++, ind++)
                    {
                      if (db->isActive(ind)) tab[ind] = c0x - correc0;
                      c1 = c0x * cxp - s0x * sxp;
                      s1 = s0x * cxp + c0x * sxp;
                      c0x = c1;
                      s0x = s1;
                    }
                  }
                }
              }
              else
              {
                scale = _computeScaleKB(param, scale) * 2;
                t = _migration(tmin, tmax, scale);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);
                t0z = t00;

                nt0 = 0;
                ind = 0;
                for (int iz = 0; iz < nz; iz++)
                {
                  t0y = t0z;
                  t0z += dzp;
                  for (int iy = 0; iy < ny; iy++)
                  {
                    t0 = t0y;
                    t0y += dyp;
                    for (int ix = 0; ix < nx; ix++, ind++)
                    {
                      if (db->isActive(ind))
                      {
                        nt0 = _rankInPoisson(nt0, t0, t);
                        tab[ind] = (2. * t0 > t[nt0 + 1] + t[nt0]) ? -vexp : vexp;
                      }
                      t0 += dxp;
                    }
                  }
                }
              }
              break;

            case ECov::E_EXPONENTIAL:
              scale *= 2.;
              t = _migration(tmin, tmax, scale);
              vexp = 1. - vexp1 + vexp2 * law_uniform(0., 1.);

              t0z = t00;
              nt0 = 0;
              ind = 0;
              for (int iz = 0; iz < nz; iz++)
              {
                t0y = t0z;
                t0z += dzp;
                for (int iy = 0; iy < ny; iy++)
                {
                  t0 = t0y;
                  t0y += dyp;
                  for (int ix = 0; ix < nx; ix++, ind++)
                  {
                    if (db->isActive(ind))
                    {
                      nt0 = _rankInPoisson(nt0, t0, t);
                      tab[ind] = (2. * t0 > t[nt0 + 1] + t[nt0]) ? -vexp : vexp;
                    }
                    t0 += dxp;
                  }
                }
              }
              break;

            case ECov::E_SPHERICAL:
              correc = sqrt(3.);
              t = _dilution(tmin, tmax, scale, &tdeb);
              tdeb /= scale;
              t0z = t00;
              ind = 0;
              for (int iz = 0; iz < nz; iz++)
              {
                t0y = t0z;
                t0z += dzp;
                for (int iy = 0; iy < ny; iy++)
                {
                  t0 = t0y;
                  t0y += dyp;
                  for (int ix = 0; ix < nx; ix++, ind++)
                  {
                    if (db->isActive(ind))
                    {
                      dt = (t0 - tdeb);
                      nt0 = (int) dt;
                      dt0 = dt - nt0;
                      tab[ind] = t[nt0] * (2. * dt0 - 1.);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;

            case ECov::E_CUBIC:
              correc = sqrt(840.);
              t = _dilution(tmin, tmax, scale, &tdeb);
              tdeb /= scale;
              t0z = t00;
              ind = 0;
              for (int iz = 0; iz < nz; iz++)
              {
                t0y = t0z;
                t0z += dzp;
                for (int iy = 0; iy < ny; iy++)
                {
                  t0 = t0y;
                  t0y += dyp;
                  for (int ix = 0; ix < nx; ix++, ind++)
                  {
                    if (db->isActive(ind))
                    {
                      dt = (t0 - tdeb);
                      nt0 = (int) dt;
                      dt0 = dt - nt0;
                      tab[ind] = t[nt0] * dt0 * (dt0 - 0.5) * (dt0 - 1.);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;

            case ECov::E_GAUSSIAN:
            case ECov::E_SINCARD:
            case ECov::E_POWER:
            case ECov::E_SPLINE_GC:
            case ECov::E_BESSEL_J:
              correc = sqrt(2.);
              switch (type.toEnum())
              {
                case ECov::E_POWER:
                  _power1D(ib, scale, param, &omega, &phi, &correc, &correc0);
                  break;

                case ECov::E_SPLINE_GC:
                  _spline1D(ib, scale, 1, &omega, &phi, &correc, &correc0);
                  break;

                case ECov::E_GAUSSIAN:
                case ECov::E_SINCARD:
                  _spectral(type, scale, 2., &omega, &phi);
                  break;

                case ECov::E_BESSEL_J:
                  _spectral(type, scale, param, &omega, &phi);
                  break;

                default:
                  break;
              }

              _getOmegaPhi(ibs, omega, phi,
                           &cxp, &sxp, &cyp, &syp, &czp, &szp, &c0z, &s0z);
              ind = 0;
              for (int iz = 0; iz < nz; iz++)
              {
                c0y = c0z;
                s0y = s0z;
                c1 = c0z * czp - s0z * szp;
                s1 = s0z * czp + c0z * szp;
                c0z = c1;
                s0z = s1;
                for (int iy = 0; iy < ny; iy++)
                {
                  c0x = c0y;
                  s0x = s0y;
                  c1 = c0y * cyp - s0y * syp;
                  s1 = s0y * cyp + c0y * syp;
                  c0y = c1;
                  s0y = s1;
                  for (int ix = 0; ix < nx; ix++, ind++)
                  {
                    if (db->isActive(ind)) tab[ind] = c0x - correc0;
                    c1 = c0x * cxp - s0x * sxp;
                    s1 = s0x * cxp + c0x * sxp;
                    c0x = c1;
                    s0x = s1;
                  }
                }
              }
              break;

            case ECov::E_LINEAR:
            case ECov::E_ORDER1_GC:
            case ECov::E_ORDER3_GC:
            case ECov::E_ORDER5_GC:
              t = _migration(tmin, tmax, theta1);
              _irfProcess(type, t, v0, v1, v2);
              correc = _irfCorrec(type, theta1, scale);

              t0z = t00;
              nt0 = 0;
              ind = 0;
              for (int iz = 0; iz < nz; iz++)
              {
                t0y = t0z;
                t0z += dzp;
                for (int iy = 0; iy < ny; iy++)
                {
                  t0 = t0y;
                  t0y += dyp;
                  for (int ix = 0; ix < nx; ix++, ind++)
                  {
                    if (db->isActive(ind))
                    {
                      nt0 = _rankInPoisson(nt0, t0, t);
                      tab[ind] = _irfProcessSample(type, nt0, t0, t, v0, v1, v2);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;

            default:
              break;
          }

          if (type != ECov::NUGGET)
            for (int iech = 0; iech < nech; iech++)
              if (db->isActive(iech))
                for (int jvar = 0; jvar < nvar; jvar++)
                  db->updSimvar(ELoc::SIMU, iech, shift + isimu, jvar, icase,
                                nbsimu, nvar, 0,
                                tab[iech] * correc * _getAIC(aic, is, jvar, ivar));
        }
  }

  /* Normation */

  for (int isimu = 0; isimu < nbsimu; isimu++)
    for (int iech = 0; iech < nech; iech++)
      for (int jvar = 0; jvar < nvar; jvar++)
        if (db->isActive(iech))
          db->updSimvar(ELoc::SIMU, iech, shift + isimu, jvar, icase, nbsimu,
                        nvar, 1, norme);

  // Set the initial seed back
  law_set_random_seed(mem_seed);
}

/*****************************************************************************/
/*!
 **  Returns the rank of the point t0 in the Poisson point process
 **
 ** \param[in]  def_rank Rank of the Poisson point
 ** \param[in]  t0 starting time
 ** \param[in]  t  Poisson point process
 **
 *****************************************************************************/
int CalcSimuTurningBands::_rankInPoisson(int def_rank, double t0, const VectorDouble& t)

{
  int it, itp, itn;

  /* First, try with the default interval then the next one and finally
   the previous one */

  int nt = (int) t.size();
  if (t0 >= t[def_rank] && t0 < t[def_rank + 1])
    return (def_rank);
  else if (def_rank < (nt - 2) && t0 >= t[def_rank + 1] && t0 < t[def_rank + 2])
    return (def_rank + 1);
  else if (def_rank > 0 && t0 >= t[def_rank - 1] && t0 < t[def_rank])
    return (def_rank - 1);

  /* The default value is not good ==> dichotomy */

  itp = 0;
  itn = nt - 1;
  while (itn - itp > 1)
  {
    it = (itn + itp) / 2;
    if (t0 >= t[it])
      itp = it;
    else
      itn = it;
  }
  return (itp);
}

/*****************************************************************************/
/*!
 **  Returns the rank of the point t0 in a regular pavement
 **
 ** \param[in]  t0    starting time
 ** \param[in]  tdeb  origin on the line
 ** \param[in]  scale scaling factor
 **
 *****************************************************************************/
int CalcSimuTurningBands::_rankRegular(double t0, double tdeb, double scale)

{
  return ((int) ((t0 - tdeb) / scale));
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
double CalcSimuTurningBands::_irfCorrec(const ECov &type, double theta1, double scale)
{
  double correc;

  /* Initializations */

  correc = 1.;

  /* Dispatch */

  switch (type.toEnum())
  {
    case ECov::E_LINEAR:
    case ECov::E_ORDER1_GC:
      correc = sqrt((4. * theta1) / scale);
      break;

    case ECov::E_ORDER3_GC:
      correc = sqrt((48. * theta1) / scale) / scale;
      break;

    case ECov::E_ORDER5_GC:
      correc = sqrt((1440. * theta1) / scale) / scale / scale;
      break;
    default:
      break;
  }

  return (correc);
}

/*****************************************************************************/
/*!
 **  Sample the Wiener-Levy (integrated) process
 **
 ** \param[in]  type   type of polynomial generalized covariance
 ** \param[in]  nt0    Rank of the Poisson point
 ** \param[in]  t0     starting time
 ** \param[in]  t      Poisson point process
 ** \param[in]  v0     Wiener-Levy process
 ** \param[in]  v1     First integration of the Wiener-Levy process
 ** \param[in]  v2     Second integration of the Wiener-Levy process
 **
 *****************************************************************************/
double CalcSimuTurningBands::_irfProcessSample(const ECov &type,
                                           int nt0,
                                           double t0,
                                           const VectorDouble &t,
                                           const VectorDouble &v0,
                                           const VectorDouble &v1,
                                           const VectorDouble &v2)
{
  double value, delta;

  /* Initializations */

  delta = t0 - t[nt0];

  /* Wiener-Levy process */

  value = v0[nt0];
  if (type == ECov::LINEAR || type == ECov::ORDER1_GC) return (value);

  /* First integration of the Wiener-Levy process */

  value = v1[nt0] + v0[nt0] * delta;
  if (type == ECov::ORDER3_GC) return (value);

  /* Second integration of the Wiener-Levy process */

  value = v2[nt0] + v1[nt0] * delta + v0[nt0] * delta * delta / 2.;
  if (type == ECov::ORDER5_GC) return (value);

  return (TEST);
}

void CalcSimuTurningBands::_getOmegaPhi(int ibs,
                                        double omega,
                                        double phi,
                                        double *cxp,
                                        double *sxp,
                                        double *cyp,
                                        double *syp,
                                        double *czp,
                                        double *szp,
                                        double *c0z,
                                        double *s0z)
{
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

/*****************************************************************************/
/*!
 **  Add the contribution of the nugget effect to the non-conditional
 **  simulations
 **
 ** \param[in]  db         Db structure
 ** \param[in]  aic        Array 'aic'
 ** \param[in]  icase      Rank of PGS or GRF
 **
 *****************************************************************************/
void CalcSimuTurningBands::_simulateNugget(Db *db, const VectorDouble& aic, int icase)
{
  int nech = db->getSampleNumber();
  int ncova = _getNCova();
  int nvar = _getNVar();
  int nbsimu = getNbSimu();

  /* Do nothing if there is no nugget effect in the model */

  bool flag_used = false;
  for (int is = 0; is < ncova && flag_used == 0; is++)
  {
    if (getModel()->getCovaType(is) == ECov::NUGGET) flag_used = true;
  }
  if (!flag_used) return;

  /* Performing the simulation */

  int mem_seed = law_get_random_seed();
  for (int isimu = 0; isimu < nbsimu; isimu++)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int is = 0; is < ncova; is++)
      {
        ECov type = getModel()->getCovaType(is);

        if (type != ECov::NUGGET) continue;
        law_set_random_seed(_getSeedBand(ivar, is, 0, isimu));

        for (int iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          double nugget = law_gaussian();
          for (int jvar = 0; jvar < nvar; jvar++)
            db->updSimvar(ELoc::SIMU, iech, isimu, jvar, icase, nbsimu, nvar, 0,
                          nugget * _getAIC(aic, is, jvar, ivar));
        }
      }

  // Set the initial seed back
  law_set_random_seed(mem_seed);
  return;
}

double CalcSimuTurningBands::_getAIC(const VectorDouble &aic,
                                 int icov,
                                 int ivar,
                                 int jvar)
{
  int nvar = _getNVar();
  return aic[jvar + nvar * (ivar + nvar * icov)];
}

/*****************************************************************************/
/*!
 **  Convert the non conditional simulations at the data points
 **  into simulation error
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  icase      Case for PGS or GRF
 ** \param[in]  flag_pgs   1 if called from PGS
 ** \param[in]  flag_gibbs 1 if called from Gibbs
 ** \param[in]  flag_dgm   1 if in the Discrete Gaussian Model
 **
 *****************************************************************************/
void CalcSimuTurningBands::_difference(Db *dbin,
                                       Model* model,
                                       int icase,
                                       bool flag_pgs,
                                       bool flag_gibbs,
                                       bool flag_dgm)
{
  int nbsimu = getNbSimu();
  int nvar = _getNVar();
  double r = 1.;
  if (flag_dgm)
  {
    const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(model->getAnam());
    r = anamH->getRCoef();
  }

  /* Transform the non conditional simulation into simulation error */

  if (! flag_pgs)
  {
    /********************************/
    /* Standard case (multivariate) */
    /********************************/

    /* Processing */

    for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
    {
      if (!dbin->isActive(iech)) continue;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        double zvar = TEST;
        if (!flag_gibbs)
        {
          zvar = dbin->getLocVariable(ELoc::Z,iech, ivar);
        }
        for (int isimu = 0; isimu < nbsimu; isimu++)
        {
          if (flag_gibbs)
          {
            zvar = dbin->getSimvar(ELoc::GAUSFAC, iech, isimu, ivar, 0, nbsimu,
                                   nvar);
            if (OptDbg::query(EDbg::SIMULATE))
              tab_printg(NULL, zvar);
          }
          double simval = dbin->getSimvar(ELoc::SIMU, iech, isimu, ivar, icase,
                                          nbsimu, nvar);
          if (flag_dgm)
          {
            simval = r * simval + sqrt(1. - r * r) * law_gaussian();
          }
          double simunc = (!FFFF(zvar) && !FFFF(simval)) ? simval - zvar : TEST;
          dbin->setSimvar(ELoc::SIMU, iech, isimu, ivar, icase, nbsimu, nvar,
                          simunc);
        }
      }
    }
  }
  else
  {

    /*********************************************************/
    /* Case of PGS: Data varies per simulation (monovariate) */
    /*********************************************************/

    /* Processing */

    for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
    {
      if (!dbin->isActive(iech)) continue;
      for (int isimu = 0; isimu < nbsimu; isimu++)
      {
        double zvar = dbin->getSimvar(ELoc::GAUSFAC, iech, isimu, 0, icase,
                                      nbsimu, 1);
        if (!FFFF(zvar))
          dbin->updSimvar(ELoc::SIMU, iech, isimu, 0, icase, nbsimu, 1, 0,
                          -zvar);
      }
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Correct for the mean in the case of non-conditional simulations
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  icase     Rank of PGS or GRF
 **
 *****************************************************************************/
void CalcSimuTurningBands::_meanCorrect(Db *dbout, int icase)
{
  int nbsimu = getNbSimu();
  int nvar = _getNVar();
  int nech = dbout->getSampleNumber();
  int ecr = 0;

  // Loop on the simulations
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {

    // Loop on the variables
    for (int ivar = 0; ivar < nvar; ivar++, ecr++)
    {

      // Loop on the samples

      for (int iech = 0; iech < nech; iech++)
      {
        if (!dbout->isActive(iech)) continue;
        dbout->updSimvar(ELoc::SIMU, iech, isimu, ivar, icase, nbsimu,
                         nvar, 0, getModel()->getMean(ivar));
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Update the conditional simulations when the target coincides
 **  with a data point
 **
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  icase     Case for PGS or GRF
 ** \param[in]  flag_pgs  1 if called from PGS
 ** \param[in]  flag_dgm  1 for the Discrete Gaussian Model
 **
 ** \remarks This migration is not performed in the case where data point
 ** \remarks coincide with the target artificially. This is the case
 ** \remarks for the Discrete Gaussian Model (DGM) where data have been
 ** \remarks migrated to the cell center to mimic a point randomized
 ** \remarks within a cell
 **
 *****************************************************************************/
void CalcSimuTurningBands::_updateData2ToTarget(Db *dbin,
                                                Db *dbout,
                                                int icase,
                                                bool flag_pgs,
                                                bool flag_dgm)
{
  if (dbin->getSampleNumber() <= 0) return;
  if (flag_dgm) return;
  int nvar = _getNVar();
  int ndim = dbin->getNDim();
  int nbsimu = getNbSimu();

  /* Calculate the field extension */

  double radius = dbin->getExtensionDiagonal();
  double eps = radius * EPSILON6;
  double eps2 = eps * eps;
  VectorDouble coor1(ndim);
  VectorDouble coor2(ndim);

  /* Dispatch according to the file type */

  if (dbout->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);

    /*********************************************/
    /* Case where the output file is a grid file */
    /*********************************************/

    for (int ip = 0; ip < dbin->getSampleNumber(); ip++)
    {
      if (!dbin->isActive(ip)) continue;
      dbin->getSampleCoordinates(ip, coor2);
      int rank = dbgrid->coordinateToRank(coor2, false, eps);
      if (rank < 0 || !dbgrid->isActive(rank)) continue;
      dbgrid->rankToCoordinateInPlace(rank, coor1);

      /* Get the distance to the target point */

      double dist = 0;
      for (int idim = 0; idim < ndim; idim++)
      {
        double delta = coor1[idim] - coor2[idim];
        dist += delta * delta;
      }
      if (dist > eps2) continue;

      /* We have found a close data point: perform the assignment */

      for (int isimu = 0; isimu < nbsimu; isimu++)
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          double valdat;
          if (!flag_pgs)
            valdat = dbin->getLocVariable(ELoc::Z,ip, ivar);
          else
            valdat = dbin->getSimvar(ELoc::GAUSFAC, ip, isimu, 0, icase, nbsimu,
                                     1);
          if (FFFF(valdat)) continue;
          dbgrid->setSimvar(ELoc::SIMU, rank, isimu, ivar, icase, nbsimu, nvar,
                           valdat);
        }
    }
  }
  else
  {

    /**********************************************/
    /* Case where the output file is a point file */
    /**********************************************/

    for (int ik = 0; ik < dbout->getSampleNumber(); ik++)
    {
      if (!dbout->isActive(ik)) continue;
      dbin->getSampleCoordinates(ik, coor1);

      /* Look for the closest data point */

      int ip_close = -1;
      for (int ip = 0; ip < dbin->getSampleNumber() && ip_close < 0; ip++)
      {
        if (!dbin->isActive(ip)) continue;
        dbin->getSampleCoordinates(ip, coor2);

        /* Get the distance to the target point */

        double dist = 0;
        for (int idim = 0; idim < ndim; idim++)
        {
          double delta = coor1[idim] - coor2[idim];
          dist += delta * delta;
        }
        if (dist <= eps2) ip_close = ip;
      }

      if (ip_close < 0) continue;

      /* We have found a close data point: perform the assignment */

      for (int isimu = 0; isimu < nbsimu; isimu++)
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          double valdat;
          if (!flag_pgs)
            valdat = dbin->getLocVariable(ELoc::Z,ip_close, ivar);
          else
            valdat = dbin->getSimvar(ELoc::GAUSFAC, ip_close, isimu, 0, icase,
                                     nbsimu, 1);
          if (FFFF(valdat)) continue;
          dbout->setSimvar(ELoc::SIMU, ik, isimu, ivar, icase, nbsimu, nvar,
                           valdat);
        }
    }
  }
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
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  icase      Case for PGS or -1
 ** \param[in]  flag_bayes 1 if the Bayes option is switched ON
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  flag_pgs   1 if called from PGS
 ** \param[in]  flag_gibbs 1 if called from Gibbs
 ** \param[in]  flag_dgm   1 if the Discrete Gaussian Model is used
 **
 *****************************************************************************/
int CalcSimuTurningBands::simulate(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   ANeighParam *neighparam,
                                   int icase,
                                   int flag_bayes,
                                   const VectorDouble &dmean,
                                   const VectorDouble &dcov,
                                   bool flag_pgs,
                                   bool flag_gibbs,
                                   bool flag_dgm)
{
  setDbin(dbin);
  setDbout(dbout);
  setModel(model);
  setNeighparam(neighparam);
  setIcase(icase);
  setFlagBayes(flag_bayes);
  setBayesMean(dmean);
  setBayesCov(dcov);
  setFlagPgs(flag_pgs);
  setFlagGibbs(flag_gibbs);
  setFlagDgm(flag_dgm);

  if (! _run()) return 1;
  return 0;
}

bool CalcSimuTurningBands::_run()
{
  law_set_random_seed(getSeed());
  bool flag_cond = hasDbin(false);
  int nbsimu = getNbSimu();

  // Initializations

  if (! _resize()) return false;
  if (_generateDirections(getDbout())) return false;
  _minmax(getDbout());
  _minmax(getDbin());
  if (_initializeSeedBands()) return false;

  // Calculate the 'aic' array

  VectorDouble aic = _createAIC();
  if (aic.empty()) return false;

  // Non conditional simulations on the data points

  if (flag_cond)
  {
    _simulatePoint(getDbin(), aic, _icase, 0);

    // Calculate the simulated error

    _difference(getDbin(), getModel(), _icase, _flagPGS, _flagGibbs, _flagDGM);
  }

  // Non conditional simulations on the grid

  if (getDbout()->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    _simulateGrid(dbgrid, aic, _icase, 0);
  }
  else
  {
    _simulatePoint(getDbout(), aic, _icase, 0);
  }

  /* Add the contribution of Nugget effect (optional) */

  _simulateNugget(getDbout(), aic, _icase);

  /* Conditional simulations */

  if (flag_cond)
  {
    if (_krigsim(getDbin(), getDbout(), getModel(), getNeighparam(),
                 _flagBayes, _bayesMean, _bayesCov, _icase,
                 nbsimu, _flagDGM)) return 1;
  }
  else
  {

    /* In non-conditional case, correct for the mean */

    _meanCorrect(getDbout(), _icase);
  }

  /* Copy value from data to coinciding grid node */

  if (flag_cond)
    _updateData2ToTarget(getDbin(), getDbout(), _icase, _flagPGS, _flagDGM);

  // Check the simulation at data location

  if (_flagCheck)
    _checkGaussianData2Grid(getDbin(), getDbout(), getModel());

  return true;
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
int CalcSimuTurningBands::simulatePotential(Db *dbiso,
                                            Db *dbgrd,
                                            Db *dbtgt,
                                            Db *dbout,
                                            Model *model,
                                            double delta)
{
  setDbout(dbout);
  setModel(model);
  if (getNbSimu() <= 0 || getNBtuba() <= 0)
  {
    messerr("You must define 'nbsimu', 'nbtuba' and the 'model' beforehand");
    return 1;
  }
  law_set_random_seed(getSeed());
  int icase = 0;

  /* Processing the Turning Bands algorithm */

  if (_generateDirections(dbout)) return 1;
  _minmax(dbout);
  _minmax(dbiso);
  _minmax(dbgrd);
  _minmax(dbtgt);
  if (_initializeSeedBands()) return 1;

  /* Calculate the 'aic' array */

  VectorDouble aic = _createAIC();
  if (aic.empty()) return 1;

  /* Non conditional simulations on the data points */

  if (dbiso != nullptr)
  {
    _simulatePoint(dbiso, aic, icase, 0);
  }

  /* Non conditional simulations on the gradient points */

  if (dbgrd != nullptr)
  {
    _simulateGradient(dbgrd, aic, delta);
  }

  /* Non conditional simulations on the tangent points */

  if (dbtgt != nullptr)
  {
    _simulateTangent(dbtgt, aic, delta);
  }

  /* Non conditional simulations on the grid */

  if (dbout->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    _simulateGrid(dbgrid, aic, icase, 0);
  }
  else
  {
    _simulatePoint(dbout, aic, icase, 0);
  }

  /* Add the contribution of nugget effect (optional) */

  _simulateNugget(dbout, aic, icase);
  return 0;
}

/****************************************************************************/
/*!
 **  Check if the Model can be simulated using Turning Bands
 **
 ** \return  True if the Model is valid; 0 otherwise
 **
 ** \param[in]  model    Model structure
 **
 *****************************************************************************/
bool CalcSimuTurningBands::isTurningBandsWorkable(const Model *model)

{
  bool workable = true;

  /* Loop on the structures */

  for (int is = 0; is < model->getCovaNumber(); is++)
  {
    ECov type = model->getCovaType(is);

    switch (type.toEnum())
    {
      case ECov::E_NUGGET:
      case ECov::E_EXPONENTIAL:
      case ECov::E_SPHERICAL:
      case ECov::E_CUBIC:
      case ECov::E_GAUSSIAN:
      case ECov::E_SINCARD:
      case ECov::E_BESSEL_J:
      case ECov::E_BESSEL_K:
      case ECov::E_STABLE:
      case ECov::E_POWER:
      case ECov::E_SPLINE_GC:
      case ECov::E_LINEAR:
      case ECov::E_ORDER1_GC:
      case ECov::E_ORDER3_GC:
      case ECov::E_ORDER5_GC:
        break;

      default:
        workable = false;
    }
  }
  return workable;
}

/****************************************************************************/
/*!
 **  Check/Show the data (gaussian) against the closest grid node
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db grid structure
 ** \param[in]  model      Model structure
 **
 ** \remark Attributes ELoc::SIMU and ELoc::GAUSFAC (for PGS) are mandatory
 ** \remark Tests have only been produced for icase=0
 **
 *****************************************************************************/
void CalcSimuTurningBands::_checkGaussianData2Grid(Db *dbin,
                                                   Db *dbout,
                                                   Model *model) const
{
  if (dbin == nullptr) return;
  if (get_LOCATOR_NITEM(dbout,ELoc::SIMU) <= 0) return;
  int nbsimu = getNbSimu();
  if (nbsimu <= 0) return;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  if (dbgrid == nullptr) return;
  int ndim = dbin->getNDim();

  mestitle(1, "Checking Gaussian of data against closest grid node");


  /* Loop on the data */

  int number = 0;
  VectorDouble coor(ndim);
  for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (!dbin->isActive(iech)) continue;

    // Find the index of the closest grid node and derive tolerance
    int jech = index_point_to_grid(dbin, iech, 0, dbgrid, coor.data());
    if (jech < 0) continue;
    double eps = model_calcul_stdev(model, dbin, iech, dbgrid, jech, 0, 2.);
    if (eps < 1.e-6) eps = 1.e-6;

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double valdat = dbin->getSimvar(ELoc::GAUSFAC, iech, 0, 0, 0, nbsimu, 1);
      double valres = dbgrid->getSimvar(ELoc::SIMU, jech, isimu, 0, 0, nbsimu, 1);
      if (ABS(valdat - valres) < eps) continue;
      number++;

      /* The data facies is different from the grid facies */

      message("Inconsistency for Simulation (%d) between :\n", isimu + 1);
      message("- Value (%lf) at Data (#%d) ", valdat, iech + 1);
      message("at (");
      for (int idim = 0; idim < ndim; idim++)
        message(" %lf", dbin->getCoordinate(iech, idim));
      message(")\n");

      message("- Value (%lf) at Grid (#%d) ", valres, jech + 1);
      message("at (");
      for (int idim = 0; idim < ndim; idim++)
        message(" %lf", dbgrid->getCoordinate(jech, idim));
      message(")\n");

      message("- Tolerance = %lf\n", eps);
    }
  }
  if (number <= 0) message("No problem found\n");
  return;
}

bool CalcSimuTurningBands::_check()
{
  if (! ACalcSimulation::_check()) return false;

  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (hasDbin(false))
  {
    if (! hasNeighParam()) return false;
  }
  int ndim = _getNDim();
  if (ndim > 3)
  {
    messerr("The Turning Band Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return false;
  }
  if (getNBtuba() <= 0)
  {
    messerr("You must define 'nbsimu' and 'nbtuba'");
    return 1;
  }

  if (_flagDGM)
  {
    if (! getDbout()->isGrid())
    {
      messerr("For DGM option, the argument 'dbout'  should be a Grid");
      return false;
    }
    if (! getModel()->hasAnam())
    {
      messerr("For DGM option, the Model must have an Anamorphosis attached");
      return false;
    }
    if (! getModel()->isChangeSupportDefined())
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

  int nvar = _getNVar();
  int nbsimu = getNbSimu();

  /* Add the attributes for storing the results */

  if (getDbin() != nullptr)
  {
    int iptr_in = _addVariableDb(1, 2, ELoc::SIMU, 0, nvar * nbsimu);
    if (iptr_in < 0) return false;
  }

  _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 0, nvar * nbsimu);
  if (_iattOut < 0) return false;

  // Centering the Data (for DGM)

  if (_flagDGM)
  {
    // Centering (only if the output file is a Grid)
    DbGrid *dbgrid = dynamic_cast<DbGrid*>(getDbout());
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

  _renameVariable(2, _getNVar(), _iattOut, String(), getNbSimu());

  if (_flagDGM)
  {
    if (!_nameCoord.empty())
      getDbin()->setLocators(_nameCoord, ELoc::X);
  }

  return true;
}

void CalcSimuTurningBands::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional simulation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure (optional)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  seed       Seed for random number generator
 ** \param[in]  nbtuba     Number of turning bands
 ** \param[in]  flag_dgm   1 for Direct Block Simulation
 ** \param[in]  flag_check 1 to check the proximity in Gaussian scale
 ** \param[in]  namconv    Naming convention
 **
 ** \remark  The arguments 'dbout' and 'neigh' are optional: they must
 ** \remark  be defined only for conditional simulations
 **
 *****************************************************************************/
int simtub(Db *dbin,
           Db *dbout,
           Model *model,
           ANeighParam *neighparam,
           int nbsimu,
           int seed,
           int nbtuba,
           bool flag_dgm,
           bool flag_check,
           const NamingConvention &namconv)
{
  // Instantiate the Calculator
  CalcSimuTurningBands situba(nbsimu, nbtuba, flag_check, seed);

  // Set the members of the Calculator
  situba.setDbin(dbin);
  situba.setDbout(dbout);
  situba.setModel(model);
  situba.setNeighparam(neighparam);
  situba.setNamingConvention(namconv);
  situba.setFlagDgm(flag_dgm);

  // Run the calculator
  int error = (situba.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional simulation
 **  with Bayesian Drift
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure (optional)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  seed       Seed for random number generator
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  nbtuba     Number of turning bands
 ** \param[in]  flag_check 1 to check the proximity in Gaussian scale
 ** \param[in]  namconv    Naming convention
 **
 ** \remark  The arguments 'dbout' and 'neigh' are optional: they must
 ** \remark  be defined only for conditional simulations
 **
 *****************************************************************************/
int simbayes(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             int nbsimu,
             int seed,
             const VectorDouble& dmean,
             const VectorDouble& dcov,
             int nbtuba,
             bool flag_check,
             const NamingConvention& namconv)
{
  // Instantiate the Calculator
  CalcSimuTurningBands situba(nbsimu, nbtuba, flag_check, seed);

  // Set the members of the Calculator
  situba.setDbin(dbin);
  situba.setDbout(dbout);
  situba.setModel(model);
  situba.setNeighparam(neighparam);
  situba.setNamingConvention(namconv);

  situba.setFlagBayes(true);
  situba.setBayesMean(dmean);
  situba.setBayesCov(dcov);

  // Run the calculator
  int error = (situba.run()) ? 0 : 1;
  return error;
}

