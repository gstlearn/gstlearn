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
#include "geoslib_old_f.h"

#include "Geometry/GeometryHelper.hpp"
#include "Fractures/FracList.hpp"
#include "Fractures/FracDesc.hpp"
#include "Fractures/FracEnviron.hpp"
#include "Fractures/FracFamily.hpp"
#include "Fractures/FracFault.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"

#include <math.h>

#define FRACSEG(ifrac,i)    (frac_segs[NBYFRAC * (ifrac) + (i)])
#define WELL(iw,i)          (well[2*(iw) + (i)])
#define WELLOUT(is,i)       (wellout[NBYWOUT * (is) + (i)])
#define TRAJ(ip,i)          (traj[2*(ip)+(i)])

FracList::FracList(int ndisc, bool flag_check, double low0, double low1, double eps)
  : AStringable(),
    _descs(),
    _layinfo(),
    _nlayers(0),
    _ndisc(ndisc),
    _flagCheck(flag_check),
    _low0(low0),
    _low1(low1),
    _xorigin(0.),
    _step(0.),
    _eps(eps),
    _verbose(false)
{
}

FracList::FracList(const FracList& r)
    : AStringable(r),
      _descs(r._descs),
      _layinfo(r._layinfo),
      _nlayers(r._nlayers),
      _ndisc(r._ndisc),
      _flagCheck(r._flagCheck),
      _low0(r._low0),
      _low1(r._low1),
      _xorigin(r._xorigin),
      _step(r._step),
      _eps(r._eps),
      _verbose(r._verbose)
{
}

FracList& FracList::operator=(const FracList& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _descs = r._descs;
    _layinfo = r._layinfo;
    _nlayers = r._nlayers;
    _ndisc = r._ndisc;
    _flagCheck = r._flagCheck;
    _low0 = r._low0;
    _low1 = r._low1;
    _xorigin = r._xorigin;
    _step = r._step;
    _eps = r._eps;
    _verbose = r._verbose;
  }
  return *this;
}

FracList::~FracList()
{
}

String FracList::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNFracs() <= 0) return sstr.str();

  mestitle(0,"Fracture Description");
  sstr << "Current number of simulated fractures = " <<  getNFracs() << std::endl;

  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
    sstr << _descs[ifrac].toString(strfmt);

  return sstr.str();
}

/****************************************************************************/
/*!
 **  Simulate a set of fractures
 **
 ** \param[in]  envir          Environ structure
 ** \param[in]  flag_sim_layer TRUE for simulating layers
 **                            FALSE if they are read
 ** \param[in]  flag_sim_fract TRUE for simulating the fractures
 **                            FALSE establish the layers and the main faults
 ** \param[in]  elevations     Array of elevations (used if flag_sim_layer=F)
 ** \param[in]  seed           Seed for the random number generator
 ** \param[in]  verbose        Verbose option
 **
 *****************************************************************************/
int FracList::simulate(const FracEnviron& envir,
                       bool flag_sim_layer,
                       bool flag_sim_fract,
                       int seed,
                       bool verbose,
                       const VectorDouble& elevations)
{
  _step = envir.getXextend() / _ndisc;
  _xorigin = -envir.getDeltax();
  _verbose = verbose;

  double y0 = 0.;
  int ninfos = 1 + envir.getNFamilies() * NPART;
  law_set_random_seed(seed);

  if (_verbose)
  {
    message("Fracture_Discretization_Count = %d \n", _ndisc);
    message("Fracture_Check_Intersect      = %d \n", _flagCheck);
    message("Fracture_Repulsion_Low0       = %lg\n", _low0);
    message("Fracture_Repulsion_Low1       = %lg\n", _low1);
  }

  /* Determine the layers */

  VectorDouble thicks;
  if (flag_sim_layer)
    thicks = _layersManage(envir, &y0);
  else
    thicks = _layersRead(elevations, &y0);
  _nlayers = (int) thicks.size();

  /* Core allocation */

  int np1 = _nlayers + 1;
  VectorDouble denstab(_ndisc);
  _layinfo = MatrixRectangular(np1, ninfos);

  /* Define the layers */

  double cote = y0;
  for (int ilayer = 0; ilayer < _nlayers; ilayer++)
  {
    double thick = thicks[ilayer];
    _setMemLayer(ilayer, cote);
    cote += thick;
  }
  _setMemLayer(_nlayers, cote);

  /* Define the main faults */

  for (int ifault = 0; ifault < envir.getNFaults(); ifault++)
  {
    double angle = envir.getFault(ifault).getOrient();

    /* Loop on the layers */

    int ifrac = -1;
    cote = y0;
    double xx = envir.getFault(ifault).faultAbscissae(cote);
    for (int ilayer = 0; ilayer < _nlayers; ilayer++)
    {
      double thick = thicks[ilayer];
      ifrac = _fracAdd(ifrac, 0, xx, cote, thick, angle, &xx);
      cote += thick;
    }
  }

  if (_verbose)
    message("Number of main faults        = %d \n", getNFracs());

  if (!flag_sim_fract) return 0;

  /* Loop on the different fault families */

  for (int ifam = 0; ifam < envir.getNFamilies(); ifam++)
  {
    if (_verbose)
      mestitle(0, "Processing Family #%d/%d", ifam + 1, envir.getNFamilies());
    const FracFamily& family = envir.getFamily(ifam);

    /* Loop on the layers */

    cote = y0;
    double thetap = 0.;
    for (int ilayer = 0; ilayer < _nlayers; ilayer++)
    {
      double thick = thicks[ilayer];
      if (_verbose)
      {
        mestitle(1, "Processing Layer #%d/%d", ilayer + 1, _nlayers);
        message("Elevation of the layer bottom     = %lf\n", cote);
        message("Thickness of the layer            = %lf\n", thick);
      }

      /* Derive the layer intensity */

      double theta1 = _layerIntensity(family, thick);

      /* Generate the density function */

      _generateDensity(envir, family, ifam, cote, denstab);

      /* Correct the density due to already existing faults */

      _correctDensity(family, ifam, cote, denstab);

      /* Extend previously existing fractures */

      double propsur = _extendFractures(family, ifam, cote, thick, denstab);

      /* Derive the layer intensity, given the survival from previous layers */

      double theta2 = _deriveIntensity(theta1, thetap, propsur);

      /* Simulate the fractures within the layer */

      int nfracs = _simulateFractures(envir, family, ifam, cote, thick, theta2, denstab);

      /* Shift the ordinate */

      cote += thick;
      thetap = theta1;
      _setMemTheta1(ilayer,ifam,theta1);
      _setMemTheta2(ilayer,ifam,theta2);
      _setMemPropsur(ilayer,ifam,propsur);
      _setMemFrac(ilayer,ifam,(double) nfracs);
      _setMemTotal(ilayer,ifam,(double) getNFracs());
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Generate the series of layers
 **
 ** \return The vector of layer thicknesses
 **
 ** \param[in]  envir      Pointer to the Environ structure
 ** \param[in]  y0         Ordinate of the origin
 **
 ** \remarks The origin is calculated so that a layer edge is located at 0.
 **
 *****************************************************************************/
VectorDouble FracList::_layersManage(const FracEnviron& envir,
                                     double *y0)
{
  VectorDouble thicks;
  double thick_min = envir.getMean() / 10.;
  double mean = envir.getMean();
  double stdev = envir.getStdev();
  double var = stdev * stdev;
  double theta = var / mean;
  bool flag_rnd = var > 0.;
  double k = (flag_rnd) ? (mean * mean) / var : 0.;

  /* Downwards */

  double total = 0.;
  while (total < envir.getDeltay())
  {
    double thick = (flag_rnd) ? theta * law_gamma(k) : mean;
    if (thick < thick_min) continue;
    total += thick;
    thicks.push_back(thick);
  }
  double cote_down = -total;

  /* Upwards */

  total = 0.;
  while (total < envir.getYmax())
  {
    double thick = (flag_rnd) ? theta * law_gamma(k) : mean;
    if (thick < thick_min) continue;
    total += thick;
    thicks.push_back(thick);
  }
  double cote_up = total;

  *y0 = cote_down;

  /* Optional information */

  if (_verbose)
  {
    int number = (int) thicks.size();
    mestitle(0, "Layer generation");
    message("Thickness law - Mean               = %lf\n", mean);
    message("Thickness law - St. Dev.           = %lf\n", stdev);
    message("Minimum simulated level            = %lf\n", cote_down);
    message("Maximum simulated level            = %lf\n", cote_up);
    message("Number of layers                   = %d \n", number);
  }

  return thicks;
}

/****************************************************************************/
/*!
 **  Read the series of layers
 **
 ** \param[in]  elevations   Array of elevations to be read
 **
 ** \param[out] y0           Ordinate of the origin
 **
 *****************************************************************************/
VectorDouble FracList::_layersRead(const VectorDouble& elevations, double *y0)
{
  VectorDouble thicks;

  /* Initializations */

  int nlayers = (int) elevations.size();
  thicks.resize(nlayers - 1);

  double cote = elevations[0];
  for (int i = 1; i < nlayers; i++)
  {
    thicks[i - 1] = elevations[i] - cote;
    cote = elevations[i];
  }
  double cote_down = (*y0) = elevations[0];
  double cote_up = elevations[nlayers - 1];

  /* Optional information */

  if (_verbose)
  {
    mestitle(0, "Layer (read)");
    message("Minimum simulated level            = %lf\n", cote_down);
    message("Maximum simulated level            = %lf\n", cote_up);
    message("Number of layers                   = %d \n", nlayers);
  }

  return thicks;
}

/****************************************************************************/
/*!
 **  Add a fracture
 **
 ** \return Rank of the fracture
 **
 ** \param[in]  ifrac        Rank of the fracture (if old) or -1 for a new one
 ** \param[in]  ifam         Rank of the family
 ** \param[in]  xx           Abscissae of the starting point
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  thick        Thickness of the layer
 ** \param[in]  orient       Orientation
 **
 ** \param[out] xp           Abscissae of the ending point
 **
 *****************************************************************************/
int FracList::_fracAdd(int ifrac,
                       int ifam,
                       double xx,
                       double cote,
                       double thick,
                       double orient,
                       double* xp)
{
  if (ifrac < 0)
  {
    addDescription();
    ifrac = getNFracs() - 1;
  }

  /* Set the attributes of the new segment */

  FracDesc& desc = _descs[ifrac];

  if (desc.getNPoint() == 0)
  {
    desc.addPoint(xx, cote);
  }
  desc.setFamily(ifam);
  desc.setOrient(orient);
  desc.addPoint(thick * tan(ut_deg2rad(orient)) + xx, thick + cote);
  *xp = desc.getXXF(desc.getNPoint() - 1);

  if (_verbose)
    message("- Adding fracture: (%lf; %lf) to (%lf; %lf)\n",
            xx, cote,
            desc.getXXF(desc.getNPoint() - 1),
            desc.getYYF(desc.getNPoint() - 1));
  _checkFractureIntersect(cote, ifrac);

  return (ifrac);
}

/****************************************************************************/
/*!
 **  Check fracture against the other ones and possibly intersect the current one
 **
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  ifrac0       Rank of the current fracture
 **
 *****************************************************************************/
void FracList::_checkFractureIntersect(double cote, int ifrac0)
{
  double xd2, yd2, xe2, ye2, x, y;

  if (! _flagCheck) return;
  FracDesc& desc = _descs[ifrac0];
  int npoint = desc.getNPoint();
  double xd1 = desc.getXXF(npoint - 2);
  double yd1 = desc.getYYF(npoint - 2);
  double xe1 = desc.getXXF(npoint - 1);
  double ye1 = desc.getYYF(npoint - 1);

  /* Check against the previous fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    if (ifrac == ifrac0) continue;
    FracDesc& desc2 = _descs[ifrac];

    /* Look for the segment belonging to the current layer */

    if (! _belongToLayer(desc2, cote, &xd2, &yd2, &xe2, &ye2)) continue;

    /* Perform the intersection between two segments */

    if (! GH::segmentIntersect(xd1, yd1, xe1, ye1, xd2, yd2, xe2, ye2, &x, &y))
      continue;

    /* Update the end-point in case of intersection */

    desc.setXXF(npoint - 1, x);
    desc.setYYF(npoint - 1, y);
  }
}

/****************************************************************************/
/*!
 **  Search for the segment of a fracture belonging to the current layer
 **
 ** \return 1 if the segment belongs to the layer; 0 otherwise
 **
 ** \param[in]  desc         Description structure
 ** \param[in]  cote         Ordinate of the fracture starting point
 **
 ** \param[in]  xd           Abscissae of the first end-point
 ** \param[in]  yd           Ordinate of the first end-point
 ** \param[in]  xe           Abscissae of the second end-point
 ** \param[in]  ye           Ordinate of the second end-point
 **
 *****************************************************************************/
bool FracList::_belongToLayer(const FracDesc& desc,
                              double cote,
                              double *xd,
                              double *yd,
                              double *xe,
                              double *ye)
{
  for (int i = 0; i < desc.getNPoint() - 1; i++)
  {
    if (ABS(desc.getYYF(i) - cote) > _eps) continue;
    *xd = desc.getXXF(i);
    *yd = desc.getYYF(i);
    *xe = desc.getXXF(i + 1);
    *ye = desc.getYYF(i + 1);
    return true;
  }
  return false;
}

/****************************************************************************/
/*!
 **  Derive the intensity of the current layer
 **
 ** \return Layer intensity
 **
 ** \param[in]  family       Family structure
 ** \param[in]  thick        Layer thickness
 **
 *****************************************************************************/
double FracList::_layerIntensity(const FracFamily& family, double thick)
{
  double theta1 = family.getTheta0() / pow(thick, family.getAlpha());

  if (_verbose)
    message("Initial Intensity                 = %lf\n", theta1);

  return (theta1);
}

/****************************************************************************/
/*!
 **  Generate the density along the scene
 **
 ** \param[in]  envir     Environ structure
 ** \param[in]  family    Family structure
 ** \param[in]  ifam      Rank of the family
 ** \param[in]  cote      Ordinate of the fracture starting point
 **
 ** \param[out] denstab   Discretization density array
 **
 *****************************************************************************/
void FracList::_generateDensity(const FracEnviron& envir,
                                const FracFamily& family,
                                int ifam,
                                double cote,
                                VectorDouble& denstab)
{
  double ratcst = family.getRatcst();
  for (int idisc = 0; idisc < _ndisc; idisc++)
    denstab[idisc] = 0.;

  /* Create the shaped density */

  if (ratcst < 1.)
  {

    /* Loop on the discretization steps */

    for (int idisc = 0; idisc < _ndisc; idisc++)
    {
      double x0 = _xorigin + _step * (idisc + 0.5);

      /* Loop on the faults */

      double value = 0.;
      for (int ifault = 0; ifault < envir.getNFaults(); ifault++)
      {
        const FracFault& fault = envir.getFault(ifault);
        if (! _sameFaultSide(envir, ifault, x0)) continue;

        /* Left side */
        value = MAX(value, _densityUpdate(fault, -1, ifam, cote, x0));

        /* Right side */
        value = MAX(value, _densityUpdate(fault, +1, ifam, cote, x0));
      }
      denstab[idisc] = value;
    }
  }

  /* Add the constant and the shaped factors */

  for (int idisc = 0; idisc < _ndisc; idisc++)
    denstab[idisc] = denstab[idisc] * (1. - ratcst) + (1. / _ndisc) * ratcst;

  if (_verbose)
    message("- Cumulated Distribution: Main Fault = %lf\n",
            _densityCumulate(denstab));

  return;
}

/****************************************************************************/
/*!
 **  Correct density due to fractures of previous families (if any)
 **
 ** \param[in]  family       Family structure
 ** \param[in]  ifam         Rank of the family
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  denstab      Discretization density array
 **
 *****************************************************************************/
void FracList::_correctDensity(const FracFamily& family,
                               int ifam,
                               double cote,
                               VectorDouble& denstab)
{
  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    FracDesc& desc = _descs[ifrac];
    if (desc.getFamily() >= ifam + 1) continue;

    /* Loop on the segments of the fault */

    int npoint = desc.getNPoint();
    for (int i = 0; i < npoint - 1; i++)
    {
      if (ABS(desc.getYYF(i) - cote) <= _eps)
        _updateRepulsion(desc.getXXF(i), family.getRange(), denstab);
      if (ABS(desc.getYYF(i+1) - cote) <= _eps)
        _updateRepulsion(desc.getXXF(i + 1), family.getRange(), denstab);
    }
  }

  if (_verbose)
    message("- Cumulated Distribution: Previous families = %lf\n",
            _densityCumulate(denstab));
}

/****************************************************************************/
/*!
 **  Derive the intensity of the current layer given the survival
 **
 ** \return Layer actual intensity
 **
 ** \param[in]  theta1       Apparent layer intensity
 ** \param[in]  thetap       Intensity of the previous layer
 ** \param[in]  propsur      Actual survival proportion
 **
 *****************************************************************************/
double FracList::_deriveIntensity(double theta1,
                                  double thetap,
                                  double propsur)
{
  double theta2 = theta1;

  if (thetap > 0) theta2 = MAX(0., theta1 - propsur * thetap);

  if (_verbose) message("Intensity corrected from survival = %lf\n", theta2);

  return (theta2);
}

/****************************************************************************/
/*!
 **  Extend already existing fractures
 **
 ** \return The proportion of survival
 **
 ** \param[in]  family        Family structure
 ** \param[in]  ifam          Rank of the family
 ** \param[in]  cote          Ordinate of the fracture starting point
 ** \param[in]  thick         Thickness of the layer
 ** \param[in]  denstab       Discretization density array
 **
 *****************************************************************************/
double FracList::_extendFractures(const FracFamily& family,
                                  int ifam,
                                  double cote,
                                  double thick,
                                  VectorDouble& denstab)
{
  int next = 0;
  int ntotal = 0;
  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    const FracDesc& desc = _descs[ifrac];
    int npoint = desc.getNPoint();
    if (desc.getFamily() != ifam + 1) continue;
    if (ABS(desc.getYYF(npoint-1) - cote) > _eps) continue;
    double angle = desc.getOrient();

    /* Evaluate the persistence probability */

    ntotal++;
    if (_fractureInterrupt(family, desc, thick)) continue;

    /* The fracture is extended */

    next++;
    double xx = desc.getXXF(npoint - 1);
    double xp;
    (void) _fracAdd(ifrac, ifam+1, xx, cote, thick, angle, &xp);

    /* Update the density to account for repulsion */

    _updateRepulsion(xx, family.getRange(), denstab);
  }

  double propsur = (ntotal > 0) ? (double) next / (double) ntotal : 0.;
  if (_verbose && ntotal > 0)
  {
    message("Survival proportion               = %lf (%d from %d)\n", propsur,
            next, ntotal);
    message("- Cumulated Distribution: After survival = %lf\n",
            _densityCumulate(denstab));
  }
  return (propsur);
}

/****************************************************************************/
/*!
 **  Check that the current discretization point is not separated from the
 **  current fault by other main faults
 **
 ** \return The discretized point is not separated by another fault
 **
 ** \param[in]  envir      Environ structure
 ** \param[in]  ifault0    Rank of the target fault
 ** \param[in]  x0         Location of the discretized point
 **
 *****************************************************************************/
bool FracList::_sameFaultSide(const FracEnviron& envir,
                              int ifault0,
                              double x0)
{
  double x1, x2;
  int ifault;

  x1 = envir.getFault(ifault0).getCoord();
  for (ifault = 0; ifault < envir.getNFaults(); ifault++)
  {
    if (ifault == ifault0) continue;
    x2 = envir.getFault(ifault).getCoord();

    if ((x0 - x2) * (x2 - x1) >= 0) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Project the contribution of a density attached to a fault
 **
 ** \return Density shape due to the Fault
 **
 ** \param[in]  fault  Fault structure
 ** \param[in]  side   -1 for left and +1 for right
 ** \param[in]  ifam   Rank of the family
 ** \param[in]  cote   Ordinate of the fracture starting point
 ** \param[in]  xx     Discretization abscissae
 **
 ** \remarks The density shape function is analoguous to the Cubic Covariance
 **
 *****************************************************************************/
double FracList::_densityUpdate(const FracFault& fault,
                                int side,
                                int ifam,
                                double cote,
                                double xx)
{
  double range;
  double theta;

  double x0 = _faultAbscissae(fault, cote);
  if (side < 0)
  {
    if (xx > x0) return (0.);
    range = fault.getRanger(ifam);
    theta = fault.getThetal(ifam);
  }
  else
  {
    if (xx < x0) return (0.);
    range = fault.getRanger(ifam);
    theta = fault.getThetar(ifam);
  }

  double dist = ABS(x0 - xx);
  double h = dist / range;
  return (theta * _cubic(h));
}

/****************************************************************************/
/*!
 **  Calculate the cumulated density
 **
 ** \returns Total cumulated density
 **
 ** \param[in]  denstab      Discretized density array
 **
 *****************************************************************************/
double FracList::_densityCumulate(const VectorDouble& denstab)
{
  double total = 0.;
  for (int idisc = 0; idisc < _ndisc; idisc++)
    total += denstab[idisc];
  return total;
}

/****************************************************************************/
/*!
 **  Calculate the cumulated density
 **
 ** \returns True if there is no more room for fractures
 **
 ** \param[in]  denstab      Discretized density array
 **
 *****************************************************************************/
bool FracList::_noRoomForMoreFracture(const VectorDouble& denstab) const
{
  double total = 0.;
  for (int idisc = 0; idisc < _ndisc; idisc++)
    total += denstab[idisc];
  bool flag_stop = (total <= _low1 * _ndisc);

  if (_verbose && flag_stop)
  {
    message("Fracture simulation interrupted: no more room available\n");
    message("Total = %lf < low1(%lf) * _ndisc(%d)\n", total, _low1, _ndisc);
  }
  return (flag_stop);
}

/****************************************************************************/
/*!
 **  Update the discretized intensity array to account for repulsion
 **
 ** \param[in]  x0           Abscissae of the fracture
 ** \param[in]  range        Range for the repulsion
 ** \param[in,out] denstab   Discretized density array
 **
 ** \remarks  We first designate the central cell of the discretized density
 ** \remarks  array where the new fault is located
 ** \remarks  The density of the central cell is set to LOW0
 ** \remarks  The density of the peripheral cells is set to LOW1
 ** \remarks  A cell is peripheral if located in the 'range' of the central one
 **
 *****************************************************************************/
void FracList::_updateRepulsion(double x0,
                                double range,
                                VectorDouble& denstab)
{
  int idisc0 = (int) ((x0 - _xorigin) / _step);
  int nrange = MAX(1, (int ) (range / _step));

  /* Central cell */
  if (_isValidDisc(idisc0)) denstab[idisc0] = _low0;

  /* Peripheral cell on the left */
  for (int i = 0; i < nrange; i++)
  {
    int idisc = idisc0 - i - 1;
    if (_isValidDisc(idisc)) denstab[idisc] = _low1;
  }

  /* Peripheral cell on the right */
  for (int i = 0; i < nrange; i++)
  {
    int idisc = idisc0 + i + 1;
    if (_isValidDisc(idisc)) denstab[idisc] = _low1;
  }
}

/******************************************************************F**********/
/*!
 **  Check if the fracture must be interrupted
 **
 ** \return  1 if the fracture must be interrupted
 **
 ** \param[in]  family       Family structure
 ** \param[in]  desc         Current fracture description
 ** \param[in]  thick        Thickness of the current layer
 **
 *****************************************************************************/
bool FracList::_fractureInterrupt(const FracFamily& family,
                                  const FracDesc& desc,
                                  double thick)
{
  /* Constant probability */
  double prop1 = family.getProp1();

  /* Length dependent probability */
  double prop2 = family.getProp2();

  /* Cumulated fracture length dependent exponent */
  double dist = _fractureExtension(desc, TEST, TEST);
  double expa = exp(-dist / family.getAterm());

  /* Thickness dependent exponent */
  double expb = exp(-thick / family.getBterm());

  /* Final probability */
  double proba = (prop1 + prop2 * (1. - expa)) * expb;

  double rndval = law_uniform(0., 1);
  return (rndval > proba);
}

/****************************************************************************/
/*!
 **  Calculate the abscissae of a fault at a given elevation
 **
 ** \return The fault abscissae
 **
 ** \param[in]  fault  Fault structure
 ** \param[in]  cote   Ordinate of the fracture starting point
 **
 *****************************************************************************/
double FracList::_faultAbscissae(const FracFault& fault, double cote)
{
  return (fault.getCoord() + cote * tan(ut_deg2rad(fault.getOrient())));
}

/****************************************************************************/
/*!
 **  Calculate the attenuation according to a cubic covariance shape
 **
 ** \return Returned value
 **
 ** \param[in]  h      Normalized distance
 **
 *****************************************************************************/
double FracList::_cubic(double h)
{
  double value, h2;

  if (h >= 1) return (0.);

  h2 = h * h;
  value = 1. - h2 * (7. - h * (35. / 4. - h2 * (7. / 2. - 3. / 4. * h2)));

  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the fracture extension
 **
 ** \return  The fracture extension
 **
 ** \param[in]  desc         Current fracture description
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 *****************************************************************************/
double FracList::_fractureExtension(const FracDesc& desc,
                                    double cote,
                                    double dcote)
{
  double dist = 0.;
  for (int i = 0; i < desc.getNPoint() - 1; i++)
  {
    double distx = desc.getXXF(i+1) - desc.getXXF(i);
    double disty = desc.getYYF(i+1) - desc.getYYF(i);
    if (!FFFF(cote) && (desc.getYYF(i) < cote - dcote ||
        desc.getYYF(i+1) < cote - dcote))
      continue;
    dist += sqrt(distx * distx + disty * disty);
  }
  return (dist);
}

/****************************************************************************/
/*!
 **  Simulate the fractures
 **
 ** \return  Number of fracture created in this layer
 **
 ** \param[in]  envir      Environ structure
 ** \param[in]  family     Family structure
 ** \param[in]  ifam       Rank of the family
 ** \param[in]  cote       Ordinate of the fracture starting point
 ** \param[in]  thick      Thickness of the layer
 ** \param[in]  theta      Intensity of the layer
 ** \param[in]  denstab    Discretized density array
 **
 *******************************************************************F**********/
int FracList::_simulateFractures(const FracEnviron& envir,
                                 const FracFamily& family,
                                 int ifam,
                                 double cote,
                                 double thick,
                                 double theta,
                                 VectorDouble& denstab)
{
  double total;
  double xp;

  double orient = family.getOrient();
  double dorient = family.getDorient();
  int nfracs = law_poisson(theta * envir.getXextend());

  /* Loop on the fractures to be simulated */

  int neff = 0;
  for (int ifrac = 0; ifrac < nfracs; ifrac++)
  {

    /* Calculate cumulative density */

    if (_noRoomForMoreFracture(denstab)) break;
    total = _densityCumulate(denstab);

    /* Simulate the location along to the regionalized density */

    double cumdens = total * law_uniform(0., 1.);

    /* Get the rank of the discretized cell */

    int idisc = _getDiscretizedRank(cumdens, denstab);
    double xdeb = _xorigin + _step * (idisc);
    double xfin = _xorigin + _step * (idisc + 1);
    double xx = law_uniform(xdeb, xfin);
    double angle = law_uniform(MAX(-90., orient - dorient),
                               MIN(+90., orient + dorient));

    /* Add a new fracture */

    (void) _fracAdd(-1, ifam + 1, xx, cote, thick, angle, &xp);

    /* Update the density to account for repulsion */

    _updateRepulsion(xx, family.getRange(), denstab);
    neff++;
  }

  if (_verbose)
    message("Number of Simulated fractures in Layer  = %d \n", neff);

  return (nfracs);
}

/****************************************************************************/
/*!
 **  Get the rank of the discretized cell
 **
 ** \return Rank of the discretized cell where the cumulated density reaches
 ** \return the target density
 **
 ** \param[in]  cumdens      Target cumulated density
 ** \param[in]  denstab      Discretized density array
 **
 *****************************************************************************/
int FracList::_getDiscretizedRank(double cumdens, const VectorDouble& denstab)
{
  double local = 0.;
  for (int idisc = 0; idisc < _ndisc; idisc++)
  {
    local += denstab[idisc];
    if (local > cumdens) return (idisc);
  }
  return (_ndisc - 1);
}

/****************************************************************************/
/*!
 **  Export the Fractures
 **
 *****************************************************************************/
MatrixRectangular FracList::fractureExport() const
{
  int ntotal = _getEndPointCount();
  MatrixRectangular segs(ntotal, NBYFRAC);

  /* Loading the fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    int npoint = _descs[ifrac].getNPoint();
    for (int ip = 0; ip < npoint; ip++)
    {
      segs.setValue(ifrac, 0, (double) ifrac+1);
      segs.setValue(ifrac, 1, (double) ip+1);
      segs.setValue(ifrac, 2, (double) _descs[ifrac].getFamily());
      segs.setValue(ifrac, 3, _descs[ifrac].getXXF(ip));
      segs.setValue(ifrac, 4, _descs[ifrac].getYYF(ip));
      segs.setValue(ifrac, 5, _descs[ifrac].getOrient());
      segs.setValue(ifrac, 6, (double) (ip == 0));
    }
  }

  return segs;
}

int FracList::_getEndPointCount() const
{
  int number = 0;
  for (int i = 0; i < getNFracs(); i++)
    number += _descs[i].getNPoint();
  return (number);
}

bool FracList::_isValidDisc(int idisc)
{
  if ((idisc) >= 0 && (idisc) < _ndisc) return true;
  return false;
}

/****************************************************************************/
/*!
 **  Import the Fractures
 **
 ** \return Pointer to the FracList structure
 **
 ** \param[in]  frac_segs    Array of fracture segments
 ** \param[in]  layinfo      Array of layer information
 ** \param[in]  nfamilies    Number of families
 **
 *****************************************************************************/
FracList* FracList::fractureImport(const VectorDouble& frac_segs,
                                   const VectorDouble& layinfo,
                                   int nfamilies)
{
  int nvalf = (int) frac_segs.size();
  if (!isMultiple(nvalf, NBYFRAC))
  {
    messerr("The number of values in 'frac_segs' (%d) should be a multiple of %d",
            nvalf, NBYFRAC);
    return nullptr;
  }
  int nseg = nvalf / NBYFRAC;

  int nvall = 0;
  int number = 0;
  int nlayer = 0;
  if (! layinfo.empty())
  {
    nvall = (int) layinfo.size();
    number = 1 + nfamilies * NPART;
    if (!isMultiple(nvall, number))
    {
      messerr("The number of values in 'layinfo' (%d) should be a multiple of %d",
              nvall, number);
      return nullptr;
    }
    nlayer = nvall / number;
  }

  /* Initializations */

  FracList* frac_list = new FracList();

  // Loop on Descriptions

  int icur = -1;
  for (int i = 0; i < nseg; i++)
  {
    int ifam = (int) FRACSEG(i, 2);
    double xx = FRACSEG(i, 3);
    double yy = FRACSEG(i, 4);
    double orient = FRACSEG(i, 5);
    bool flag_new = (int) FRACSEG(i, 6);

    if (flag_new || icur < 0)
    {
      frac_list->addDescription();
      icur = frac_list->getNFracs() - 1;
    }
    frac_list->setFamily(icur, ifam);
    frac_list->setOrient(icur, orient);
    frac_list->addPoint(icur, xx, yy);
  }

  // Loop on the layer informations

  if (! layinfo.empty())
  {
    frac_list->_layinfo.reset(nlayer, number);
    int ecr = 0;
    for (int i = 0; i < number; i++)
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
        frac_list->_layinfo.setValue(ilayer, i, layinfo[ecr++]);
  }
  return (frac_list);
}

void FracList::addDescription(const FracDesc& description)
{
  _descs.push_back(description);
}

/****************************************************************************/
/*!
 **  Plunge a subset of the simulated fractures into a block
 **
 **  \return Error return code
 **
 ** \param[in,out] dbgrid    Db structure
 ** \param[in]  xmax         Maximum extension along horizontal axis
 ** \param[in]  permtab      Permabilities per family (starting from 0)
 ** \param[in]  perm_mat     Permability for the matrix
 ** \param[in]  perm_bench   Permability along the bench edge
 ** \param[in]  ndisc        Number of discretization steps
 ** \param[in]  verbose      Verbose flag
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int FracList::fractureToBlock(DbGrid *dbgrid,
                              double xmax,
                              VectorDouble& permtab,
                              double perm_mat,
                              double perm_bench,
                              int ndisc,
                              bool verbose,
                              const NamingConvention& namconv)
{
  if (dbgrid->getNDim() != 2)
  {
    messerr("This application is limited to 2-D grid");
    return 1;
  }
  double dmin = 1.e30;
  for (int idim = 0; idim < dbgrid->getNDim(); idim++)
  {
    if (dbgrid->getDX(idim) < dmin) dmin = dbgrid->getDX(idim);
  }
  double delta = dmin / (double) ndisc;

  // Allocate the new variable

  int iptr = dbgrid->addColumnsByConstant(1, perm_mat);

  // Plunge the environment

  if (perm_bench > 0)
  {
    for (int ilayer = 0; ilayer < _nlayers; ilayer++)
    {
      double cote = _getMemLayer(ilayer);
      _plungeSegment(dbgrid, iptr, delta, perm_bench, 0., cote, xmax, cote);
    }
  }

  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    FracDesc &desc = _descs[ifrac];
    int ifam = desc.getFamily();
    double perm = 0.;
    if (ifam < (int) permtab.size()) perm = permtab[ifam];
    if (perm <= 0.) continue;

    /* Loop on the segments */

    int npoint = desc.getNPoint();

    if (verbose)
      message("Fracture %d/%d (Number of end-points = %d)\n", ifrac+1, getNFracs(), npoint);

    double permval = permtab[desc.getFamily()];
    for (int ip = 0; ip < npoint - 1; ip++)
    {
      int jp = ip + 1;
      _plungeSegment(dbgrid, iptr, delta, permval,
                     desc.getXXF(ip), desc.getYYF(ip),
                     desc.getXXF(jp), desc.getYYF(jp));
    }
  }

  namconv.setNamesAndLocators(dbgrid, iptr);
  return 0;
}

/****************************************************************************/
/*!
 **  Plunge a segment in a block Db
 **
 ** \param[in]  dbgrid    Db structure
 ** \param[in]  iptr      Pointer to the Perm variable
 ** \param[in]  delta     Increment
 ** \param[in]  value     Value assigned to the fracture segment
 ** \param[in]  x1        First coordinate of the first point
 ** \param[in]  y1        Second coordinate of the first point
 ** \param[in]  x2        First coordinate of the second point
 ** \param[in]  y2        Second coordinate of the second point
 **
 *****************************************************************************/
void FracList::_plungeSegment(DbGrid *dbgrid,
                              int iptr,
                              double delta,
                              double value,
                              double x1,
                              double y1,
                              double x2,
                              double y2)
{
  VectorDouble coor(2);
  double deltax = x2 - x1;
  double deltay = y2 - y1;
  double dist = sqrt(deltax * deltax + deltay * deltay);
  int number = (int) floor(dist / delta);

  // Loop on the discretization lags

  for (int i = 0; i <= number; i++)
  {
    coor[0] = x1 + deltax * i / number;
    coor[1] = y1 + deltay * i / number;

    int iech = dbgrid->coordinateToRank(coor);
    if (iech >= 0) dbgrid->updArray(iech, iptr, EOperator::MAX, value);
  }
}

/****************************************************************************/
/*!
 **  Plunge a well line in a set of fractures
 **
 **  \return Array of the intersections (should be checked against NULL)
 **
 ** \param[in]  nval         Number of well information
 ** \param[in]  well         Array giving the well trajectory
 ** \param[in]  xmax         Maximum extension along horizontal axis
 ** \param[in]  permtab      Permabilities per family (starting from 0) Optional
 **
 ** \param[out] nint         Number of intersections
 ** \param[out] ncol         Number of attributes per intersection
 **
 ** \remark Output array must be freed by the calling function
 **
 *****************************************************************************/
VectorDouble FracList::fractureToWell(int nval,
                                      const VectorDouble& well,
                                      double xmax,
                                      const VectorDouble& permtab,
                                      int *nint,
                                      int *ncol)
{
  double x, y;
  VectorDouble wellout;

  /* Preliminary checks */

  if (!isMultiple(nval, 2))
  {
    messerr("The number of values read from 'well' should be a multiple of 2");
    return VectorDouble();
  }
  int nw_xy = nval / 2;
  if (nw_xy < 2)
  {
    messerr("Number of end points for the well line must not be less than 2");
    return VectorDouble();
  }

  /* Loop on the line segments */

  for (int iw = 0; iw < nw_xy - 1; iw++)
  {
    double x1 = WELL(iw, 0);
    double y1 = WELL(iw, 1);
    double x2 = WELL(iw + 1, 0);
    double y2 = WELL(iw + 1, 1);

    /* Store the starting point */

    _welloutAdd(wellout, x1, y1, -1, -1, 0, 0.);

    /* Loop on the fractures */

    for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
    {
      FracDesc& desc = _descs[ifrac];
      int ifam = desc.getFamily();
      int npoint = desc.getNPoint();

      /* Loop on the segments */

      for (int ip = 0; ip < npoint - 1; ip++)
      {
        int jp = ip + 1;
        if (! GH::segmentIntersect(x1, y1, x2, y2,
                                   desc.getXXF(ip), desc.getYYF(ip),
                                   desc.getXXF(jp), desc.getYYF(jp),
                                   &x, &y)) continue;
        if (y <= -_eps || x <= -_eps || x >= xmax + _eps) continue;

        /* Store the new intersection */

        double perm = (permtab.empty()) ? 0. : permtab[ifam];
        _welloutAdd(wellout, x, y, ifrac, ip, ifam, perm);
      }
    }

    /* Store the ending point */

    _welloutAdd(wellout, x2, y2, -1, -1, 0, 0.);
  }

  (*nint) = (int) wellout.size() / NBYWOUT;
  (*ncol) = NBYWOUT;
  return wellout;
}

/****************************************************************************/
/*!
 **  Add an end point to the wellout array
 **
 ** \param[in]  wellout      Array provided as input
 ** \param[in]  x            Ascissae of the intersection
 ** \param[in]  y            Ordinate of the intersection
 ** \param[in]  ifrac        Rank of the fracture (starting from 1)
 ** \param[in]  ip           Rank of the segment (starting from 1)
 ** \param[in]  family       Family to which the fracture belongs
 ** \param[in]  perm         Assigned permeability or TEST
 **
 *****************************************************************************/
void FracList::_welloutAdd(VectorDouble& wellout,
                           double x,
                           double y,
                           int ifrac,
                           int ip,
                           int family,
                           double perm)
{
  int nloc;

  nloc = (int) wellout.size() / NBYWOUT;
  wellout.resize(NBYWOUT * (nloc+1));

  WELLOUT(nloc,0) = x; /* First coordinate */
  WELLOUT(nloc,1) = y; /* Second coordinate */
  WELLOUT(nloc,2) = (double) ifrac + 1; /* Rank of the fracture */
  WELLOUT(nloc,3) = (double) ip + 1; /* Rank of segment in frature */
  WELLOUT(nloc,4) = (double) family; /* Rank of family */
  WELLOUT(nloc,5) = perm; /* Starting permeability */
  WELLOUT(nloc,6) = perm; /* Ending permeability */
  WELLOUT(nloc,7) = TEST; /* Range of permeability change */
}

/****************************************************************************/
/*!
 **  Extract the array of fracture lengths
 **
 ** \return The returned array
 **
 ** \param[in]  ifam       Rank of the family or ITEST for all
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 *****************************************************************************/
VectorDouble FracList::fractureExtractLength(int ifam,
                                             double cote,
                                             double dcote)
{
  VectorDouble tab(getNFracs(), 0);

  /* Loading the fractures */

  int ecr = 0;
  for (int i = 0; i < getNFracs(); i++)
  {
    FracDesc& desc = _descs[i];

    /* Selection according to the family criterion */

    if (!IFFFF(ifam) && desc.getFamily() != ifam) continue;

    /* Calculate the fracture length */

    double value = desc.fractureExtension(cote, dcote);
    if (value <= 0.) continue;
    tab[ecr++] = value;
  }

  /* Sort the intersections */

  VH::sortInPlace(tab, true, ecr);

  return tab;
}

/****************************************************************************/
/*!
 **  Extract the fracture interdistances
 **
 ** \return The returned array
 **
 ** \param[in]  ifam         Rank of the family or ITEST for all
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 *****************************************************************************/
VectorDouble FracList::fractureExtractDist(int ifam, double cote, double dcote)
{
  int number = _getEndPointCount();
  VectorDouble tab(number, 0);

  /* Loading the abcissae of the fracture (with the target layer) */

  int ecr = 0;
  for (int i = 0; i < getNFracs(); i++)
  {
    FracDesc& desc = _descs[i];
    int npoint = desc.getNPoint();

    /* Loop on the end points */

    for (int ip = 0; ip < npoint; ip++)
    {

      /* Selection according to the family criterion */

      if (!IFFFF(ifam) && desc.getFamily() != ifam) continue;

      /* Selection according to the layer */

      if (!FFFF(cote) && !FFFF(dcote) &&
      ABS(cote - desc.getYYF(ip)) > dcote) continue;

      tab[ecr++] = desc.getXXF(ip);
    }
  }
  int ndist = ecr;
  if (ndist <= 0) return VectorDouble();

  /* Sort the intersections */

  VH::sortInPlace(tab, true, ndist);

  /* Calculate the inter-fracture distances */

  double valref = tab[0];
  for (int i = 1; i < ndist; i++)
  {
    double value = tab[i];
    tab[i - 1] = value - valref;
    valref = value;
  }

  /* Sort the fracture inter-distances */

  ndist--;
  VH::sortInPlace(tab, true, ndist);

  return tab;
}

/****************************************************************************/
/*!
 **  Plunge a line trajectory and the modified permeabilities within an
 **  existing block
 **
 **  \return Error return code
 **
 ** \param[in,out] dbgrid    Db structure
 ** \param[in]  col_perm     Existing attribute for permeability (or ITEST)
 ** \param[in]  col_fluid    Existing attribute for fluid (or ITEST)
 ** \param[in]  flag_fluid   1 for performing the Fluid filling
 ** \param[in]  val_fluid    Value assigned to the fluid
 ** \param[in]  wellout      Pointer to the new wellout information
 ** \param[in]  nval         Number of values for wellout informations
 ** \param[in]  ndisc        Number of discretization steps
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
int FracList::fractureWellToBlock(DbGrid *dbgrid,
                                  int col_perm,
                                  int col_fluid,
                                  int flag_fluid,
                                  double val_fluid,
                                  const VectorDouble& wellout,
                                  int nval,
                                  int ndisc,
                                  bool verbose)
{
  VectorDouble traj;
  int iptr_perm = -1;
  int iptr_fluid = -1;

  /* Preliminary checks */

  if (!isMultiple(nval, NBYWOUT))
  {
    messerr("The number of values read (%d) in well should be a multiple of %d",
            nval, NBYWOUT);
    return 1;
  }
  int nw_xy = nval / NBYWOUT;
  if (dbgrid->getNDim() != 2)
  {
    messerr("This application is limited to 2-D grid");
    return 1;
  }
  double dmin = 1.e30;
  for (int idim = 0; idim < dbgrid->getNDim(); idim++)
  {
    if (dbgrid->getDX(idim) < dmin) dmin = dbgrid->getDX(idim);
  }
  double delta = dmin / (double) ndisc;

  /* Allocate the new variable */

  iptr_perm = dbgrid->addColumnsByConstant(1, 0);
  if (!IFFFF(col_perm)) db_attribute_copy(dbgrid, col_perm, iptr_perm);

  if (flag_fluid)
  {
    iptr_fluid = dbgrid->addColumnsByConstant(1, 0.);
    if (!IFFFF(col_fluid)) db_attribute_copy(dbgrid, col_fluid, iptr_fluid);
  }

  /* Verbose option */

  if (verbose)
    print_matrix("Well information", 0, 0, NBYWOUT, nw_xy, NULL, wellout.data());

  /* Paint the fluid (optional) */

  if (flag_fluid)
  {
    for (int iw = 0; iw < nw_xy - 1; iw++)
    {
      double x1 = WELLOUT(iw, 0);
      double y1 = WELLOUT(iw, 1);
      double x2 = WELLOUT(iw + 1, 0);
      double y2 = WELLOUT(iw + 1, 1);
      _plungeSegment(dbgrid, iptr_fluid, delta, val_fluid, x1, y1, x2, y2);
    }
  }

  /* Loop on the end points of the wellout */

  for (int iw = 0; iw < nw_xy; iw++)
  {
    double xx = WELLOUT(iw, 0);
    double yy = WELLOUT(iw, 1);
    int ifrac = (int) WELLOUT(iw, 2) - 1;
    int ip = (int) WELLOUT(iw, 3) - 1;
    double perm1 = WELLOUT(iw, 5);
    double perm2 = WELLOUT(iw, 6);
    double range = WELLOUT(iw, 7);
    if (ifrac < 0) continue;

    if (ifrac >= getNFracs())
    {
      messerr("The well information (line %d/%d) is invalid:", iw + 1, nw_xy);
      messerr("- Fracture number (%d) should lie within [1,%d]", ifrac + 1,
              getNFracs());
      continue;
    }

    FracDesc &desc = _descs[ifrac];
    if (ip < 0 || ip >= desc.getNPoint())
    {
      messerr("The well information (line %d/%d) is invalid:", iw + 1, nw_xy);
      messerr("- For fracture (%d, Segment number (%d) is not within [1,%d]",
              ifrac + 1, ip + 1, desc.getNPoint());
      continue;
    }
    int npoint = desc.getNPoint();

    /* Paint the blocks with the stationary permeability */

    for (int i = 0; i < npoint - 1; i++)
      _plungeSegment(dbgrid, iptr_perm, delta, perm2,
                     desc.getXXF(i), desc.getYYF(i),
                     desc.getXXF(i + 1), desc.getYYF(i + 1));

    /* Allocate the trajectory */

    traj.resize(2 * (npoint + 1), 0);

    /* Trajectory upwards */

    _trajAdd(traj, xx, yy);
    for (int i = ip; i >= 0; i--)
      _trajAdd(traj, desc.getXXF(ip), desc.getYYF(ip));
    _plungeSegmentGradual(dbgrid, iptr_perm, delta, traj, perm1, perm2, range);

    /* Trajectory downwards */

    _trajAdd(traj, xx, yy);
    for (int i = ip; i < npoint; i++)
      _trajAdd(traj, desc.getXXF(ip), desc.getYYF(ip));
    _plungeSegmentGradual(dbgrid, iptr_perm, delta, traj, perm1, perm2, range);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Add a sample to the trajectory
 **
 ** \param[in]  traj         Array for storing the trajectory
 ** \param[in]  x            First coordinate
 ** \param[in]  y            Second coordinate
 **
 *****************************************************************************/
void FracList::_trajAdd(VectorDouble& traj, double x, double y)
{
  int nloc = (int) traj.size();
  traj.resize(2 * (nloc+1));
  TRAJ(nloc,0) = x;
  TRAJ(nloc,1) = y;
}

/****************************************************************************/
/*!
 **  Plunge a segment painted with gradually changing permeability in a block Db
 **
 ** \param[in]  dbgrid    Db structure
 ** \param[in]  iptr      Pointer to the Perm variable
 ** \param[in]  delta     Increment
 ** \param[in]  traj      Vector describing the trajectory
 ** \param[in]  perm1     Permeability at origin of curvilinear distance
 ** \param[in]  perm2     Permeability at range of curvilinear distance
 ** \param[in]  range     Range of the permeability change
 **
 *****************************************************************************/
void FracList::_plungeSegmentGradual(DbGrid *dbgrid,
                                     int iptr,
                                     double delta,
                                     VectorDouble& traj,
                                     double perm1,
                                     double perm2,
                                     double range)
{
  VectorDouble coor(2);

  int ntraj = (int) traj.size();
  double total = 0.;
  for (int ip = 0; ip < ntraj - 1; ip++)
  {
    double x1 = TRAJ(ip, 0);
    double y1 = TRAJ(ip, 1);
    double x2 = TRAJ(ip + 1, 0);
    double y2 = TRAJ(ip + 1, 1);
    double deltax = x2 - x1;
    double deltay = y2 - y1;
    double dist = sqrt(deltax * deltax + deltay * deltay);
    int number = (int) floor(dist / delta);
    double incx = deltax / number;
    double incy = deltay / number;
    double incr = dist / number;

    /* Loop on the discretization lags */

    for (int i = 0; i <= number; i++)
    {
      coor[0] = x1 + i * incx;
      coor[1] = y1 + i * incy;
      total += incr;

      int iech = dbgrid->coordinateToRank(coor);
      if (iech < 0) continue;

      double perm;
      if (FFFF(range))
        perm = perm2;
      else
      {
        double h = total / range;
        double red = _cubic(h);
        perm = perm2 + (perm1 - perm2) * red;
      }
      dbgrid->setArray(iech, iptr, perm);
    }
  }
}

