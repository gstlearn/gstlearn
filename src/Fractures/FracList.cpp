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
#include "Fractures/FracList.hpp"
#include "Fractures/Description.hpp"
#include "Fractures/Environ.hpp"
#include "Fractures/Family.hpp"
#include "Fractures/Fault.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Geometry.hpp"
#include "Basic/Law.hpp"
#include "Matrix/MatrixRectangular.hpp"

#include <math.h>

FracList::FracList(int ndisc, bool flag_check, double low0, double low1, double eps)
  : AStringable(),
    _descs(),
    _layinfo(),
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
 ** \param[in]  environ        Environ structure
 ** \param[in]  flag_sim_layer TRUE for simulating layers
 **                            FALSE if they are read
 ** \param[in]  flag_sim_fract TRUE for simulating the fractures
 **                            FALSE establish the layers and the main faults
 ** \param[in]  nlayers_in     Number of input layers (used if flag_sim_layer=F)
 ** \param[in]  elevations     Array of elevations (used if flag_sim_layer=F)
 ** \param[in]  seed           Seed for the random number generator
 ** \param[in]  verbose        Verbose option
 **
 *****************************************************************************/
int FracList::simulate(const Environ* environ,
                       bool flag_sim_layer,
                       bool flag_sim_fract,
                       int seed,
                       bool verbose,
                       int nlayers_in,
                       const VectorDouble& elevations)
{
  _step = environ->getXextend() / _ndisc;
  _xorigin = -environ->getDeltax();
  _verbose = verbose;

  double y0 = 0.;
  int ninfos = 1 + environ->getNFamilies() * NPART;
  law_set_random_seed(seed);

  if (_verbose)
  {
    message("Options set by the keypair mechanism:\n");
    message("Fracture_Discretization_Count = %d \n", _ndisc);
    message("Fracture_Check_Intersect      = %d \n", _flagCheck);
    message("Fracture_Repulsion_Low0       = %lg\n", _low0);
    message("Fracture_Repulsion_Low1       = %lg\n", _low1);
  }

  /* Determine the layers */

  VectorDouble thicks;
  if (flag_sim_layer)
    thicks = _layersManage(*environ, &y0);
  else
    thicks = _layersRead(nlayers_in, elevations, &y0);
  int nlayers = (int) thicks.size();

  /* Core allocation */

  int np1 = nlayers + 1;
  VectorDouble denstab(_ndisc);
  _layinfo = MatrixRectangular(np1, ninfos);

  /* Define the layers */

  double cote = y0;
  for (int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    double thick = thicks[ilayer];
    setMemLayer(ilayer, cote);
    cote += thick;
  }
  setMemLayer(nlayers, cote);

  /* Define the main faults */

  for (int ifault = 0; ifault < environ->getNFaults(); ifault++)
  {
    double angle = environ->getFault(ifault).getOrient();

    /* Loop on the layers */

    int ifrac = -1;
    cote = y0;
    double xx = environ->getFault(ifault).faultAbscissae(cote);
    for (int ilayer = 0; ilayer < nlayers; ilayer++)
    {
      double thick = thicks[ilayer];
      ifrac = _fracAdd(ifrac, 0, xx, cote, thick, angle, &xx);
      cote += thick;
    }
  }

  if (_verbose)
    message("Total number of main faults        = %d \n", getNFracs());

  if (!flag_sim_fract) goto label_suite;

  /* Loop on the different fault families */

  for (int ifam = 0; ifam < environ->getNFamilies(); ifam++)
  {
    if (_verbose)
      mestitle(0, "Processing Family #%d/%d", ifam + 1, environ->getNFamilies());
    const Family& family = environ->getFamily(ifam);

    /* Loop on the layers */

    cote = y0;
    double thetap = 0.;
    for (int ilayer = 0; ilayer < nlayers; ilayer++)
    {
      double thick = thicks[ilayer];
      if (_verbose)
      {
        mestitle(1, "Processing Layer #%d/%d", ilayer + 1, nlayers);
        message("Elevation of the layer bottom     = %lf\n", cote);
        message("Thickness of the layer            = %lf\n", thick);
      }

      /* Derive the layer intensity */

      double theta1 = _layerIntensity(family, thick);

      /* Generate the density function */

      _generateDensity(*environ, family, ifam, cote, denstab);

      /* Correct the density due to already existing faults */

      _correctDensity(family, ifam, cote, denstab);

      /* Extend previously existing fractures */

      double propsur = _extendFractures(family, ifam, cote, thick, denstab);

      /* Derive the layer intensity, given the survival from previous layers */

      double theta2 = _deriveIntensity(theta1, thetap, propsur);

      /* Simulate the fractures abscissa */

      int nfracs = _simulateFractures(*environ, family, ifam, cote, thick, theta2, denstab);

      /* Shift the ordinate */

      cote += thick;
      thetap = theta1;
      setMemTheta1(ilayer,ifam,theta1);
      setMemTheta2(ilayer,ifam,theta2);
      setMemPropsur(ilayer,ifam,propsur);
      setMemFrac(ilayer,ifam,(double) nfracs);
      setMemTotal(ilayer,ifam,(double) getNFracs());
    }
  }

  label_suite:
  return (0);
}

/****************************************************************************/
/*!
 **  Generate the series of layers
 **
 ** \return The vector of layer thicknesses
 **
 ** \param[in]  environ      Pointer to the Frac_Environ structure
 ** \param[in]  y0           Ordinate of the origin
 **
 ** \remarks The origin is calculated so that a layer edge is located at 0.
 **
 *****************************************************************************/
VectorDouble FracList::_layersManage(const Environ& environ,
                                     double *y0)
{
  VectorDouble thicks;
  double thick_min = environ.getMean() / 10.;
  double mean = environ.getMean();
  double stdev = environ.getStdev();
  double var = stdev * stdev;
  double theta = var / mean;
  bool flag_rnd = var > 0.;
  int k = (flag_rnd) ? (mean * mean) / var : 0.;

  /* Allocation */

  int number = 0;

  /* Downwards */

  double total = 0.;
  while (total < environ.getDeltay())
  {
    double thick = (flag_rnd) ? theta * law_gamma(k) : mean;
    if (thick < thick_min) continue;
    total += thick;
    thicks.push_back(thick);
  }
  double cote_down = -total;

  /* Upwards */

  total = 0.;
  while (total < environ.getYmax())
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
 ** \param[in]  nlayers_in   Number of elevations to be read
 ** \param[in]  elevations   Array of elevations to be read
 **
 ** \param[out] y0           Ordinate of the origin
 **
 *****************************************************************************/
VectorDouble FracList::_layersRead(int nlayers_in,
                                   const VectorDouble& elevations,
                                   double *y0)
{
  VectorDouble thicks;

  /* Initializations */

  thicks.resize(nlayers_in - 1);

  double cote = elevations[0];
  for (int i = 1; i < nlayers_in; i++)
  {
    thicks[i - 1] = elevations[i] - cote;
    cote = elevations[i];
  }
  double cote_down = (*y0) = elevations[0];
  double cote_up = elevations[nlayers_in - 1];

  /* Optional information */

  if (_verbose)
  {
    mestitle(0, "Layer (read)");
    message("Minimum simulated level            = %lf\n", cote_down);
    message("Maximum simulated level            = %lf\n", cote_up);
    message("Number of layers                   = %d \n", nlayers_in);
  }

  return thicks;
}

/****************************************************************************/
/*!
 **  Manage the Frac_List structure
 **
 ** \return  Pointer to the Frac_List structure
 **
 ** \param[in]  mode         0 initialization; 1 allocation; -1 deallocation
 **
 *****************************************************************************/
void FracList::_manage(int mode)
{
  int number = getNFracs();
  switch (mode)
  {
    case 0:
      _descs.clear();
      break;

    case 1:

      _descs.resize(number + 1);
      _descs[number] = Description();
      break;

    case -1:
      _descs.clear();
      break;
  }
  return;
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
    _manage(1);
    ifrac = getNFracs() - 1;
  }

  /* Set the attributes of the new segment */

  Description& desc = _descs[ifrac];

  if (desc.getNPoint() == 0)
  {
    desc.setFamily(ifam);
    desc.setOrient(orient);
    desc.addPoint(xx, cote);
  }
  desc.setFamily(ifam);
  desc.setOrient(orient);
  desc.addPoint(thick * tan(ut_deg2rad(orient)) + xx, thick + cote);
  *xp = desc.getXXF(desc.getNPoint() - 1);

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
  Description& desc = _descs[ifrac0];
  int npoint = desc.getNPoint();
  double xd1 = desc.getXXF(npoint - 2);
  double yd1 = desc.getYYF(npoint - 2);
  double xe1 = desc.getXXF(npoint - 1);
  double ye1 = desc.getYYF(npoint - 1);

  /* Check against the previous fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    if (ifrac == ifrac0) continue;
    Description& desc2 = _descs[ifrac];

    /* Look for the segment belonging to the current layer */

    if (! _belongToLayer(desc2, cote, &xd2, &yd2, &xe2, &ye2)) continue;

    /* Perform the intersection between two segments */

    if (segment_intersect(xd1, yd1, xe1, ye1, xd2, yd2, xe2, ye2, &x, &y))
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
 ** \param[in]  desc         Frac_Desc structure
 ** \param[in]  cote         Ordinate of the fracture starting point
 **
 ** \param[in]  xd           Abscissae of the first end-point
 ** \param[in]  yd           Ordinate of the first end-point
 ** \param[in]  xe           Abscissae of the second end-point
 ** \param[in]  ye           Ordinate of the second end-point
 **
 *****************************************************************************/
bool FracList::_belongToLayer(const Description& desc,
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
double FracList::_layerIntensity(const Family& family, double thick)
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
 ** \param[in]  environ     Environ structure
 ** \param[in]  family      Family structure
 ** \param[in]  ifam        Rank of the family
 ** \param[in]  cote        Ordinate of the fracture starting point
 **
 ** \param[out] denstab      Discretization density array
 **
 *****************************************************************************/
void FracList::_generateDensity(const Environ& environ,
                                const Family& family,
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
      for (int ifault = 0; ifault < environ.getNFaults(); ifault++)
      {
        const Fault& fault = environ.getFault(ifault);
        if (! _sameFaultSide(environ, ifault, x0)) continue;

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
  {
    double total;
    (void) _densityCumulate("Main Fault         ", true, denstab, &total);
  }

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
void FracList::_correctDensity(const Family& family,
                               int ifam,
                               double cote,
                               VectorDouble& denstab)
{
  double total;

  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    Description& desc = _descs[ifrac];
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
    (void) _densityCumulate("Previous families  ", true, denstab, &total);
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
double FracList::_extendFractures(const Family& family,
                                  int ifam,
                                  double cote,
                                  double thick,
                                  VectorDouble& denstab)
{
  int next = 0;
  int ntotal = 0;
  for (int ifrac = 0; ifrac < getNFracs(); ifrac++)
  {
    const Description& desc = _descs[ifrac];
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
    (void) _fracAdd(ifrac, ifam, xx, cote, thick, angle, &xp);

    /* Update the density to account for repulsion */

    _updateRepulsion(xx, family.getRange(), denstab);
  }

  double propsur = (ntotal > 0) ? (double) next / (double) ntotal : 0.;
  if (_verbose && ntotal > 0)
  {
    message("Survival proportion               = %lf (%d from %d)\n", propsur,
            next, ntotal);
    double total;
    (void) _densityCumulate("After survival     ", true, denstab, &total);
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
 ** \param[in]  environ      Environ structure
 ** \param[in]  ifault0      Rank of the target fault
 ** \param[in]  x0           Location of the discretized point
 **
 *****************************************************************************/
bool FracList::_sameFaultSide(const Environ& environ,
                              int ifault0,
                              double x0)
{
  double x1, x2;
  int ifault;

  x1 = environ.getFault(ifault0).getCoord();
  for (ifault = 0; ifault < environ.getNFaults(); ifault++)
  {
    if (ifault == ifault0) continue;
    x2 = environ.getFault(ifault).getCoord();

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
double FracList::_densityUpdate(const Fault& fault,
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
 ** \returns 1 if there is no more room for fractures
 **
 ** \param[in]  title        Title of the printout
 ** \param[in]  flag_print   1 if the message must be printed
 ** \param[in]  denstab      Discretized density array
 **
 ** \param[out] totdens      Cumulated density
 **
 *****************************************************************************/
bool FracList::_densityCumulate(const char *title,
                               bool flag_print,
                               const VectorDouble& denstab,
                               double *totdens)
{
  double total = 0.;
  for (int idisc = 0; idisc < _ndisc; idisc++)
    total += denstab[idisc];

  (*totdens) = total;
  bool flag_stop = (total <= _low1 * _ndisc);

  if (_verbose)
  {
    if (flag_stop)
      message("Fracture simulation interrupted: no more room available\n");
    else
    {
      if (flag_print)
        message("- Cumulated Distribution: %s = %lf\n", title, total);
    }
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
bool FracList::_fractureInterrupt(const Family& family,
                                  const Description& desc,
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
 ** \param[in]  fault  Frac_Fault structure
 ** \param[in]  cote   Ordinate of the fracture starting point
 **
 *****************************************************************************/
double FracList::_faultAbscissae(const Fault& fault, double cote)
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
double FracList::_fractureExtension(const Description& desc,
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
 ** \param[in]  environ      Environ structure
 ** \param[in]  family       Family structure
 ** \param[in]  ifam         Rank of the family
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  thick        Thickness of the layer
 ** \param[in]  theta        Intensity of the layer
 ** \param[in]  denstab      Discretized density array
 **
 *******************************************************************F**********/
int FracList::_simulateFractures(const Environ& environ,
                                 const Family& family,
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
  int nfracs = law_poisson(theta * environ.getXextend());

  /* Loop on the fractures to be simulated */

  int neff = 0;
  for (int ifrac = 0; ifrac < nfracs; ifrac++)
  {

    /* Calculate cumulative density */

    if (_densityCumulate("Fracture generation", false, denstab, &total)) break;

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
  {
    message("New Simulated fractures in Bench  = %d \n", neff);
    message("Total count of fractures          = %d \n", getNFracs());
  }
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
  int nbyfrac = 7;
  int ntotal = _getEndPointCount();
  MatrixRectangular segs(ntotal, nbyfrac);

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

