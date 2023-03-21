/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include <Geometry/GeometryHelper.hpp>
#include "geoslib_old_f.h"

#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Space/ASpace.hpp"

#include <math.h>

DirParam::DirParam(int npas,
                   double dpas,
                   double toldis,
                   double tolang,
                   int opt_code,
                   int idate,
                   double bench,
                   double cylrad,
                   double tolcode,
                   const VectorDouble& breaks,
                   const VectorDouble& codir,
                   const VectorInt& grincr,
                   const ASpace* space)
    : ASpaceObject(space),
      _nPas(npas),
      _optionCode(opt_code),
      _idate(idate),
      _definedForGrid(false),
      _dPas(dpas),
      _bench(bench),
      _cylRad(cylrad),
      _tolDist(toldis),
      _tolAngle(tolang),
      _tolCode(tolcode),
      _breaks(breaks),
      _codir(codir),
      _grincr(grincr)
{
  _completeDefinition();
}

DirParam::DirParam(const DirParam& r)
    : ASpaceObject(r),
      _nPas(r._nPas),
      _optionCode(r._optionCode),
      _idate(r._idate),
      _definedForGrid(r._definedForGrid),
      _dPas(r._dPas),
      _bench(r._bench),
      _cylRad(r._cylRad),
      _tolDist(r._tolDist),
      _tolAngle(r._tolAngle),
      _tolCode(r._tolCode),
      _breaks(r._breaks),
      _codir(r._codir),
      _grincr(r._grincr)
{
}

DirParam& DirParam::operator=(const DirParam& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _nPas = r._nPas;
    _optionCode = r._optionCode;
    _idate = r._idate;
    _definedForGrid = r._definedForGrid;
    _dPas = r._dPas;
    _bench = r._bench;
    _cylRad = r._cylRad;
    _tolDist = r._tolDist;
    _tolAngle = r._tolAngle;
    _tolCode = r._tolCode;
    _breaks = r._breaks;
    _codir = r._codir;
    _grincr = r._grincr;
  }
  return *this;
}

DirParam::~DirParam()
{
}

DirParam* DirParam::create(int npas,
                           double dpas,
                           double toldis,
                           double tolang,
                           int opt_code,
                           int idate,
                           double bench,
                           double cylrad,
                           double tolcode,
                           const VectorDouble& breaks,
                           const VectorDouble& codir,
                           const ASpace* space)
{
  return new DirParam(npas, dpas, toldis, tolang, opt_code, idate,
                      bench, cylrad, tolcode, breaks, codir, VectorInt(), space);
}

DirParam* DirParam::createOmniDirection(int npas,
                                        double dpas,
                                        double toldis,
                                        int opt_code,
                                        int idate,
                                        double bench,
                                        double cylrad,
                                        double tolcode,
                                        const VectorDouble& breaks,
                                        const ASpace* space)
{
  return new DirParam(npas, dpas, toldis, 90.1, opt_code, idate,
                      bench, cylrad, tolcode, breaks, VectorDouble(), VectorInt(), space);
}

DirParam* DirParam::createFromGrid(int npas,
                                   const VectorInt &grincr,
                                   const ASpace *space)
{
  return new DirParam(npas, 0., 0.5, 90., 0, 0, TEST, TEST, 0.,
                      VectorDouble(), VectorDouble(), grincr, space);
}

double DirParam::getBreak(int i) const
{
  if (i < 0 || i >= (int)_breaks.size())
  {
    mesArg("Break Index",i,(int) _breaks.size());
    return TEST;
  }
  return _breaks[i];
}

double DirParam::getCodir(int i) const
{
  if (i < 0 || i >= (int)_codir.size())
  {
    mesArg("Codir Index",i,(int) _codir.size());
    return TEST;
  }
  return _codir[i];
}

void DirParam::_completeDefinition()
{
  if (! _breaks.empty())
  {
    if (_breaks.size() < 2) _breaks.clear();
  }

  int ndim = getNDim();

  bool flagPoint = true;
  if (_codir.empty())
  {
    flagPoint = false;
    _codir.resize(ndim,0.);
    _codir[0] = 1.;
  }

  bool flagGrid = true;
  if (_grincr.empty())
  {
    flagGrid = false;
    _grincr.resize(ndim,0);
    _grincr[0] = 1;
  }

  if (flagPoint)
    _definedForGrid = false;
  if (flagGrid)
    _definedForGrid = true;
  // Correction for the particular case of Omni-Direction definition
  if (! flagPoint && ! flagGrid) _definedForGrid = false;
}

bool DirParam::isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= (int) getNDim())
  {
    mesArg("Space Dimension",idim,getNDim());
    return false;
  }
  return true;
}

bool DirParam::isLagValid(int ilag, bool flagAsym) const
{
  int nlag = getLagNumber();
  if (flagAsym) nlag = 2 * nlag + 1;
  if (ilag < 0 || ilag >= nlag)
  {
    mesArg("Lag Index",ilag,nlag);
    return false;
  }
  return true;
}

/**
 * Set the value of the lag as computed from the Db (Grid organized)
 * @param db Db structure
 */
void DirParam::setDPas(const DbGrid* db)
{
  if (_grincr.empty()) return;
  double dpas = 0;
  for (int idim = 0; idim < (int) getNDim(); idim++)
  {
    double delta = _grincr[idim] * db->getDX(idim);
    dpas += delta * delta;
  }
  _dPas = sqrt(dpas);
}

int DirParam::getGrincr(int idim) const
{
  if (_grincr.empty()) return 0;
  if (! isDimensionValid(idim)) return 0;
  return _grincr[idim];
}

double DirParam::getMaximumDistance() const
{
  double maxdist;

  if (getFlagRegular())
    maxdist = getDPas() * (getLagNumber() + getTolDist());
  else
    maxdist = getBreak(getLagNumber());
  return (maxdist);
}

bool DirParam::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String DirParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ndim = getNDim();

  if (getLagNumber() > 0)
    sstr << "Number of lags              = " << getLagNumber() << std::endl;

  if (! _definedForGrid)
  {

    // Case of a Direction defined for a non-grid Db

    sstr << toVector("Direction coefficients      = ", _codir);
    if (ndim > 1)
    {
      VectorDouble angles(ndim);
      (void) GH::rotationGetAngles(_codir,angles);
      sstr << toVector("Direction angles (degrees)  = ", angles);
    }

    if (! FFFF(_tolAngle))
      sstr << "Tolerance on direction      = " << toDouble(_tolAngle)
      << " (degrees)" << std::endl;
  }

  if (! FFFF(_bench)  && _bench > 0.)
    sstr << "Slice bench                 = " << toDouble(_bench) << std::endl;
  if (! FFFF(_cylRad) && _cylRad > 0.)
    sstr << "Slice radius                = " << toDouble(_cylRad) << std::endl;

  if (getFlagRegular())
  {
    if (getDPas() > .0)
    {
      sstr << "Calculation lag             = " << toDouble(getDPas()) << std::endl;
      sstr << "Tolerance on distance       = " << toDouble(100. * getTolDist())
                     << " (Percent of the lag value)" << std::endl;
    }
  }
  else
  {
    sstr << "Calculation intervals       = " << std::endl;
    for (int i = 0; i < getBreakNumber(); i++)
    {
      sstr << " - Interval " << i + 1 << " = ["
          << toInterval(getBreak(i), getBreak(i + 1)) << "]" << std::endl;
    }
  }

  if (_definedForGrid)
  {

    // Case of a variogram defined on a Grid db

    sstr << toVector("Grid Direction coefficients = ", _grincr);
  }

  /* Selection on the 'code' */

  if (getOptionCode() == 1)
    sstr << "Selection if Codes are close enough (" << getTolCode() << ")"
         << std::endl;
  if (getOptionCode() == 2)
    sstr << "Selection if Codes are different" << std::endl;

  return sstr.str();
}

std::vector<DirParam> DirParam::createMultiple(int ndir,
                                               int npas,
                                               double dpas,
                                               double toldis,
                                               const ASpace* space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VectorDouble angles = VectorDouble(1);
  VectorDouble codir  = VectorDouble(ndim);
  std::vector<DirParam> dirs;
  for (int idir = 0; idir < ndir; idir++)
  {
    angles[0] = 180. * (double) idir / (double) ndir;
    (void) GH::rotationGetDirection(ndim, 1, angles,codir);
    double tolang = 90. / (double) ndir;
    DirParam dirparam = DirParam(npas, dpas, toldis, tolang, 0, 0, TEST, TEST, 0.,
                                 VectorDouble(), codir, VectorInt(), space);
    dirs.push_back(dirparam);
  }
  return dirs;
}

std::vector<DirParam> DirParam::createMultipleFromGrid(int npas, const ASpace* space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VectorInt grincr = VectorInt(ndim);
  std::vector<DirParam> dirs;
  int ndir = ndim;
  for (int idir = 0; idir < ndir; idir++)
  {
    VH::fill(grincr, 0);
    grincr[idir] = 1;
    DirParam* dirparam = DirParam::createFromGrid(npas, grincr, space);
    dirs.push_back(*dirparam);
  }
  return dirs;
}
