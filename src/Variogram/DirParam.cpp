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
#include <Geometry/GeometryHelper.hpp>

#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Space/ASpace.hpp"

#include <math.h>

DirParam::DirParam(int nlag,
                   double dlag,
                   double toldis,
                   double tolang,
                   int opt_code,
                   int idate,
                   double bench,
                   double cylrad,
                   double tolcode,
                   const VectorDouble& breaks,
                   const VectorDouble& codir,
                   double angle2D,
                   const ASpaceSharedPtr& space)
    : ASpaceObject(space),
      _nLag(nlag),
      _optionCode(opt_code),
      _idate(idate),
      _dLag(dlag),
      _bench(bench),
      _cylRad(cylrad),
      _tolDist(toldis),
      _tolAngle(tolang),
      _tolCode(tolcode),
      _breaks(breaks),
      _codir(codir),
      _grincr()
{
  _completeDefinition(angle2D);
}

DirParam::DirParam(const DbGrid *dbgrid,
                   int nlag,
                   const VectorInt &grincr,
                   const ASpaceSharedPtr& space)
    : ASpaceObject(space),
      _nLag(nlag),
      _optionCode(0),
      _idate(0),
      _dLag(0),
      _bench(TEST),
      _cylRad(TEST),
      _tolDist(0.5),
      _tolAngle(0.),
      _tolCode(0),
      _breaks(),
      _codir(),
      _grincr(grincr)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  _codir = dbgrid->getCodir(grincr);
  double dlag = 0.;
  for (int idim = 0; idim < ndim; idim++)
    dlag += _codir[idim] * _codir[idim];
  _dLag = sqrt(dlag);
  VH::normalize(_codir);
}

DirParam::DirParam(const DirParam& r)
    : ASpaceObject(r),
      _nLag(r._nLag),
      _optionCode(r._optionCode),
      _idate(r._idate),
      _dLag(r._dLag),
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
    _nLag = r._nLag;
    _optionCode = r._optionCode;
    _idate = r._idate;
    _dLag = r._dLag;
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

DirParam* DirParam::create(int nlag,
                           double dlag,
                           double toldis,
                           double tolang,
                           int opt_code,
                           int idate,
                           double bench,
                           double cylrad,
                           double tolcode,
                           const VectorDouble& breaks,
                           const VectorDouble& codir,
                           double angle2D,
                           const ASpaceSharedPtr& space)
{
  return new DirParam(nlag, dlag, toldis, tolang, opt_code, idate,
                      bench, cylrad, tolcode, breaks, codir, angle2D, space);
}

DirParam* DirParam::createOmniDirection(int nlag,
                                        double dlag,
                                        double toldis,
                                        int opt_code,
                                        int idate,
                                        double bench,
                                        double cylrad,
                                        double tolcode,
                                        const VectorDouble& breaks,
                                        const ASpaceSharedPtr& space)
{
  return new DirParam(nlag, dlag, toldis, 90.1, opt_code, idate,
                      bench, cylrad, tolcode, breaks, VectorDouble(), TEST, space);
}

DirParam* DirParam::createFromGrid(const DbGrid* dbgrid,
                                   int nlag,
                                   const VectorInt &grincr,
                                   const ASpaceSharedPtr& space)
{
  return new DirParam(dbgrid, nlag, grincr, space);
}

double DirParam::getBreak(int i) const
{
  if (!checkArg("Break Index", i, (int)_breaks.size())) return TEST;
  return _breaks[i];
}

double DirParam::getCodir(int i) const
{
  if (!checkArg("Codir Index", i, (int)_codir.size())) return TEST;
  return _codir[i];
}

void DirParam::_completeDefinition(double angle2D)
{
  if (! _breaks.empty())
  {
    if (_breaks.size() < 2) _breaks.clear();
  }

  int ndim = getNDim();

  if (! FFFF(angle2D))
  {
    _codir.resize(ndim,0.);
    _codir[0] = cos(angle2D * GV_PI / 180.);
    _codir[1] = sin(angle2D * GV_PI / 180.);
  }

  if (_codir.empty())
  {
    _codir.resize(ndim, 0.);
    _codir[0] = 1.;
  }

  // Capping the tolerance on angles
  if (_tolAngle > 90.) _tolAngle = 90.;
}

void DirParam::setTolAngle(double tolang)
{
  _tolAngle = tolang;
  if (_tolAngle > 90.) _tolAngle = 90.;
}

bool DirParam::isDimensionValid(int idim) const
{
  return checkArg("Space Dimension", idim, getNDim());
}

bool DirParam::isLagValid(int ilag, bool flagAsym, bool flagCheck) const
{
  if (!flagCheck) return true;
  int nlag = getNLag();
  if (flagAsym) nlag = 2 * nlag + 1;
  return checkArg("Lag Index", ilag, nlag);
}

/**
 * Set the value of the lag as computed from the Db (Grid organized)
 * @param db Db structure
 */
void DirParam::setDPas(const DbGrid* db)
{
  if (_grincr.empty()) return;
  double dlag = 0;
  for (int idim = 0; idim < (int) getNDim(); idim++)
  {
    double delta = _grincr[idim] * db->getDX(idim);
    dlag += delta * delta;
  }
  _dLag = sqrt(dlag);
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
    maxdist = getDPas() * (getNLag() + getTolDist());
  else
    maxdist = getBreak(getNLag());
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

  if (getNLag() > 0)
    sstr << "Number of lags              = " << getNLag() << std::endl;

  sstr << toVector("Direction coefficients      = ", _codir);
  if (ndim > 1)
  {
    VectorDouble angles(ndim);
    (void) GH::rotationGetAnglesFromCodirInPlace(_codir,angles);
    if (ndim > 2)
      sstr << toVector("Direction angles (degrees)  = ", angles);
    else
      sstr << "Direction angles (degrees)  = " << toDouble(angles[0]) << std::endl;
  }

  if (! FFFF(_tolAngle))
    sstr << "Tolerance on direction      = " << toDouble(_tolAngle)
    << " (degrees)" << std::endl;

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
    for (int i = 0; i < getNBreak(); i++)
    {
      sstr << " - Interval " << i + 1 << " = ["
          << toInterval(getBreak(i), getBreak(i + 1)) << "]" << std::endl;
    }
  }

  if (! _grincr.empty())
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
                                               int nlag,
                                               double dlag,
                                               double toldis,
                                               double angref,
                                               const ASpaceSharedPtr& space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VectorDouble angles = VectorDouble(1);
  VectorDouble codir  = VectorDouble(ndim, 0.);
  std::vector<DirParam> dirs;
  for (int idir = 0; idir < ndir; idir++)
  {
    angles[0] = 180. * (double) idir / (double) ndir + angref;
    (void) GH::rotationGetDirection2D(angles,codir);
    double tolang = 90. / (double) ndir;
    DirParam dirparam = DirParam(nlag, dlag, toldis, tolang, 0, 0, TEST, TEST, 0.,
                                 VectorDouble(), codir, TEST, space);
    dirs.push_back(dirparam);
  }
  return dirs;
}

std::vector<DirParam> DirParam::createSeveral2D(const VectorDouble &angles,
                                                int nlag,
                                                double dlag,
                                                double toldis,
                                                double tolang,
                                                const ASpaceSharedPtr& space)
{
  std::vector<DirParam> dirs;
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();
  if (ndim != 2)
  {
    messerr("This method is limited to 2D sapce");
    return dirs;
  }

  VectorDouble anglesloc = VectorDouble(1);
  VectorDouble codir  = VectorDouble(ndim);
  int ndir = (int) angles.size();
  if (FFFF(tolang)) tolang = 90. / ndir;
  for (int idir = 0; idir < ndir; idir++)
  {
    anglesloc[0] = angles[idir];
    (void) GH::rotationGetDirection2D(anglesloc,codir);
    DirParam dirparam = DirParam(nlag, dlag, toldis, tolang, 0, 0, TEST, TEST, 0.,
                                 VectorDouble(), codir, TEST, space);
    dirs.push_back(dirparam);
  }
  return dirs;

}

/**
 * Create a set of calculation directions based on the Space definition
 * - one direction per space axis
 * - the other parameters are applied to each direction, such as:
 * @param nlag Number of lags
 * @param dlag Dimension for the lag
 * @param space Pointer to the Space definition
 * @return
 *
 * @remark: the angular tolerance is set equal to 0
 */
std::vector<DirParam> DirParam::createMultipleInSpace(int nlag, double dlag, const ASpaceSharedPtr& space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VectorDouble codir = VectorDouble(ndim);
  std::vector<DirParam> dirs;
  for (int idim = 0; idim < ndim; idim++)
  {
    VH::fill(codir, 0);
    codir[idim] = 1;
    DirParam* dirparam = DirParam::create(nlag, dlag, 0.5, 0., 0, 0, TEST, TEST, 0., VectorDouble(), codir, TEST, space);
    dirs.push_back(*dirparam);
    delete dirparam;
  }
  return dirs;
}

/****************************************************************************/
/*!
 **  Return the rank of the lag
 **
 ** \return  Rank of the lag or ITEST
 **
 ** \param[in]  dist         Distance
 **
 *****************************************************************************/
int DirParam::getLagRank(double dist) const
{
  double distloc = ABS(dist);

  /* Determine the rank of the lag */

  int ilag = -1;
  if (getFlagRegular())
  {
    ilag = (int) floor(distloc / getDPas() + 0.5);
    if (ABS(distloc - ilag * getDPas()) > getTolDist() * getDPas()) return (ITEST);
  }
  else
  {
    ilag = -1;
    for (int k = 0; k < getNLag() && ilag < 0; k++)
      if (distloc > getBreaks()[k] && distloc <= getBreaks()[k + 1]) ilag = k;
  }
  if (ilag < 0 || ilag >= getNLag()) return (ITEST);

  return ilag;
}
