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
#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

DirParam::DirParam(int ndim,
                   int npas,
                   double dpas,
                   double toldis,
                   double tolang,
                   int opt_code,
                   int idate,
                   double bench,
                   double cylrad,
                   double tolcode,
                   VectorDouble breaks,
                   VectorDouble codir,
                   VectorInt grincr)
    : _ndim(ndim),
      _nPas(npas),
      _optionCode(opt_code),
      _idate(idate),
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

DirParam::DirParam(int ndim, int npas, const VectorInt& grincr)
    : _ndim(ndim),
      _nPas(npas),
      _optionCode(0),
      _idate(0),
      _dPas(TEST),
      _bench(0.),
      _cylRad(0.),
      _tolDist(0.),
      _tolAngle(0.),
      _tolCode(0.),
      _breaks(),
      _codir(),
      _grincr(grincr)
{
  _completeDefinition();
}

DirParam::DirParam(const DirParam& r)
    : _ndim(r._ndim),
      _nPas(r._nPas),
      _optionCode(r._optionCode),
      _idate(r._idate),
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
  _completeDefinition();
}

DirParam& DirParam::operator=(const DirParam& r)
{
  if (this != &r)
  {
    _ndim = r._ndim;
    _nPas = r._nPas;
    _optionCode = r._optionCode;
    _idate = r._idate;
    _dPas = r._dPas;
    _bench = r._bench;
    _cylRad = r._cylRad;
    _tolDist = r._tolDist;
    _tolAngle = r._tolAngle;
    _tolCode = r._tolCode;
    _breaks = r._breaks;
    _codir = r._codir;
    _grincr = r._grincr;

    _completeDefinition();
  }
  return *this;
}

DirParam::~DirParam()
{
}

void DirParam::init(int ndim,
                    int npas,
                    double dpas,
                    double toldis,
                    double tolang,
                    int opt_code,
                    int idate,
                    double bench,
                    double cylrad,
                    double tolcode,
                    VectorDouble breaks,
                    VectorDouble codir,
                    VectorInt grincr)
{
  _ndim = ndim;
  _nPas = npas;
  _optionCode = opt_code;
  _idate = idate;
  _dPas = dpas;
  _tolDist = toldis;
  _tolAngle = tolang;
  _bench = bench;
  _cylRad = cylrad;
  _tolCode = tolcode;
  _breaks = breaks;
  _codir  = codir;
  _grincr = grincr;

  _completeDefinition();
}

void DirParam::_completeDefinition()
{
  if (! _breaks.empty())
  {
    if (_breaks.size() < 2) _breaks.clear();
  }
  if (_codir.empty())
  {
    _codir.resize(_ndim,0.);
    _codir[0] = 1.;
  }
  if (_grincr.empty())
  {
    _grincr.resize(_ndim,0);
    _grincr[0] = 1;
  }
}

bool DirParam::isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= getDimensionNumber())
  {
    mesArg("Space Dimension",idim,getDimensionNumber());
    return false;
  }
  return true;
}

bool DirParam::isLagValid(int ilag) const
{
  if (ilag < 0 || ilag >= getLagNumber())
  {
    mesArg("Lag Index",ilag,getLagNumber());
    return false;
  }
  return true;
}

/**
 * Set the value of the lag as computed from the Db (Grid organized)
 * @param db Db structure
 */
void DirParam::setDPas(const Db* db)
{
  if (! db->isGrid()) return;
  if (_grincr.empty()) return;
  double dpas = 0;
  for (int idim = 0; idim < _ndim; idim++)
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
    maxdist = getBreaks(getLagNumber());
  return (maxdist);
}

String DirParam::toString(int level) const
{
  std::stringstream sstr;

  sstr << "Number of lags              = " << getLagNumber() << std::endl;
  int ndim = getDimensionNumber();

  if (_grincr.empty())
  {

    // Case of a Direction defined for a non-grid Db

    sstr << toVector("Direction coefficients      = ", _codir);
    if (ndim > 1)
    {
      VectorDouble angles(ndim);
      (void) ut_angles_from_codir(ndim,1,_codir,angles);
      sstr << toVector("Direction angles (degrees)  = ", angles);
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
      sstr << "Calculation lag             = " << toDouble(getDPas()) << std::endl;
      sstr << "Tolerance on distance       = " << toDouble(100. * getTolDist())
             << " (Percent of the lag value)" << std::endl;
    }
    else
    {
      sstr << "Calculation intervals       = " << std::endl;
      for (int i = 0; i < getBreakNumber(); i++)
      {
        sstr << " - Interval " << i + 1 << " = ["
            << toInterval(getBreaks(i), getBreaks(i + 1)) << "]" << std::endl;
      }
    }
  }
  else
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

std::vector<DirParam> generateMultipleDirs(int ndim,
                                           int ndir,
                                           int npas,
                                           double dpas,
                                           double toldis)
{
  VectorDouble angles = VectorDouble(1);
  VectorDouble codir  = VectorDouble(ndim);
  std::vector<DirParam> dirs;
  for (int idir = 0; idir < ndir; idir++)
  {
    angles[0] = 180. * (double) idir / (double) ndir;
    (void) ut_angles_to_codir(ndim, 1, angles,codir);
    DirParam dirparam = DirParam(ndim, npas, dpas, toldis);
    dirparam.setTolAngle(90. / (double) ndir);
    dirparam.setCodir(codir);
    dirs.push_back(dirparam);
  }
  return dirs;
}

std::vector<DirParam> generateMultipleGridDirs(int ndim, int npas)
{
  VectorInt grincr = VectorInt(ndim);
  std::vector<DirParam> dirs;
  int ndir = ndim;
  for (int idir = 0; idir < ndir; idir++)
  {
    ut_ivector_fill(grincr, 0);
    grincr[idir] = 1;
    DirParam dirparam = DirParam(ndim, npas, grincr);
    dirs.push_back(dirparam);
  }
  return dirs;
}

