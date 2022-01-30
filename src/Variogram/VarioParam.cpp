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
#include "Variogram/VarioParam.hpp"
#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

VarioParam::VarioParam(double scale,
                       VectorDouble dates)
  : AStringable()
  , IClonable()
  , _scale(scale)
  , _dates(dates)
  , _dirparams()
{
}

VarioParam::VarioParam(const VarioParam& VarioParam, const VectorInt& dircols)
    : AStringable(),
      IClonable(),
      _scale(),
      _dates(),
      _dirparams()
{
    _scale = VarioParam.getScale();
    _dates = VarioParam.getDates();

    for (int idir = 0; idir < (int) dircols.size(); idir++)
    {
      _dirparams.push_back(VarioParam.getDirParam(dircols[idir]));
    }
}

VarioParam::VarioParam(const VarioParam& r)
    : AStringable(r),
      _scale(r._scale),
      _dates(r._dates),
      _dirparams(r._dirparams)
{
}

VarioParam& VarioParam::operator=(const VarioParam& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _scale = r._scale;
    _dates = r._dates;
    _dirparams  = r._dirparams;
  }
  return *this;
}

VarioParam::~VarioParam()
{
}

void VarioParam::addDirs(const DirParam& dirparam)
{
  _dirparams.push_back(dirparam);
}

void VarioParam::addMultiDirs(const std::vector<DirParam>& dirparams)
{
  for (int i = 0; i < (int) dirparams.size(); i++)
    _dirparams.push_back(dirparams[i]);
}

void VarioParam::delDir(int rank)
{
  if (rank < 0 || rank >= getDirectionNumber()) return;
  _dirparams.erase(_dirparams.begin() + rank);
}

void VarioParam::delAllDirs()
{
  _dirparams.clear();
}

String VarioParam::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  // Print the Main part

  sstr << toStringMain(strfmt);

  /* Loop on the directions */

  for (int idir=0; idir<getDirectionNumber(); idir++)
  {
    sstr << toTitle(1,"Direction #%d",idir+1);
    sstr << _dirparams[idir].toString(strfmt);
  }

  return sstr.str();
}

String VarioParam::toStringMain(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ndir = getDirectionNumber();

  /* General parameters */

  sstr << "Number of direction(s)      = " << ndir << std::endl;
  sstr << "Space dimension             = " << getDimensionNumber() << std::endl;

  if (hasDate())
  {
    sstr << "Number of Date Intervals    = " << getDateNumber() << std::endl;
    sstr << toMatrix("Matrix of Bounds for Data Intervals",VectorString(),VectorString(),
                 false,2,getDateNumber(),getDates());
  }
  return sstr.str();
}

double VarioParam::getDate(int idate, int icas) const
{
  if (!_isDateValid(idate)) return 0.;
  return _dates[2 * idate + icas];
}

int VarioParam::getLagNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirparams[idir].getLagNumber();
}

VectorDouble VarioParam::getCodir(int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirparams[idir].getCodir();
}

bool VarioParam::_isDirectionValid(int idir) const
{
  if (idir < 0 || idir >= getDirectionNumber())
  {
    mesArg("Direction Index",idir,getDirectionNumber());
    return false;
  }
  return true;
}

bool VarioParam::_isDateValid(int idate) const
{
  if (! hasDate()) return false;
  if (idate < 0 || idate >= getDateNumber())
  {
    mesArg("Date Index",idate,getDateNumber());
    return false;
  }
  return true;
}

VectorDouble VarioParam::_getDirectionInterval(int idir) const
{
  VectorDouble bounds(2);
  if (idir < 0 || idir >= getDimensionNumber())
  {
    bounds[0] = 0;
    bounds[1] = getDirectionNumber();
  }
  else
  {
    bounds[0] = idir;
    bounds[1] = idir + 1;
  }
  return bounds;
}

void VarioParam::setDPas(int idir,const DbGrid* db)
{
  if (! _isDirectionValid(idir)) return;
  _dirparams[idir].setDPas(db);
}

void VarioParam::setGrincr(int idir, const VectorInt& grincr)
{
  if (! _isDirectionValid(idir)) return;
  _dirparams[idir].setGrincr(grincr);

}
