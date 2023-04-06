/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f_private.h"

#include "Variogram/VarioParam.hpp"
#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Stats/Classical.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"

VarioParam::VarioParam(double scale,
                       const VectorDouble& dates,
                       const Faults* faults)
  : AStringable()
  , ICloneable()
  , _scale(scale)
  , _dates(dates)
  , _dirparams()
  , _faults(faults)
{
}

VarioParam::VarioParam(const VarioParam &VarioParam,
                       const VectorInt &dircols,
                       const Faults *faults)
    : AStringable(),
      ICloneable(),
      _scale(),
      _dates(),
      _dirparams(),
      _faults(faults)
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
      _dirparams(r._dirparams),
      _faults(r._faults)
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
    _faults = r._faults;
  }
  return *this;
}

VarioParam::~VarioParam()
{
}

bool VarioParam::isDefinedForGrid() const
{
  int ndir = getDirectionNumber();
  if (ndir <= 0) return false;
  return _dirparams[0].isDefinedForGrid();
}

bool VarioParam::_validDefinedFromGrid(const DirParam& dirparam) const
{
  int ndir = getDirectionNumber();
  bool definedFromGrid = dirparam.isDefinedForGrid();
  if (ndir > 0)
  {
    for (int idir = 0; idir < ndir; idir++)
    {
      if (_dirparams[idir].isDefinedForGrid() != definedFromGrid)
      {
        messerr("The current 'dirParam' cannot be added to 'varioParam'");
        if (_dirparams[idir].isDefinedForGrid())
          messerr("Element (%d) is defined using Grid definition",idir+1);
        else
          messerr("Element(%d) is defined NOT using Grid definition",idir+1);

        if (definedFromGrid)
          messerr("Current 'dirparam' is defined using Grid definition");
        else
          messerr("Current 'dirparam' is defined NOT using Grid definition");
        return false;
      }
    }
  }
  return true;
}

VarioParam* VarioParam::createOmniDirection(int npas,
                                            double dpas,
                                            double toldis,
                                            int opt_code,
                                            int idate,
                                            double bench,
                                            double cylrad,
                                            double tolcode,
                                            const VectorDouble& breaks,
                                            double scale,
                                            const VectorDouble& dates,
                                            const ASpace* space)
{
  DirParam* dir = DirParam::createOmniDirection(npas, dpas, toldis,
                                                opt_code, idate, bench, cylrad,
                                                tolcode, breaks, space);
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addDir(*dir);
  return varioparam;
}

VarioParam* VarioParam::createMultiple(int ndir,
                                       int npas,
                                       double dpas,
                                       double toldis,
                                       double scale,
                                       const VectorDouble &dates,
                                       const ASpace* space)
{
  std::vector<DirParam> dirs = DirParam::createMultiple(ndir, npas, dpas,
                                                        toldis, space);
  if (dirs.empty()) return nullptr;
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addMultiDirs(dirs);
  return varioparam;
}

VarioParam* VarioParam::createMultipleFromGrid(int npas,
                                               double scale,
                                               const VectorDouble &dates,
                                               const ASpace* space)
{
  std::vector<DirParam> dirs = DirParam::createMultipleFromGrid(npas, space);
  if (dirs.empty()) return nullptr;
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addMultiDirs(dirs);
  return varioparam;
}

VarioParam* VarioParam::createFromSpaceDimension(int npas,
                                                 double dpas,
                                                 double toldis,
                                                 double tolang,
                                                 double scale,
                                                 const VectorDouble &dates,
                                                 const ASpace *space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VarioParam* varioparam = new VarioParam(scale, dates);

  for (int idim = 0; idim < ndim; idim++)
  {
    DirParam dirparam(npas, dpas, toldis, tolang);
    VectorDouble codir(ndim,0.);
    codir[idim] = 1.;
    dirparam.setCodir(codir);

    varioparam->addDir(dirparam);
  }
  return varioparam;
}

void VarioParam::addDir(const DirParam& dirparam)
{
  if (! _validDefinedFromGrid(dirparam)) return;
  _dirparams.push_back(dirparam);
}

void VarioParam::addMultiDirs(const std::vector<DirParam>& dirparams)
{
  for (int i = 0; i < (int) dirparams.size(); i++)
  {
    if (! _validDefinedFromGrid(dirparams[i])) return;
    _dirparams.push_back(dirparams[i]);
  }
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
  if (getDirectionNumber() <= 0) return sstr.str();

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

VectorDouble VarioParam::getCodirs(int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirparams[idir].getCodirs();
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

int VarioParam::getDimensionNumber() const
{
  if (getDirectionNumber() <= 0) return 0;
  return _dirparams[0].getNDim();
}
