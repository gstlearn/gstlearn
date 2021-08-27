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

VarioParam::VarioParam(double scale, bool flagSample, VectorDouble dates)
  : AStringable()
  , IClonable()
  , _calculName("undefined")
  , _nDim(0)
  , _flagSample(flagSample)
  , _scale(scale)
  , _dates(dates)
  , _dirparams()
{
}

VarioParam::VarioParam(const VarioParam& VarioParam, const VectorInt& dircols)
    : AStringable(),
      IClonable(),
      _calculName("undefined"),
      _nDim(0),
      _flagSample(),
      _scale(),
      _dates(),
      _dirparams()
{
    _calculName = VarioParam.getCalculName();
    _nDim = VarioParam.getDimensionNumber();
    _flagSample = VarioParam.getFlagSample();
    _scale = VarioParam.getScale();
    _dates = VarioParam.getDates();

    _dirparams = std::vector<DirParam>();
    for (int idir = 0; idir < (int) dircols.size(); idir++)
    {
      _dirparams.push_back(VarioParam.getDirParam(dircols[idir]));
    }
}

VarioParam::VarioParam(const VarioParam& r)
    : _calculName(r._calculName),
      _nDim(r._nDim),
      _flagSample(r._flagSample),
      _scale(r._scale),
      _dates(r._dates),
      _dirparams(r._dirparams)
{
}

VarioParam& VarioParam::operator=(const VarioParam& r)
{
  if (this != &r)
  {
    _calculName   = r._calculName;
    _nDim = r._nDim;
    _flagSample = r._flagSample;
    _scale = r._scale;
    _dates = r._dates;
    _dirparams  = r._dirparams;
  }
  return *this;
}

VarioParam::~VarioParam()
{
}

IClonable* VarioParam::clone() const
{
  return new VarioParam(*this);
}

void VarioParam::addDirs(const DirParam& dirparam)
{
  _dirparams.push_back(dirparam);
}

void VarioParam::addDirs(const std::vector<DirParam> dirparams)
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

String VarioParam::toString(int level) const
{
  std::stringstream sstr;
  int ndir = getDirectionNumber();
  int ndim = getDimensionNumber();

  /* General parameters */

  switch (getCalculType())
  {
    case CALCUL_UNDEFINED:
      sstr << toTitle(0,"Undefined");
      break;

    case CALCUL_VARIOGRAM:
      sstr << toTitle(0,"Variogram characteristics");
      break;

    case CALCUL_MADOGRAM:
      sstr << toTitle(0,"Madogram characteristics");
      break;

    case CALCUL_RODOGRAM:
      sstr << toTitle(0,"Rodogram characteristics");
      break;

    case CALCUL_POISSON:
      sstr << toTitle(0,"Poisson variogram characteristics");
      break;

    case CALCUL_COVARIANCE:
      sstr << toTitle(0,"Covariance characteristics");
      break;

    case CALCUL_COVARIANCE_NC:
      sstr << toTitle(0,"Non-centered Covariance characteristics");
      break;

    case CALCUL_COVARIOGRAM:
      sstr << toTitle(0,"Transitive Covariogram characteristics");
      break;

    case CALCUL_GENERAL1:
      sstr << toTitle(0,"Generalized Variogram of order 1 characteristics");
      break;

    case CALCUL_GENERAL2:
      sstr << toTitle(0,"Generalized Variogram of order 2 characteristics");
      break;

    case CALCUL_GENERAL3:
      sstr << toTitle(0,"Generalized Variogram of order 3 characteristics");
      break;

    case CALCUL_ORDER4:
      sstr << toTitle(0,"Order-4 Variogram");
      break;

    case CALCUL_TRANS1:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/G1");
      break;

    case CALCUL_TRANS2:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/G2");
      break;

    case CALCUL_BINORMAL:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/sqrt(G1*G2)");
      break;
  }
  sstr << "Number of direction(s)      = " << ndir << std::endl;
  sstr << "Space dimension             = " << ndim << std::endl;

  if (hasDate())
  {
    sstr << "Number of Date Intervals    = " << getDateNumber() << std::endl;
    sstr << toMatrix("Matrix of Bounds for Data Intervals",VectorString(),VectorString(),
                 false,2,getDateNumber(),getDates());
  }

  if (getCalculType() == CALCUL_UNDEFINED) return sstr.str();

  /* Loop on the directions */

  sstr << std::endl;
  for (int idir=0; idir<ndir; idir++)
  {
    sstr << toTitle(1,"Direction #%d",idir+1);
    sstr << _dirparams[idir].toString(level);
  }

  return sstr.str();
}
/**
 * Convert the Calculation Name into a Calculation Type (enum)
 * @param calcul_name Input calculation name to be identified
 * @return
 */
int identifyCalculTypeN(const String& calcul_name)
{
  int calcul_type;

  if (!strcmp(calcul_name.c_str(), "undefined"))
    calcul_type = CALCUL_UNDEFINED;
  else if (!strcmp(calcul_name.c_str(), "vg"))
    calcul_type = CALCUL_VARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "cov"))
    calcul_type = CALCUL_COVARIANCE;
  else if (!strcmp(calcul_name.c_str(), "covnc"))
    calcul_type = CALCUL_COVARIANCE_NC;
  else if (!strcmp(calcul_name.c_str(), "covg"))
    calcul_type = CALCUL_COVARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "mado"))
    calcul_type = CALCUL_MADOGRAM;
  else if (!strcmp(calcul_name.c_str(), "rodo"))
    calcul_type = CALCUL_RODOGRAM;
  else if (!strcmp(calcul_name.c_str(), "poisson"))
    calcul_type = CALCUL_POISSON;
  else if (!strcmp(calcul_name.c_str(), "general1"))
    calcul_type = CALCUL_GENERAL1;
  else if (!strcmp(calcul_name.c_str(), "general2"))
    calcul_type = CALCUL_GENERAL2;
  else if (!strcmp(calcul_name.c_str(), "general3"))
    calcul_type = CALCUL_GENERAL3;
  else if (!strcmp(calcul_name.c_str(), "order4"))
    calcul_type = CALCUL_ORDER4;
  else if (!strcmp(calcul_name.c_str(), "trans1"))
    calcul_type = CALCUL_TRANS1;
  else if (!strcmp(calcul_name.c_str(), "trans2"))
    calcul_type = CALCUL_TRANS2;
  else if (!strcmp(calcul_name.c_str(), "binormal"))
    calcul_type = CALCUL_BINORMAL;
  else
  {
    messerr("Invalid variogram calculation name : %s", calcul_name.c_str());
    messerr("The only valid names are:");
    messerr("vg       : Variogram");
    messerr("cov      : Covariance");
    messerr("covnc    : Non-centered ergodic covariance");
    messerr("covg     : Covariogram");
    messerr("mado     : Madogram");
    messerr("rodo     : Rodogram");
    messerr("poisson  : Poisson");
    messerr("general1 : Generalized variogram of order 1");
    messerr("general2 : Generalized variogram of order 2");
    messerr("general3 : Generalized variogram of order 3");
    messerr("order4   : Variogram of order 4");
    messerr("trans1   : Cross-to-Simple Variogram G12/G1");
    messerr("trans2   : Cross-to-Simple Variogram G12/G1");
    messerr("binormal : Cross-to-Simple Variogram G12/sqrt(G1*G2)");

    calcul_type = CALCUL_UNDEFINED;
  }
  return calcul_type;
}

int identifyFlagAsymN(const String& calcul_name)
{
  int flagAsym = 0;

  switch (identifyCalculType(calcul_name))
  {
    case CALCUL_VARIOGRAM:
    case CALCUL_MADOGRAM:
    case CALCUL_RODOGRAM:
    case CALCUL_POISSON:
    case CALCUL_GENERAL1:
    case CALCUL_GENERAL2:
    case CALCUL_GENERAL3:
    case CALCUL_ORDER4:
    case CALCUL_TRANS1:
    case CALCUL_TRANS2:
    case CALCUL_BINORMAL:
      flagAsym = 0;
      break;

    case CALCUL_COVARIANCE:
    case CALCUL_COVARIANCE_NC:
    case CALCUL_COVARIOGRAM:
      flagAsym = 1;
      break;
  }
  return flagAsym;
}

double VarioParam::getDates(int idate, int icas) const
{
  if (!_isDateValid(idate)) return 0.;
  return _dates[2 * idate + icas];
}

int VarioParam::getLagNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirparams[idir].getLagNumber();
}

int VarioParam::getLagTotalNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirparams[idir].getLagTotalNumber();
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
