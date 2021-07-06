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
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "geoslib_f.h"

Vario::Vario(const String& calculName,
             bool flagSample,
             double scale,
             VectorDouble dates)
  : _calculName(calculName)
  , _nDim(0)
  , _nVar(0)
  , _flagSample(flagSample)
  , _scale(scale)
  , _means()
  , _vars()
  , _dates(dates)
  , _dirs()
{
}

Vario::Vario(const Vario& r)
    : _calculName(r._calculName),
      _nDim(r._nDim),
      _nVar(r._nVar),
      _flagSample(r._flagSample),
      _scale(r._scale),
      _means(r._means),
      _vars(r._vars),
      _dates(r._dates),
      _dirs(r._dirs)
{
}

Vario& Vario::operator=(const Vario& r)
{
  if (this != &r)
  {
    _calculName   = r._calculName;
    _nDim = r._nDim;
    _nVar = r._nVar;
    _flagSample = r._flagSample;
    _scale = r._scale;
    _means = r._means;
    _vars = r._vars;
    _dates = r._dates;
    _dirs = r._dirs;
  }
  return *this;
}

Vario::~Vario()
{
}

void Vario::resize(int ndim, int nvar)
{
  _nDim = ndim;
  _nVar = nvar;

  // for backwards compatibility, these arrays are updated only if their dimension
  // is not consistent with the current dimension
  if ((int) _means.size() != _nVar) _means.resize(_nVar);
  if ((int) _vars.size() != _nVar * _nVar) _vars.resize(_nVar * _nVar);

  for (int idir = 0; idir < getDirectionNumber(); idir++)
    _dirs[idir].resize(nvar, getFlagAsym());
}

double Vario::getHmax(int ivar, int jvar) const
{
  double hmax = 0.;
  for (int idir = 0; idir < (int) getDirectionNumber(); idir++)
  {
    double hloc = _dirs[idir].getHmax(ivar,jvar);
    if (hloc > hmax) hmax = hloc;
  }
  return hmax;
}

void Vario::addDirs(const Dir& dir)
{
  _dirs.push_back(dir);
}

void Vario::addDirs(const std::vector<Dir> dirs)
{
  for (int i = 0; i < (int) dirs.size(); i++)
    _dirs.push_back(dirs[i]);
}

void Vario::delDir(int rank)
{
  if (rank < 0 || rank >= getDirectionNumber()) return;
  _dirs.erase(_dirs.begin() + rank);
}

void Vario::delAllDirs()
{
  _dirs.clear();
}

int Vario::compute(Db *db,
                   const VectorDouble& means,
                   const VectorDouble& vars,
                   bool flag_grid,
                   bool flag_gen,
                   bool flag_sample,
                   bool verr_mode,
                   bool flag_model,
                   Model *model,
                   bool verbose)
{
  int ndim = db->getNDim();
  int nvar = db->getVariableNumber();
  resize(ndim, nvar);

  setMeans(means);
  setVars(vars);

  int error = variogram_compute(db, this, means, vars, flag_grid, flag_gen,
                                flag_sample, verr_mode, flag_model, model,
                                verbose);
  return error;
}

bool Vario::isCalculated() const
{
  if (getDimensionNumber() <= 0) return false;
  if (getVariableNumber() <= 0) return false;
  if (getDirectionNumber() <= 0) return false;

  for (int idir = 0; idir < (int) getDirectionNumber(); idir++)
  {
    if (! _dirs[idir].isCalculated()) return false;
  }
  return true;
}

String Vario::toString(int level) const
{
  std::stringstream sstr;
  int nvar = getVariableNumber();
  int ndir = getDirectionNumber();
  int ndim = getDimensionNumber();

  /* General parameters */

  switch (getCalculType())
  {
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
  sstr << "Number of variable(s)       = " << nvar << std::endl;
  sstr << "Number of direction(s)      = " << ndir << std::endl;
  sstr << "Space dimension             = " << ndim << std::endl;

  if (hasDate())
  {
    sstr << "Number of Date Intervals    = " << getDateNumber() << std::endl;
    sstr << toMatrix("Matrix of Bounds for Data Intervals",VectorString(),VectorString(),
                 false,2,getDateNumber(),getDates());
  }

  // Print the variance matrix

  sstr << toMatrix("Variance-Covariance Matrix",VectorString(),VectorString(),
                    0,nvar,nvar,getVars());

  /* Loop on the directions */

  sstr << std::endl;
  for (int idir=0; idir<ndir; idir++)
  {
    sstr << toTitle(1,"Direction #%d",idir+1);
    sstr << _dirs[idir].toString(level);
  }

  return sstr.str();
}
/**
 * Convert the Calculation Name into a Calculation Type (enum)
 * @param calcul_name Input calculation name to be identified
 * @return
 */
int identifyCalculType(const String& calcul_name)
{
  int calcul_type;

  if (!strcmp(calcul_name.c_str(), "vg"))
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

int identifyFlagAsym(const String& calcul_name)
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

double Vario::getVars(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _vars[_getAddress(ivar, jvar)];
}

double Vario::getVars(int i) const
{
  if (! _isBivariableValid(i)) return TEST;
  return _vars[i];
}

double Vario::getDates(int idate, int icas) const
{
  if (!_isDateValid(idate)) return 0.;
  return _dates[2 * idate + icas];
}

void Vario::setVars(int i, double value)
{
  if (! _isBivariableValid(i)) return;
  _vars[i] = value;
}

void Vario::setVars(int ivar, int jvar, double value)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  _vars[_getAddress(ivar,jvar)] = value;
}

int Vario::getLagNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirs[idir].getLagNumber();
}

VectorDouble Vario::getGg(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirs[idir].getGg(ivar,jvar);
}

VectorDouble Vario::getHh(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirs[idir].getHh(ivar,jvar);
}

VectorDouble Vario::getSw(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirs[idir].getSw(ivar,jvar);
}

VectorDouble Vario::getCodir(int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirs[idir].getCodir();
}

int Vario::_getAddress(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return ITEST;
  return ivar + _nVar * jvar;
}

bool Vario::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= _nVar)
  {
    mesArg("Variable Index",ivar,_nVar);
    return false;
  }
  return true;
}

bool Vario::_isBivariableValid(int i) const
{
  if (i < 0 || i >= _nVar * _nVar)
  {
    mesArg("Multivariate Index",i,_nVar * _nVar);
    return false;
  }
  return true;
}

bool Vario::_isDirectionValid(int idir) const
{
  if (idir < 0 || idir >= getDirectionNumber())
  {
    mesArg("Direction Index",idir,getDirectionNumber());
    return false;
  }
  return true;
}

bool Vario::_isDateValid(int idate) const
{
  if (! hasDate()) return false;
  if (idate < 0 || idate >= getDateNumber())
  {
    mesArg("Date Index",idate,getDateNumber());
    return false;
  }
  return true;
}
