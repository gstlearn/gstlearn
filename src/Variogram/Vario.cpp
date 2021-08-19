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
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

Vario::Vario(double scale,
             bool flagSample,
             VectorDouble dates)
  : AStringable()
  , ASerializable()
  , IClonable()
  , _calculName("undefined")
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

/**
 * Constructing a Variogram by extracting a subset of variables and/or directions
 * from an Input variogram
 *
 * @param vario The input variogram
 * @param varcols Vector of variable ranks (starting from 0)
 * @param dircols Vector of direction ranks (starting from 0)
 * @param flagVario True if the output should be a variogram (rather than a covariance)
 */
Vario::Vario(const Vario& vario,
             const VectorInt& varcols,
             const VectorInt& dircols,
             bool flagVario)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _calculName(),
      _nDim(0),
      _nVar(0),
      _flagSample(false),
      _scale(0.),
      _means(),
      _vars(),
      _dates(),
      _dirs()
{
  if (flagVario)
    _calculName = "vg";
  else
    _calculName = vario.getCalculName();
  _nDim = vario.getDimensionNumber();
  _flagSample = vario.getFlagSample();
  _scale = vario.getScale();
  _dates = vario.getDates();

  VectorInt selvars;
  if (varcols.empty())
    selvars = ut_ivector_sequence(vario.getVariableNumber());
  else
  {
    selvars = varcols;
    for (int i = 0; i < (int) varcols.size(); i++)
    {
      if (selvars[i] < 0 || selvars[i] >= vario.getVariableNumber())
      {
        messerr(
            "Argument 'varcols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, selvars[i], vario.getVariableNumber());
        my_throw("Error when extracting a variogram");
      }
    }
  }

  VectorInt seldirs;
  if (seldirs.empty())
    seldirs = ut_ivector_sequence(vario.getDirectionNumber());
  else
  {
    seldirs = dircols;
    for (int i = 0; i < (int) vario.getDimensionNumber(); i++)
    {
      if (seldirs[i] < 0 || seldirs[i] >= vario.getDirectionNumber())
      {
        messerr(
            "Argument 'dircols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, seldirs[i], vario.getDirectionNumber());
        my_throw("Error when extracting a variogram");
      }
    }
  }

  // Dimension the new 'vario'
  _nDim = vario.getDimensionNumber();
  _nVar = (int) selvars.size();

  // Reset Mean and variance arrays (only if variable number has been modified)
  if (_nVar != vario.getDimensionNumber())
  {
    _means.resize(_nVar);
    for (int ivar = 0; ivar < _nVar; ivar++)
      setMeans(ivar, vario.getMeans(selvars[ivar]));

    _vars.resize(_nVar * _nVar);
    for (int ivar = 0; ivar < _nVar; ivar++)
      for (int jvar = 0; jvar < _nVar; jvar++)
        setVars(ivar, jvar, vario.getVars(selvars[ivar], selvars[jvar]));
  }
  else
  {
    _means = vario.getMeans();
    _vars = vario.getVars();
  }

  // Add the directions

  for (int idir = 0; idir < (int) seldirs.size(); idir++)
  {
    const Dir& dirFrom = vario.getDirs(seldirs[idir]);

    // Add the direction to the new Variogram
    addDirs(dirFrom);

    // Resize it to the correct number of variables
    _dirs[idir].internalResize(_nVar, getFlagAsym());

    // Load the relevant information
    int npas = _dirs[idir].getNPas();
    for (int ivar = 0; ivar < _nVar; ivar++)
    {
      int ivarp = selvars[ivar];
      for (int jvar = 0; jvar < _nVar; jvar++)
      {
        int jvarp = selvars[jvar];
        for (int ipas = 0; ipas < npas; ipas++)
        {
          int iadFrom = dirFrom.getAddress(ivarp, jvarp, ipas, false, 1);
          double sw = dirFrom.getSw(iadFrom);
          double gg = dirFrom.getGg(iadFrom);
          double hh = dirFrom.getHh(iadFrom);

          if (flagVario && dirFrom.getFlagAsym())
          {
            int iadFrom = dirFrom.getAddress(ivarp, jvarp, ipas, false, -1);
            sw = (sw + dirFrom.getSw(iadFrom)) / 2.;
            gg = (gg + dirFrom.getGg(iadFrom)) / 2.;
            hh = (ABS(hh) + ABS(dirFrom.getHh(iadFrom))) / 2.;
          }

          setSw(idir, ivar, jvar, ipas, sw);
          setGg(idir, ivar, jvar, ipas, getVars(ivar, jvar) - gg);
          setHh(idir, ivar, jvar, ipas, hh);
        }
      }
    }
  }
}

Vario::Vario(const String& neutralFileName,
             bool verbose)
  : AStringable()
  , ASerializable()
  , IClonable()
  , _calculName()
  , _nDim(0)
  , _nVar(0)
  , _flagSample(false)
  , _scale(0.)
  , _means()
  , _vars()
  , _dates()
  , _dirs()
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
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

void Vario::internalResize(int ndim, int nvar, const String& calculName)
{
  if (ndim <= 0 || nvar <= 0 || identifyCalculType(calculName) == CALCUL_UNDEFINED)
  {
    messerr("The Internal Dimension assigned to the variogram cannot be calculated:");
    messerr("- Space Dimension = %d",ndim);
    messerr("- Number of variables = %d",nvar);
    messerr("- The Calculation type has not been defined yet");
    my_throw("Error in Internal Variogram dimensioning");
  }

  _nDim = ndim;
  _nVar = nvar;
  _calculName = calculName;

  // for backwards compatibility, these arrays are updated only if their dimension
  // is not consistent with the current dimension
  if ((int) _means.size() != _nVar) _means.resize(_nVar);
  if ((int) _vars.size() != _nVar * _nVar) _initVars();

  for (int idir = 0; idir < getDirectionNumber(); idir++)
    _dirs[idir].internalResize(_nVar, getFlagAsym());
}

IClonable* Vario::clone() const
{
  return new Vario(*this);
}

double Vario::getHmax(int ivar, int jvar, int idir) const
{
  VectorDouble ivb = _getVariableInterval(ivar);
  VectorDouble jvb = _getVariableInterval(jvar);
  VectorDouble idb = _getDirectionInterval(idir);

  double hmax = 0.;
  for (int idir = idb[0]; idir < idb[1]; idir++)
    for (int ivar = ivb[0]; ivar < ivb[1]; ivar++)
      for (int jvar = jvb[0]; jvar < jvb[1]; jvar++)
      {
        double hloc = _dirs[idir].getHmax(ivar,jvar);
        if (hloc > hmax) hmax = hloc;
      }
  return hmax;
}

VectorDouble Vario::getHRange(int ivar, int jvar, int idir) const
{
  VectorDouble ivb = _getVariableInterval(ivar);
  VectorDouble jvb = _getVariableInterval(jvar);
  VectorDouble idb = _getDirectionInterval(idir);

  VectorDouble vec(2);
  vec[0] =  1.e30;
  vec[1] = -1.e30;
  for (int idir = idb[0]; idir < idb[1]; idir++)
    for (int ivar = ivb[0]; ivar < ivb[1]; ivar++)
      for (int jvar = jvb[0]; jvar < jvb[1]; jvar++)
      {
        double hloc = _dirs[idir].getHmax(ivar, jvar);
        if (hloc < vec[0]) vec[0] = hloc;
        if (hloc > vec[1]) vec[1] = hloc;
      }
  return vec;
}

double Vario::getGmax(int ivar, int jvar, int idir, bool flagAbs, bool flagSill) const
{
  VectorDouble ivb = _getVariableInterval(ivar);
  VectorDouble jvb = _getVariableInterval(jvar);
  VectorDouble idb = _getDirectionInterval(idir);

  double gmax = 0.;

  for (int idir = idb[0]; idir < idb[1]; idir++)
    for (int ivar = ivb[0]; ivar < ivb[1]; ivar++)
      for (int jvar = jvb[0]; jvar < jvb[1]; jvar++)
      {
        double gloc = _dirs[idir].getGmax(ivar, jvar);
        if (flagAbs) gloc = ABS(gloc);
        if (gmax < gloc) gmax = gloc;
        if (flagSill)
        {
          double sill = ABS(getVars(ivar, jvar));
          if (gmax < sill) gmax = sill;
        }
      }
  return gmax;
}

VectorDouble Vario::getGRange(int ivar, int jvar, int idir, bool flagSill) const
{
  VectorDouble ivb = _getVariableInterval(ivar);
  VectorDouble jvb = _getVariableInterval(jvar);
  VectorDouble idb = _getDirectionInterval(idir);

  VectorDouble vec(2);
  vec[0] =  1.e30;
  vec[1] = -1.e30;

  for (int idir = idb[0]; idir < idb[1]; idir++)
    for (int ivar = ivb[0]; ivar < ivb[1]; ivar++)
      for (int jvar = jvb[0]; jvar < jvb[1]; jvar++)
      {
        double gloc = _dirs[idir].getGmax(ivar, jvar);
        if (gloc < vec[0]) vec[0] = gloc;
        if (gloc > vec[1]) vec[1] = gloc;
        if (flagSill)
        {
          double sill = getVars(ivar, jvar);
          if (sill < vec[0]) vec[0] = sill;
          if (sill > vec[1]) vec[1] = sill;
        }
      }
  return vec;
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
                   const String& calculName,
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
  internalResize(ndim, nvar, calculName);

  setMeans(means);
  setVars(vars);

  return _variogram_compute(db, this, means, vars, flag_grid, flag_gen,
                            flag_sample, verr_mode, flag_model, model, verbose);
}

/**
 * Compute the Variogram of the Indicator of the input variable
 * It is assumed that the Single input variable is already a Facies
 * (starting from 1)
 */
int Vario::computeIndic(Db *db,
                        const String& calculName,
                        const VectorDouble& means,
                        const VectorDouble& vars,
                        bool flag_grid,
                        bool flag_gen,
                        bool flag_sample,
                        bool verr_mode,
                        bool flag_model,
                        Model *model,
                        bool verbose,
                        int nfacmax)
{
  int ndim = db->getNDim();
  int nvar = db->getVariableNumber();
  if (nvar != 1)
  {
    messerr("This method is only considered for a Single input Variable");
    return 1;
  }

  // Calculate the number of Facies in 'Db'
  VectorDouble props = dbStatisticsFacies(db);
  int nclass = static_cast<int> (props.size());
  if (nclass <= 0 || nclass > nfacmax)
  {
    messerr("The input variable should exhibit Facies");
    messerr("Number of Facies (%d) should be positive and smaller than 'nfacmax'",
            nclass);
    messerr("Note: the value of 'nfacmax'(%d) can be changed in argument list",
            nfacmax);
    return 1;
  }

  // Translate the 'Facies' into 'categories'
  Limits limits = Limits(nclass);
  int iatt = db->getAttribute(LOC_Z, 0);
  if (limits.toIndicator(db, iatt))
  {
    messerr("Problem when translating Facies into Categories");
    return 1;
  }

  internalResize(ndim, nclass, calculName);

  // Calculate the variogram of Indicators
  if (compute(db, calculName, props, VectorDouble(), flag_grid, flag_gen,
              flag_sample, verr_mode, flag_model, model, verbose))
  {
    messerr("Error when calculating the Variogram of Indicators");
    return 1;
  }

  // Delete the Indicators (created locally)
  db->deleteFieldByLocator(LOC_Z);
  db->setLocatorByAttribute(iatt, LOC_Z);

  return 0;
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

  if (getCalculType() == CALCUL_UNDEFINED) return sstr.str();

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

double Vario::getMeans(int ivar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  return _means[ivar];
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

void Vario::_initMeans()
{
  _means.resize(_nVar);
  for (int ivar = 0; ivar < _nVar; ivar++)
    _means[ivar] = 0.;
}

void Vario::setMeans(const VectorDouble& means)
{
  if (_means.empty()) _initMeans();
  if (! means.empty() && (int) means.size() == _nVar)
    _means = means;
}

void Vario::setMeans(int ivar, double mean)
{
  if (_means.empty()) _initMeans();
  if (! _isVariableValid(ivar)) return;
  _means[ivar] = mean;
}

void Vario::_initVars()
{
  _vars.resize(_nVar * _nVar);
  int ecr = 0;
  for (int ivar = 0; ivar < _nVar; ivar++)
    for (int jvar = 0; jvar < _nVar; jvar++)
      _vars[ecr++] = (ivar == jvar);
}

void Vario::setVars(const VectorDouble& vars)
{
  if (_vars.empty()) _initVars();
  if (! vars.empty() && (int) vars.size() == _nVar * _nVar)
    _vars = vars;
}

void Vario::setVars(int i, double value)
{
  if (_vars.empty()) _initVars();
  if (! _isBivariableValid(i)) return;
  _vars[i] = value;
}

void Vario::setVars(int ivar, int jvar, double value)
{
  if (_vars.empty()) _initVars();
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  _vars[_getAddress(ivar,jvar)] = value;
}

int Vario::getLagNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirs[idir].getLagNumber();
}

int Vario::getLagTotalNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirs[idir].getLagTotalNumber();
}

int Vario::getSize(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirs[idir].getSize();
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
int Vario::getCenter(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return ITEST;
  return _dirs[idir].getCenter(ivar,jvar);
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

int Vario::deSerialize(const String& filename, bool verbose)
{
  int ndim, nvar, ndir, npas, opt_code, flag_calcul, flag_regular;
  double dpas, tolang, scale, tolcode, toldis;
  VectorDouble codir, grloc, vars;
  VectorInt grincr;

  // Open the Neutral File

  if (_fileOpen(filename, "Vario", "r", verbose)) return 1;

  /* Create the Vario structure */

  if (_recordRead("Space Dimension", "%d", &ndim)) goto label_end;
  if (_recordRead("Number of Variables", "%d", &nvar)) goto label_end;
  if (_recordRead("Number of Variogram Directions", "%d", &ndir)) goto label_end;
  if (_recordRead("Scale", "%lf", &scale)) goto label_end;

  /* Read the variances (optional) */

  if (_recordRead("Variogram calculation Option", "%d", &flag_calcul))  goto label_end;
  vars.resize(nvar * nvar);
  if (flag_calcul)
  {
    int ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++, ecr++)
        if (_recordRead("Experimental Variance term", "%lf", &vars[ecr])) goto label_end;
  }

  /* Initialize the variogram structure */

  internalResize(ndim, nvar, "vg");
  setScale(scale);
  setVars(vars);

  /* Reading the variogram calculation directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    if (_recordRead("Regular Variogram Calculation", "%d", &flag_regular)) goto label_end;
    if (_recordRead("Number of Variogram Lags", "%d", &npas)) goto label_end;
    if (_recordRead("Variogram Code Option", "%d", &opt_code)) goto label_end;
    if (_recordRead("Tolerance on Code", "%lf", &tolcode)) goto label_end;
    if (_recordRead("Lag Value", "%lf", &dpas)) goto label_end;
    if (_recordRead("Tolerance on Distance", "%lf", &toldis)) goto label_end;
    if (_recordRead("Tolerance on Direction", "%lf", &tolang)) goto label_end;
    codir.resize(ndim);
    grincr.resize(ndim);
    grloc.resize(ndim);
    for (int idim = 0; idim < ndim; idim++)
      if (_recordRead("Direction vector", "%lf", &codir[idim]))
        goto label_end;
    for (int idim = 0; idim < ndim; idim++)
      if (_recordRead("Grid Increment", "%lf", &grloc[idim]))
        goto label_end;

    Dir dir = Dir(ndim);
    for (int idim = 0; idim < ndim; idim++) grincr[idim] = (int) grloc[idim];
    dir.init(ndim, npas, dpas, toldis, tolang, getFlagAsym(),
             opt_code, 0, TEST, TEST, tolcode, VectorDouble(), codir, grincr);


    /* Read the arrays of results (optional) */

    if (flag_calcul)
    {
      dir.internalResize(nvar, getFlagAsym());
      for (int i = 0; i < dir.getSize(); i++)
      {
        double sw, hh, gg;
        if (_recordRead("Experimental Variogram Weight", "%lf", &sw)) goto label_end;
        dir.setSw(i, sw);
        if (_recordRead("Experimental Variogram Distance", "%lf", &hh)) goto label_end;
        dir.setHh(i, hh);
        if (_recordRead("Experimental Variogram Value", "%lf", &gg)) goto label_end;
        dir.setGg(i, gg);
      }
    }
    addDirs(dir);
  }

  label_end:

  /* Close the file */
  _fileClose(verbose);

  return 0;
}

int Vario::serialize(const String& filename, bool verbose) const
{
  double value;
  static int flag_calcul = 1;

  if (_fileOpen(filename, "Vario", "w", verbose)) return 1;

  /* Write the Vario structure */

  _recordWrite("%d", getDimensionNumber());
  _recordWrite("#", "Space Dimension");
  _recordWrite("%d", getVariableNumber());
  _recordWrite("#", "Number of variables");
  _recordWrite("%d", getDirectionNumber());
  _recordWrite("#", "Number of directions");
  _recordWrite("%lf", getScale());
  _recordWrite("#", "Scale");
  _recordWrite("%d", flag_calcul);
  _recordWrite("#", "Calculation Flag");

  /* Dumping the Variances */

  if (flag_calcul)
  {
    _recordWrite("#", "Variance");
    for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    {
      for (int jvar = 0; jvar < getVariableNumber(); jvar++)
        _recordWrite("%lf", getVars(ivar,jvar));
      _recordWrite("\n");
    }
  }

  /* Loop on the directions */

  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    _recordWrite("#", "Direction characteristics");
    const Dir& dir = getDirs(idir);
    _recordWrite("%d", dir.getLagRegular());
    _recordWrite("#", "Regular lags");
    _recordWrite("%d", dir.getNPas());
    _recordWrite("#", "Number of lags");
    _recordWrite("%d", dir.getOptionCode());
    _recordWrite("%lf", dir.getTolCode());
    _recordWrite("#", "Code selection: Option - Tolerance");
    _recordWrite("%lf", dir.getDPas());
    _recordWrite("#", "Lag value");
    _recordWrite("%lf", dir.getTolDist());
    _recordWrite("#", "Tolerance on distance");
    _recordWrite("%lf", dir.getTolAngle());
    _recordWrite("#", "Tolerance on angle");

    for (int idim = 0; idim < getDimensionNumber(); idim++)
      _recordWrite("%lf", dir.getCodir(idim));
    _recordWrite("#", "Direction coefficients");

    for (int idim = 0; idim < getDimensionNumber(); idim++)
      _recordWrite("%lf", (double) dir.getGrincr(idim));
    _recordWrite("#", "Direction increments on grid");

    if (!flag_calcul) continue;
    _recordWrite("#", "Variogram results (Weight, Distance, Variogram)");
    for (int i = 0; i < dir.getSize(); i++)
    {
      value = FFFF(dir.getSw(i)) ? 0. : dir.getSw(i);
      _recordWrite("%lf", value);
      value = FFFF(dir.getHh(i)) ? 0. : dir.getHh(i);
      _recordWrite("%lf", value);
      value = FFFF(dir.getGg(i)) ? 0. : dir.getGg(i);
      _recordWrite("%lf", value);
      _recordWrite("\n");
    }
  }

  /* Close the file */

  _fileClose(verbose);

  return 0;
}

int Vario::fill(int idir,
                int nvar,
                int flagAsym,
                const VectorDouble& sw,
                const VectorDouble& gg,
                const VectorDouble& hh)
{
  return _dirs[idir].fill(nvar, flagAsym, sw, gg, hh);
}

VectorDouble Vario::_getVariableInterval(int ivar) const
{
  VectorDouble bounds(2);
  if (ivar < 0 || ivar >= getVariableNumber())
  {
    bounds[0] = 0;
    bounds[1] = getVariableNumber();
  }
  else
  {
    bounds[0] = ivar;
    bounds[1] = ivar + 1;
  }
  return bounds;
}
VectorDouble Vario::_getDirectionInterval(int idir) const
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
