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
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_f.h"
#include "geoslib_define.h"
#include "geoslib_f_private.h"

/**
 * Build a Vario object by calculating the experimental variogram
 * @param varioparam VarioParam structure
 * @param db         Db structure (optional)
 * @param means      Array of variable means
 * @param vars       Array of variable variances
 */
Vario::Vario(const VarioParam* varioparam,
             Db* db,
             const VectorDouble& means,
             const VectorDouble& vars)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(means),
      _vars(vars),
      _calculName( ),
      _flagSample( ),
      _db(nullptr),
      _sw(),
      _gg(),
      _hh(),
      _utilize(),
      _flagAsym(false)
{
  _varioparam = *varioparam;
  attachDb(db,vars,means);
}

/**
 * Constructing a Variogram by extracting a subset of variables and/or directions
 * from an Input variogram
 *
 * @param vario_in The input variogram
 * @param varcols Vector of variable ranks (starting from 0)
 * @param dircols Vector of direction ranks (starting from 0)
 * @param asSymmetric Turn the result into as Symmetrical function (i.e. variogram)
 */
Vario::Vario(const Vario& vario_in,
             const VectorInt& varcols,
             const VectorInt& dircols,
             bool asSymmetric)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _calculName(vario_in._calculName),
      _flagSample(vario_in._flagSample),
      _db(nullptr),
      _sw(),
      _gg(),
      _hh(),
      _utilize(),
      _flagAsym(false)
{
  VectorInt selvars;
  VectorInt seldirs;
  int nvar_in = vario_in.getVariableNumber();
  int ndir_in = vario_in._varioparam.getDirectionNumber();

  // Checking arguments
  if (varcols.empty())
    selvars = ut_ivector_sequence(nvar_in);
  else
  {
    selvars = varcols;
    for (int i = 0; i < (int) varcols.size(); i++)
    {
      if (selvars[i] < 0 || selvars[i] >= nvar_in)
      {
        messerr(
            "Argument 'varcols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, selvars[i], nvar_in);
        my_throw("Error when extracting a variogram");
      }
    }
  }
  _nVar = (int) selvars.size();
  if (_nVar <= 0)
  {
    my_throw("The number of variables extracted cannot be zero");
  }

  if (seldirs.empty())
    seldirs = ut_ivector_sequence(ndir_in);
  else
  {
    seldirs = dircols;
    for (int i = 0; i < ndir_in; i++)
    {
      if (seldirs[i] < 0 || seldirs[i] >= ndir_in)
      {
        messerr(
            "Argument 'dircols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, seldirs[i], ndir_in);
        my_throw("Error when extracting a variogram");
      }
    }
  }
  int ndir = (int) seldirs.size();
  if (ndir <= 0)
  {
    my_throw("The number of directions extracted cannot be zero");
  }

  // Extract the sub-part of VarioParam

  _varioparam = VarioParam(vario_in._varioparam,seldirs);
  bool flagMakeSym = false;
  if (asSymmetric)
  {
    if (vario_in.getFlagAsym()) flagMakeSym = true;
    setCalculName("vg");
  }

  // Reset Mean and variance arrays (only if variable number has been modified)
  if (_nVar != nvar_in)
  {
    if (! vario_in.getMeans().empty())
    {
      _means.resize(_nVar);
      for (int ivar = 0; ivar < _nVar; ivar++)
        setMeans(ivar, vario_in.getMeans(selvars[ivar]));
    }

    if (! vario_in.getVars().empty())
    {
      _vars.resize(_nVar * _nVar);
      for (int ivar = 0; ivar < _nVar; ivar++)
        for (int jvar = 0; jvar < _nVar; jvar++)
          setVars(ivar, jvar, vario_in.getVars(selvars[ivar], selvars[jvar]));
    }
  }
  else
  {
    _means = vario_in.getMeans();
    _vars = vario_in.getVars();
  }

  // Add the directions

  internalDirectionResize(ndir_in,true);

  for (int idir0 = 0; idir0 < ndir; idir0++)
  {
    int idir = seldirs[idir0];
    for (int ipas = 0; ipas < getLagNumber(idir); ipas++)
    {
      for (int ivar0 = 0; ivar0 < _nVar; ivar0++)
        for (int jvar0 = 0; jvar0 < _nVar; jvar0++)
        {
          int ivar = selvars[ivar0];
          int jvar = selvars[jvar0];

          int iadto = getDirAddress(idir0,ivar0,jvar0,ipas);

          if (! flagMakeSym)
          {
            int iadfrom = vario_in.getDirAddress(idir,ivar,jvar,ipas,false,0);
            _sw[idir][iadto] = vario_in.getSw(idir0,iadfrom);
            _gg[idir][iadto] = vario_in.getGg(idir0,iadfrom);
            _hh[idir][iadto] = vario_in.getHh(idir0,iadfrom);
            _utilize[idir][iadto] = vario_in.getUtilize(idir0,iadfrom);
          }
          else
          {
            int iadf1 = vario_in.getDirAddress(idir,ivar,jvar,ipas,false,-1);
            int iadf2 = vario_in.getDirAddress(idir,ivar,jvar,ipas,false,1);
            _sw[idir][iadto] = (vario_in.getSw(idir0, iadf1)
                + vario_in.getSw(idir0, iadf2)) / 2.;
            _gg[idir][iadto] = (vario_in.getGg(idir0, iadf1)
                + vario_in.getGg(idir0, iadf2)) / 2.;
            _hh[idir][iadto] = (ABS(vario_in.getHh(idir0, iadf1))
                + ABS(vario_in.getHh(idir0, iadf2))) / 2.;
            _utilize[idir][iadto] = (vario_in.getUtilize(idir0, iadf1)
                + vario_in.getUtilize(idir0, iadf2)) / 2.;
            if (flagMakeSym)
            {
              double c0 = vario_in.getVars(ivar,jvar);
              _gg[idir][iadto] = c0 - _gg[idir][iadto];
            }
          }
        }
    }
  }
}

Vario::Vario(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _calculName("undefined"),
      _flagSample(false),
      _db(nullptr),
      _sw(),
      _gg(),
      _hh(),
      _utilize(),
      _flagAsym(false)
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
}

Vario::Vario(const Vario& r)
    : _nVar(r._nVar),
      _varioparam(r._varioparam),
      _means(r._means),
      _vars(r._vars),
      _calculName(r._calculName),
      _flagSample(r._flagSample),
      _db(r._db),
      _sw(r._sw),
      _gg(r._gg),
      _hh(r._hh),
      _utilize(r._utilize),
      _flagAsym(r._flagAsym)
{
}

Vario& Vario::operator=(const Vario& r)
{
  if (this != &r)
  {
    _nVar = r._nVar;
    _varioparam = r._varioparam;
    _means = r._means;
    _vars  = r._vars;
    _calculName = r._calculName;
    _flagSample = r._flagSample;
    _db = r._db;
    _sw = r._sw;
    _gg = r._gg;
    _hh = r._hh;
    _utilize = r._utilize;
    _flagAsym = r._flagAsym;
  }
  return *this;
}

Vario::~Vario()
{
}

int Vario::compute(const String& calcul_name,
                   bool flag_grid,
                   bool flag_gen,
                   bool flag_sample,
                   bool verr_mode,
                   Model *model,
                   bool verbose)
{
  if (_db == nullptr)
  {
    messerr("The 'Db' must have been attached beforehand");
    return 1;
  }
  _nVar = _db->getVariableNumber();
  if (_nVar <= 0)
  {
    messerr("The 'db' must contain at least one variable defined");
    return 1;
  }

  // Preparation

  setCalculName(calcul_name);
  _setDPasFromGrid(flag_grid);
  if (internalVariableResize()) return 1;
  internalDirectionResize();

  if (_variogram_compute(_db, this, _means, _vars, flag_grid, flag_gen,
                         flag_sample, verr_mode, model, verbose))
  {
    messerr("Error when calculating the Variogram");
    return 1;
  }
  return 0;
}

int Vario::computeIndic(const String& calcul_name,
                        bool flag_grid,
                        bool flag_gen,
                        bool flag_sample,
                        bool verr_mode,
                        Model *model,
                        bool verbose,
                        int nfacmax)
{
  // Preliminary checks
  if (_db == nullptr)
  {
    messerr("The 'Db' must have been attached beforehand");
    return 1;
  }
  int nvar = _db->getVariableNumber();
  if (nvar != 1)
  {
    messerr("The 'db' must contain ONE variable defined");
    return 1;
  }

  // Calculate the number of Facies in 'Db'
  VectorDouble props = dbStatisticsFacies(_db);
  int nclass = props.size();
  if (nclass <= 0 || (nfacmax > 0 && nclass > nfacmax))
  {
    messerr("The input variable should exhibit Facies");
    messerr(
        "Number of Facies (%d) should be positive and smaller than 'nfacmax'",
        nclass);
    messerr("Note: the value of 'nfacmax'(%d) can be changed in argument list",
            nfacmax);
    return 1;
  }

  // Translate the 'Facies' into 'categories'   VectorDouble props =
  Limits limits = Limits(nclass);
  int iatt = _db->getAttribute(LOC_Z, 0);
  if (limits.toIndicator(_db, iatt))
  {
    messerr("Problem when translating Facies into Categories");
    return 1;
  }

  // Preparation

  setCalculName(calcul_name);
  _nVar  = nclass;
  _means = props;
  _vars  = VectorDouble();
  _setDPasFromGrid(flag_grid);
  if (internalVariableResize()) return 1;
  internalDirectionResize();

  // Calculate the variogram of indicators
  if (_variogram_compute(_db, this, _means, _vars, flag_grid, flag_gen,
                         flag_sample, verr_mode, model, verbose))
  {
    messerr("Error when calculating the Variogram of Indicators");
    return 1;
  }

  // Delete the Indicators (created locally)
  _db->deleteFieldByLocator(LOC_Z);
  _db->setLocatorByAttribute(iatt, LOC_Z);

  return 0;
}

int Vario::attachDb(Db* db, const VectorDouble& vars, const VectorDouble& means)
{
  _db = db;
  if (db != (Db *) NULL)
  {
    int nvar = _db->getVariableNumber();
    if (nvar <= 0)
    {
      messerr("Some variables must be defined in the 'Db'");
      return 1;
    }
  }
  _means = means;
  _vars  = vars;
  return 0;
}

int Vario::internalVariableResize()
{
  if (! _means.empty())
  {
    int nloc = _means.size();
    if (nloc != _nVar)
    {
      messerr("Invalid dimension for 'means' (%d)",nloc);
      messerr("It should match the number of variables in 'Db' (%d)",_nVar);
      return 1;
    }
  }
  else
  {
    _initMeans();
  }

  if (! _vars.empty())
  {
    int nloc = _vars.size();
    if (nloc != _nVar * _nVar)
    {
      messerr("Invalid dimension for 'vars' (%d)",nloc);
      messerr("It should match the number of variables in 'Db' (squared) (%d)",
              _nVar * _nVar);
      return 1;
    }
  }
  else
  {
    _initVars();
  }
  return 0;
}

void Vario::internalDirectionResize(int ndir, bool flagDirs)
{
  if (ndir <= 0) ndir = getDirectionNumber();
  _sw.resize(ndir);
  _gg.resize(ndir);
  _hh.resize(ndir);
  _utilize.resize(ndir);

  if (flagDirs)
    for (int idir = 0; idir < getDirectionNumber(); idir++)
      _directionResize(idir);
}

void Vario::_directionResize(int idir)
{
  int size = getDirSize(idir);
  _sw[idir].resize(size);
  _gg[idir].resize(size);
  _hh[idir].resize(size);
  _utilize[idir].resize(size,1.); // By default, all lags are usable
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
        VectorDouble hh = getHhVec(idir, ivar, jvar);
        double hhloc = ut_vector_max(hh);
        if (hhloc > hmax) hmax = hhloc;
      }
  return hmax;
}

/**
 * Returns a vector with Minimum-Maximum of the Hh array
 * @param ivar Target variable (or -1)
 * @param jvar Target variable (or -1)
 * @param idir Target Direction (or -1)
 * @return
 */
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
        VectorDouble hh = getHhVec(idir, ivar, jvar);
        double hhmin = ut_vector_min(hh);
        double hhmax = ut_vector_max(hh);
        if (hhmin < vec[0]) vec[0] = hhmin;
        if (hhmax > vec[1]) vec[1] = hhmax;
      }
  return vec;
}

double Vario::getGmax(int ivar,
                       int jvar,
                       int idir,
                       bool flagAbs,
                       bool flagSill) const
{
  VectorDouble ivb = _getVariableInterval(ivar);
  VectorDouble jvb = _getVariableInterval(jvar);
  VectorDouble idb = _getDirectionInterval(idir);

  double gmax = 0.;

  for (int idir = idb[0]; idir < idb[1]; idir++)
    for (int ivar = ivb[0]; ivar < ivb[1]; ivar++)
      for (int jvar = jvb[0]; jvar < jvb[1]; jvar++)
      {
        VectorDouble gg = getGgVec(idir, ivar, jvar);
        double ggloc = ut_vector_max(gg, flagAbs);
        if (ggloc > gmax) gmax = ggloc;
        if (flagSill)
        {
          double sill = ABS(getVars(ivar, jvar));
          if (gmax < sill) gmax = sill;
        }
      }
  return gmax;
}

VectorDouble Vario::getGRange(int ivar,
                               int jvar,
                               int idir,
                               bool flagSill) const
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
        VectorDouble gg = getGgVec(idir, ivar, jvar);
        double ggmin = ut_vector_min(gg);
        double ggmax = ut_vector_max(gg);
        if (ggmin < vec[0]) vec[0] = ggmin;
        if (ggmax > vec[1]) vec[1] = ggmax;
        if (flagSill)
        {
          double sill = getVars(ivar, jvar);
          if (sill < vec[0]) vec[0] = sill;
          if (sill > vec[1]) vec[1] = sill;
        }
      }
  return vec;
}

String Vario::toString(int level) const
{
  std::stringstream sstr;

  // Print the calculation type

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
  if (getCalculType() == CALCUL_UNDEFINED) return sstr.str();
  sstr << "Number of variable(s)       = " << _nVar << std::endl;

  // Print the environment

  sstr << _varioparam.toStringMain(level);

  // Print the variance matrix

  sstr << toMatrix("Variance-Covariance Matrix",VectorString(),VectorString(),
                    0,_nVar,_nVar,getVars());

  if (getCalculType() == CALCUL_UNDEFINED) return sstr.str();

  /* Loop on the directions (only if the resulting arrays have been defined) */

  if (!_sw.empty())
  {
    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {
      sstr << toTitle(1,"Direction #%d",idir+1);
      sstr << getDirParam(idir).toString(level);
      sstr << _toStringByDirection(level,idir);
    }
  }
  return sstr.str();
}

String Vario::_toStringByDirection(int level, int idir) const
{
  std::stringstream sstr;

  const DirParam dirparam = _varioparam.getDirParam(idir);

  /* Print the variogram contents */

  for (int ivar = 0; ivar < _nVar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      sstr << std::endl;
      if (ivar == jvar)
        sstr << "For variable " << ivar + 1 << std::endl;
      else
        sstr << "For variables " << ivar + 1 << " and " << jvar + 1
             << std::endl;
      sstr << toStr("Rank");
      sstr << toStr("Npairs");
      sstr << toStr("Distance");
      sstr << toStr("Value");
      sstr << std::endl;

      for (int i = 0; i < getLagTotalNumber(idir); i++)
      {
        int j = getDirAddress(idir, ivar, jvar, i, true, 0);
        if (_sw[idir][j] <= 0) continue;
        int rank = (!getFlagAsym()) ? i : i - dirparam.getLagNumber();
        sstr << toInt(rank);
        sstr << toDouble(_sw[idir][j]);
        sstr << toDouble(_hh[idir][j]);
        sstr << toDouble(_gg[idir][j]);
        sstr << std::endl;
      }
    }
  return sstr.str();
}

/**
 * Convert the Calculation Name into a Calculation Type (enum)
 * @return
 */
int Vario::getCalculType() const
{
  int calcul_type;

  if (_calculName == "undefined")
    calcul_type = CALCUL_UNDEFINED;
  else if (_calculName == "vg")
    calcul_type = CALCUL_VARIOGRAM;
  else if (_calculName == "cov")
    calcul_type = CALCUL_COVARIANCE;
  else if (_calculName == "covnc")
    calcul_type = CALCUL_COVARIANCE_NC;
  else if (_calculName == "covg")
    calcul_type = CALCUL_COVARIOGRAM;
  else if (_calculName =="mado")
    calcul_type = CALCUL_MADOGRAM;
  else if (_calculName =="rodo")
    calcul_type = CALCUL_RODOGRAM;
  else if (_calculName =="poisson")
    calcul_type = CALCUL_POISSON;
  else if (_calculName =="general1")
    calcul_type = CALCUL_GENERAL1;
  else if (_calculName =="general2")
    calcul_type = CALCUL_GENERAL2;
  else if (_calculName =="general3")
    calcul_type = CALCUL_GENERAL3;
  else if (_calculName =="order4")
    calcul_type = CALCUL_ORDER4;
  else if (_calculName =="trans1")
    calcul_type = CALCUL_TRANS1;
  else if (_calculName =="trans2")
    calcul_type = CALCUL_TRANS2;
  else if (_calculName =="binormal")
    calcul_type = CALCUL_BINORMAL;
  else
  {
    messerr("Invalid variogram calculation name : %s", _calculName.c_str());
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

double Vario::getMeans(int ivar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  return _means[ivar];
}

double Vario::getVars(int ivar, int jvar) const
{
  int iad = getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return TEST;
  return _vars[iad];
}

double Vario::getVars(int ijvar) const
{
  if (! _isBivariableValid(ijvar)) return TEST;
  return _vars[ijvar];
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

void Vario::setVars(int ijvar, double value)
{
  if (_vars.empty()) _initVars();
  if (! _isBivariableValid(ijvar)) return;
  _vars[ijvar] = value;
}

void Vario::setVars(int ivar, int jvar, double value)
{
  if (_vars.empty()) _initVars();
  int iad = getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return;
  _vars[iad] = value;
}

double Vario::getGg(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _gg[idir][i];
}

double Vario::getHh(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _hh[idir][i];
}

double Vario::getSw(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _sw[idir][i];
}

double Vario::getUtilize(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _utilize[idir][i];
}

void Vario::setGg(int idir, int i, double gg)
{
  if (! _isAddressValid(idir, i)) return;
  _gg[idir][i] = gg;
}

void Vario::setHh(int idir, int i, double hh)
{
  if (! _isAddressValid(idir, i)) return;
  _hh[idir][i] = hh;
}

void Vario::setSw(int idir, int i, double sw)
{
  if (! _isAddressValid(idir, i)) return;
  _sw[idir][i] = sw;
}

void Vario::setUtilize(int idir, int i, double utilize)
{
  if (! _isAddressValid(idir, i)) return;
  _utilize[idir][i] = utilize;
}

void Vario::setSw(int idir, int ivar, int jvar, int ipas, double sw)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return;
  _sw[idir][iad] = sw;
}

void Vario::setHh(int idir, int ivar, int jvar, int ipas, double hh)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return;
  _hh[idir][iad] = hh;
}

void Vario::setGg(int idir, int ivar, int jvar, int ipas, double gg)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return;
  _gg[idir][iad] = gg;
}

void Vario::setUtilize(int idir, int ivar, int jvar, int ipas, double utilize)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return;
  _utilize[idir][iad] = utilize;
}

void Vario::updSw(int idir, int i, double sw)
{
  if (! _isAddressValid(idir, i)) return;
  _sw[idir][i] += sw;
}

void Vario::updHh(int idir, int i, double hh)
{
  if (! _isAddressValid(idir, i)) return;
  _hh[idir][i] += hh;
}

void Vario::updGg(int idir, int i, double gg)
{
  if (! _isAddressValid(idir, i)) return;
  _gg[idir][i] += gg;
}

double Vario::getGg(int idir,
                    int ivar,
                    int jvar,
                    int ipas,
                    bool flagCov,
                    bool flagNorm) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  double val = _gg[idir][iad];
  if (flagCov || flagNorm)
  {
    double c0 = getVars(ivar, jvar);
    if (flagCov && ! getFlagAsym())  val = c0 - val;
    if (flagNorm) val /= c0;
  }
  return val;
}

double Vario::getHh(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _hh[idir][iad];
}

double Vario::getSw(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _sw[idir][iad];
}

double Vario::getUtilize(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _utilize[idir][iad];
}

VectorDouble Vario::getGgVec(int idir,
                             int ivar,
                             int jvar,
                             bool asCov,
                             bool flagNorm) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble gg;
  double c0 = 0.;
  if (asCov || flagNorm) c0 = getVars(ivar, jvar);
  int npas = dirparam.getLagNumber();

  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] <= 0.) continue;
    double val = _gg[idir][iad];
    if (asCov && ! getFlagAsym())  val = c0 - val;
    if (flagNorm) val /= c0;
    gg.push_back(val);
  }
  return gg;
}

VectorDouble Vario::getHhVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble hh;
  int npas = dirparam.getLagNumber();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] > 0.) hh.push_back(_hh[idir][iad]);
  }
  return hh;
}

VectorDouble Vario::getSwVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble sw;
  int npas = dirparam.getLagNumber();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    sw.push_back(_sw[idir][iad]);
  }
  return sw;
}

VectorDouble Vario::getUtilizeVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble utilize;
  int npas = dirparam.getLagNumber();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] > 0.) utilize.push_back(_hh[idir][iad]);
  }
  return utilize;
}

const VectorDouble& Vario::getGgVec(int idir) const
{
  return _gg[idir];
}

const VectorDouble& Vario::getHhVec(int idir) const
{
  return _hh[idir];
}

const VectorDouble& Vario::getSwVec(int idir) const
{
  return _sw[idir];
}

const VectorDouble& Vario::getUtilizeVec(int idir) const
{
  return _utilize[idir];
}

int Vario::getCenter(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return ITEST;
  if (! getFlagAsym()) return ITEST;
  for (int ivar=0; ivar<_nVar; ivar++)
    for (int jvar=0; jvar<=ivar; jvar++)
    {
      int i = getDirAddress(idir, ivar, jvar, 0, false, 0);
      return i;
    }
  return ITEST;
}

int Vario::getVarAddress(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar))  return ITEST;
  if (! _isVariableValid(jvar))  return ITEST;
  return ivar + _nVar * jvar;
}

int Vario::getDirAddress(int idir,
                         int ivar,
                         int jvar,
                         int ipas,
                         bool flag_abs,
                         int sens) const
{
  if (!_isDirectionValid(idir)) return ITEST;
  if (!_isVariableValid(ivar)) return ITEST;
  if (!_isVariableValid(jvar)) return ITEST;

  int rank;

  /* Get the order of the variables */

  if (ivar > jvar)
    rank = ivar * (ivar + 1) / 2 + jvar;
  else
    rank = jvar * (jvar + 1) / 2 + ivar;

  /* Get the position in the array */

  int iad = 0;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  if (! getFlagAsym())
  {
    if (! dirparam.isLagValid(ipas)) return ITEST;
    iad = ipas;
  }
  else
  {
    if (flag_abs)
    {
      iad = ipas;
    }
    else
    {
      if (! dirparam.isLagValid(ipas)) return ITEST;
      int npas = dirparam.getLagNumber();
      switch (sens)
      {
        case 1:
          iad = npas + ipas + 1;
          break;

        case -1:
          iad = npas - ipas - 1;
          break;

        case 0:
          iad = npas;
          break;
      }
    }
  }
  iad += rank * getLagTotalNumber(idir);
  return (iad);
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

bool Vario::_isBivariableValid(int ijvar) const
{
  if (ijvar < 0 || ijvar >= _nVar * _nVar)
  {
    mesArg("Multivariate Index",ijvar,_nVar * _nVar);
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

bool Vario::_isAddressValid(int idir, int i) const
{
  if (! _isDirectionValid(idir)) return false;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (i < 0 || i >= getDirSize(idir)) return false;
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
        if (_recordRead("Experimental Variance term", "%lf", &vars[ecr]))
          goto label_end;
  }

  /* Initialize the variogram structure */

  _nVar = nvar;
  internalDirectionResize(ndir,false);
  setVars(vars);
  setCalculName("vg");
  setScale(scale);

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

    DirParam dirparam = DirParam(ndim);
    for (int idim = 0; idim < ndim; idim++)
      grincr[idim] = (int) grloc[idim];
    dirparam.init(ndim, npas, dpas, toldis, tolang, opt_code, 0,
                  TEST, TEST, tolcode, VectorDouble(), codir, grincr);
    _varioparam.addDirs(dirparam);

    /* Read the arrays of results (optional) */

    if (flag_calcul)
    {
      _directionResize(idir);
      for (int i = 0; i < getDirSize(idir); i++)
      {
        double sw, hh, gg;
        if (_recordRead("Experimental Variogram Weight", "%lf", &sw)) goto label_end;
        setSw(idir, i, sw);
        if (_recordRead("Experimental Variogram Distance", "%lf", &hh)) goto label_end;
        setHh(idir, i, hh);
        if (_recordRead("Experimental Variogram Value", "%lf", &gg)) goto label_end;
        setGg(idir, i, gg);
      }
    }
  }

  label_end:

  _fileClose(verbose);

  return 0;
}

int Vario::serialize(const String& filename, bool verbose) const
{
  double value;
  static int flag_calcul = 1;

  if (_fileOpen(filename, "Vario", "w", verbose)) return 1;

  /* Write the Vario structure */

  _recordWrite("%d", _varioparam.getDimensionNumber());
  _recordWrite("#", "Space Dimension");
  _recordWrite("%d", getVariableNumber());
  _recordWrite("#", "Number of variables");
  _recordWrite("%d", getDirectionNumber());
  _recordWrite("#", "Number of directions");
  _recordWrite("%lf", _varioparam.getScale());
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
    const DirParam dirparam = _varioparam.getDirParam(idir);
    _recordWrite("#", "Direction characteristics");
    _recordWrite("%d", dirparam.getFlagRegular());
    _recordWrite("#", "Regular lags");
    _recordWrite("%d", dirparam.getLagNumber());
    _recordWrite("#", "Number of lags");
    _recordWrite("%d", dirparam.getOptionCode());
    _recordWrite("%lf", dirparam.getTolCode());
    _recordWrite("#", "Code selection: Option - Tolerance");
    _recordWrite("%lf", dirparam.getDPas());
    _recordWrite("#", "Lag value");
    _recordWrite("%lf", dirparam.getTolDist());
    _recordWrite("#", "Tolerance on distance");
    _recordWrite("%lf", dirparam.getTolAngle());
    _recordWrite("#", "Tolerance on angle");

    for (int idim = 0; idim < dirparam.getDimensionNumber(); idim++)
      _recordWrite("%lf", dirparam.getCodir(idim));
    _recordWrite("#", "Direction coefficients");

    for (int idim = 0; idim < dirparam.getDimensionNumber(); idim++)
      _recordWrite("%lf", (double) dirparam.getGrincr(idim));
    _recordWrite("#", "Direction increments on grid");

    if (!flag_calcul) continue;
    _recordWrite("#", "Variogram results (Weight, Distance, Variogram)");
    for (int i = 0; i < getDirSize(idir); i++)
    {
      value = FFFF(getSw(idir, i)) ? 0. : getSw(idir, i);
      _recordWrite("%lf", value);
      value = FFFF(getHh(idir, i)) ? 0. : getHh(idir, i);
      _recordWrite("%lf", value);
      value = FFFF(getGg(idir, i)) ? 0. : getGg(idir, i);
      _recordWrite("%lf", value);
      _recordWrite("\n");
    }
  }

  /* Close the file */

  _fileClose(verbose);

  return 0;
}

void Vario::patchCenter(int idir, int nech, double rho)
{
  if (! getFlagAsym()) return;
  for (int ivar=0; ivar<_nVar; ivar++)
    for (int jvar=0; jvar<=ivar; jvar++)
    {
      // Get the central address
      int iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
      if (IFFFF(iad)) continue;
      setSw(idir, iad, (double) nech);
      setHh(idir, iad, 0.);
      if (ivar == jvar)
        setGg(idir, iad, 1.);
      else
        setGg(idir, iad, rho);
    }
}

int Vario::fill(int idir,
                 const VectorDouble& sw,
                 const VectorDouble& gg,
                 const VectorDouble& hh)
{
  if (! _isDirectionValid(idir)) return 1;
  int size = getDirSize(idir);
  if (size != (int) sw.size() ||
      size != (int) hh.size() ||
      size != (int) gg.size())
  {
    messerr("The argument do not have correct dimension");
    return 1;
  }
  for (int i=0; i<(int) size; i++)
  {
    setSw(idir, i, sw[i]);
    setHh(idir, i, hh[i]);
    setGg(idir, i, gg[i]);
  }
  return 0;
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
  if (idir < 0 || idir >= getDirectionNumber())
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

/**
 * Returns the Dimension of internal arrays
 * after the number of variables has been defined
 * @param idir Target direction
 * @return
 */
int Vario::getDirSize(int idir) const
{
  return (getLagTotalNumber(idir) * _nVar * (_nVar + 1) / 2);
}

/**
 * Get the Asymmetrical flag
 * @return
 */

void Vario::_setFlagAsym()
{
  switch (getCalculType())
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
      _flagAsym = false;
      break;

    case CALCUL_COVARIANCE:
    case CALCUL_COVARIANCE_NC:
    case CALCUL_COVARIOGRAM:
      _flagAsym = true;
      break;
  }
}

/**
 * Derive the number of variables from arguments
 * @param db Db structure (optional)
 * @return Error return code
 */
int Vario::_getNVar(const Db* db)
{
  if (db != nullptr)
  {
    _nVar = db->getVariableNumber();
    return 0;
  }
  else if (!_means.empty())
  {
    _nVar = _means.size();
    return 0;
  }

  messerr("Cannot determine the Number of Variables from arguments");
  return 1;
}

int Vario::getLagTotalNumber(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  int npas = getDirParam(idir).getLagNumber();
  return ((_flagAsym) ? 2 * npas + 1 : npas);
}

void Vario::setCalculName(const String calcul_name)
{
  _calculName = calcul_name;
  _setFlagAsym();
}

/**
 * In the case of Grid, set the lag value according to the Db grid parameters
 * @param flag_grid Flag for using the Grid organization
 */
void Vario::_setDPasFromGrid(bool flag_grid)
{
  if (_db->isGrid() && flag_grid)
  {
    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {
      _varioparam.setDPas(idir, _db);
    }
  }
  else
  {
    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {
      _varioparam.setGrincr(idir, VectorInt());
    }
  }
}
