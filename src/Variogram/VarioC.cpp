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
#include "Variogram/VarioC.hpp"
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

VarioC::VarioC(bool flagSample, VectorDouble dates)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
}

/**
 * Constructing a Variogram by extracting a subset of variables and/or directions
 * from an Input variogram
 *
 * @param VarioC The input variogram
 * @param varcols Vector of variable ranks (starting from 0)
 * @param dircols Vector of direction ranks (starting from 0)
 */
VarioC::VarioC(const VarioC& VarioC,
               const VectorInt& varcols,
               const VectorInt& dircols)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  VectorInt selvars;
  VectorInt seldirs;
  int nvar = VarioC.getVariableNumber();
  int ndir = VarioC._varioparam.getDirectionNumber();

  // Checking arguments
  if (varcols.empty())
    selvars = ut_ivector_sequence(nvar);
  else
  {
    selvars = varcols;
    for (int i = 0; i < (int) varcols.size(); i++)
    {
      if (selvars[i] < 0 || selvars[i] >= nvar)
      {
        messerr(
            "Argument 'varcols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, selvars[i], nvar);
        my_throw("Error when extracting a variogram");
      }
    }
  }

  if (seldirs.empty())
    seldirs = ut_ivector_sequence(ndir);
  else
  {
    seldirs = dircols;
    for (int i = 0; i < ndir; i++)
    {
      if (seldirs[i] < 0 || seldirs[i] >= ndir)
      {
        messerr(
            "Argument 'dircols' is incorrect: item(%d)=%d should lie between 0 and %d",
            i + 1, seldirs[i], ndir);
        my_throw("Error when extracting a variogram");
      }
    }
  }

  // Dimension the new 'VarioC'
  _nVar = (int) selvars.size();

  // Extract the sub-part of VarioParam

  _varioparam = VarioParam(VarioC._varioparam,seldirs);

  // Reset Mean and variance arrays (only if variable number has been modified)
  if (_nVar != nvar)
  {
    _means.resize(_nVar);
    for (int ivar = 0; ivar < _nVar; ivar++)
      setMeans(ivar, VarioC.getMeans(selvars[ivar]));

    _vars.resize(_nVar * _nVar);
    for (int ivar = 0; ivar < _nVar; ivar++)
      for (int jvar = 0; jvar < _nVar; jvar++)
        setVars(ivar, jvar, VarioC.getVars(selvars[ivar], selvars[jvar]));
  }
  else
  {
    _means = VarioC.getMeans();
    _vars = VarioC.getVars();
  }

  // Add the directions

  _sw = std::vector<VectorDouble>();
  _gg = std::vector<VectorDouble>();
  _hh = std::vector<VectorDouble>();
  _utilize = std::vector<VectorDouble>();

  for (int idir = 0; idir < (int) seldirs.size(); idir++)
  {
    _sw.push_back(VarioC.getSwVec(seldirs[idir]));
    _gg.push_back(VarioC.getGgVec(seldirs[idir]));
    _hh.push_back(VarioC.getHhVec(seldirs[idir]));
    _utilize.push_back(VarioC.getUtilizeVec(seldirs[idir]));
  }
}

VarioC::VarioC(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      IClonable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
}

VarioC::VarioC(const VarioC& r)
    : _nVar(r._nVar),
      _varioparam(r._varioparam),
      _means(r._means),
      _vars(r._vars),
      _sw(r._sw),
      _gg(r._gg),
      _hh(r._hh),
      _utilize(r._utilize)
{
}

VarioC& VarioC::operator=(const VarioC& r)
{
  if (this != &r)
  {
    _nVar = r._nVar;
    _varioparam = r._varioparam;
    _means = r._means;
    _vars  = r._vars;
    _sw = r._sw;
    _gg = r._gg;
    _hh = r._hh;
    _utilize = r._utilize;
  }
  return *this;
}

VarioC::~VarioC()
{
}

void VarioC::internalResize(int nvar)
{
  if (nvar <= 0)
  {
    messerr("The Internal Dimension assigned to the variogram cannot be calculated:");
    messerr("- Number of variables = %d",nvar);
    my_throw("Error in Internal Variogram dimensioning");
  }

  // Set the number of variables
  _nVar = nvar;

  // For backwards compatibility, these arrays are updated only if their dimension
  // is not consistent with the current dimension
  if ((int) _means.size() != _nVar) _initMeans();
  if ((int) _vars.size() != _nVar * _nVar) _initVars();

  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    const DirParam dirparam = _varioparam.getDirParam(idir);
    int size = getDirSize(idir);
    _sw[idir].resize(size);
    _gg[idir].resize(size);
    _hh[idir].resize(size);
    _utilize[idir].resize(size);
  }
}

void VarioC::directionResize(int idir)
{
  const DirParam dirparam = _varioparam.getDirParam(idir);
  int size = getDirSize(idir);
  _sw[idir].resize(size);
  _gg[idir].resize(size);
  _hh[idir].resize(size);
  _utilize[idir].resize(size);
}

IClonable* VarioC::clone() const
{
  return new VarioC(*this);
}

double VarioC::getHmax(int ivar, int jvar, int idir) const
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
VectorDouble VarioC::getHRange(int ivar, int jvar, int idir) const
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

double VarioC::getGmax(int ivar,
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

VectorDouble VarioC::getGRange(int ivar,
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

int VarioC::compute(Db *db,
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
  int nvar = db->getVariableNumber();
  internalResize(nvar);

  setMeans(means);
  setVars(vars);

//  return _variogram_compute(db, this, means, vars, flag_grid, flag_gen,
//                            flag_sample, verr_mode, flag_model, model, verbose);
  return 0;
}

/**
 * Compute the Variogram of the Indicator of the input variable
 * It is assumed that the Single input variable is already a Facies
 * (starting from 1)
 */
int VarioC::computeIndic(Db *db,
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

  internalResize(nclass);

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

String VarioC::toString(int level) const
{
  std::stringstream sstr;

  // Print the environment

  sstr << _varioparam.toString(level);

  // Print the variance matrix

  sstr << toMatrix("Variance-Covariance Matrix",VectorString(),VectorString(),
                    0,_nVar,_nVar,getVars());

  if (getCalculType() == CALCUL_UNDEFINED) return sstr.str();

  /* Loop on the directions */

  sstr << std::endl;
  for (int idir=0; idir<getDirectionNumber(); idir++)
  {
    sstr << toTitle(1,"Direction #%d",idir+1);
    const DirParam dirparam = _varioparam.getDirParam(idir);

    /* Print the variogram contents */

    for (int ivar = 0; ivar < _nVar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        sstr << std::endl;
        if (ivar == jvar)
          sstr << "For variable " << ivar + 1 << std::endl;
        else
          sstr << "For variables " << ivar + 1 << " and " << jvar + 1 << std::endl;
        sstr << toStr("Rank");
        sstr << toStr("Npairs");
        sstr << toStr("Distance");
        sstr << toStr("Value");
        sstr << std::endl;

        for (int i = 0; i < dirparam.getLagTotalNumber(); i++)
        {
          int j = _getDirAddress(idir, ivar, jvar, i, true, 0);
          if (_sw[idir][j] <= 0) continue;
          int rank = (! dirparam.getFlagAsym()) ? i : i - dirparam.getNPas();
          sstr << toInt(rank);
          sstr << toDouble(_sw[idir][j]);
          sstr << toDouble(_hh[idir][j]);
          sstr << toDouble(_gg[idir][j]);
          sstr << std::endl;
        }
      }
  }
  return sstr.str();
}

/**
 * Convert the Calculation Name into a Calculation Type (enum)
 * @param calcul_name Input calculation name to be identified
 * @return
 */
int identifyCalculTypeC(const String& calcul_name)
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

double VarioC::getMeans(int ivar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  return _means[ivar];
}

double VarioC::getVars(int ivar, int jvar) const
{
  int iad = _getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return TEST;
  return _vars[iad];
}

double VarioC::getVars(int ijvar) const
{
  if (! _isBivariableValid(ijvar)) return TEST;
  return _vars[ijvar];
}

void VarioC::_initMeans()
{
  _means.resize(_nVar);
  for (int ivar = 0; ivar < _nVar; ivar++)
    _means[ivar] = 0.;
}

void VarioC::setMeans(const VectorDouble& means)
{
  if (_means.empty()) _initMeans();
  if (! means.empty() && (int) means.size() == _nVar)
    _means = means;
}

void VarioC::setMeans(int ivar, double mean)
{
  if (_means.empty()) _initMeans();
  if (! _isVariableValid(ivar)) return;
  _means[ivar] = mean;
}

void VarioC::_initVars()
{
  _vars.resize(_nVar * _nVar);
  int ecr = 0;
  for (int ivar = 0; ivar < _nVar; ivar++)
    for (int jvar = 0; jvar < _nVar; jvar++)
      _vars[ecr++] = (ivar == jvar);
}

void VarioC::setVars(const VectorDouble& vars)
{
  if (_vars.empty()) _initVars();
  if (! vars.empty() && (int) vars.size() == _nVar * _nVar)
    _vars = vars;
}

void VarioC::setVars(int ijvar, double value)
{
  if (_vars.empty()) _initVars();
  if (! _isBivariableValid(ijvar)) return;
  _vars[ijvar] = value;
}

void VarioC::setVars(int ivar, int jvar, double value)
{
  if (_vars.empty()) _initVars();
  int iad = _getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return;
  _vars[iad] = value;
}

double VarioC::getGg(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _gg[idir][i];
}

double VarioC::getHh(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _hh[idir][i];
}

double VarioC::getSw(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _sw[idir][i];
}

double VarioC::getUtilize(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _utilize[idir][i];
}

void VarioC::setGg(int idir, int i, double gg)
{
  if (! _isAddressValid(idir, i)) return;
  _gg[idir][i] = gg;
}

void VarioC::setHh(int idir, int i, double hh)
{
  if (! _isAddressValid(idir, i)) return;
  _hh[idir][i] = hh;
}

void VarioC::setSw(int idir, int i, double sw)
{
  if (! _isAddressValid(idir, i)) return;
  _sw[idir][i] = sw;
}

void VarioC::setUtilize(int idir, int i, double utilize)
{
  if (! _isAddressValid(idir, i)) return;
  _utilize[idir][i] = utilize;
}

double VarioC::getGg(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (! dirparam.isLagValid(ipas)) return TEST;
  int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _gg[idir][iad];
}

double VarioC::getHh(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (! dirparam.isLagValid(ipas)) return TEST;
  int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _hh[idir][iad];
}

double VarioC::getSw(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (! dirparam.isLagValid(ipas)) return TEST;
  int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _sw[idir][iad];
}

double VarioC::getUtilize(int idir, int ivar, int jvar, int ipas) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (! dirparam.isLagValid(ipas)) return TEST;
  int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  return _utilize[idir][iad];
}

VectorDouble VarioC::getGgVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble gg;
  int npas = dirparam.getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] > 0.) gg.push_back(_hh[idir][iad]);
  }
  return gg;
}

VectorDouble VarioC::getHhVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble hh;
  int npas = dirparam.getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] > 0.) hh.push_back(_hh[idir][iad]);
  }
  return hh;
}

VectorDouble VarioC::getSwVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble sw;
  int npas = dirparam.getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    sw.push_back(_sw[idir][iad]);
  }
  return sw;
}

VectorDouble VarioC::getUtilizeVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble utilize;
  int npas = dirparam.getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = _getDirAddress(idir,ivar,jvar,ipas,true,0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] > 0.) utilize.push_back(_hh[idir][iad]);
  }
  return utilize;
}

const VectorDouble& VarioC::getGgVec(int idir) const
{
  return _gg[idir];
}

const VectorDouble& VarioC::getHhVec(int idir) const
{
  return _hh[idir];
}

const VectorDouble& VarioC::getSwVec(int idir) const
{
  return _sw[idir];
}

const VectorDouble& VarioC::getUtilizeVec(int idir) const
{
  return _utilize[idir];
}

int VarioC::getCenter(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return ITEST;
  if (! _varioparam.getDirParam(idir).getFlagAsym()) return ITEST;
  for (int ivar=0; ivar<_nVar; ivar++)
    for (int jvar=0; jvar<=ivar; jvar++)
    {
      int i = _getDirAddress(idir, ivar, jvar, 0, false, 0);
      return i;
    }
  return ITEST;
}

int VarioC::_getVarAddress(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar))  return ITEST;
  if (! _isVariableValid(jvar))  return ITEST;
  return ivar + _nVar * jvar;
}

int VarioC::_getDirAddress(int idir,
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
  if (!dirparam.getFlagAsym())
  {
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
      int npas = dirparam.getNPas();
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
  iad += rank * dirparam.getLagTotalNumber();
  return (iad);
}

bool VarioC::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= _nVar)
  {
    mesArg("Variable Index",ivar,_nVar);
    return false;
  }
  return true;
}

bool VarioC::_isBivariableValid(int ijvar) const
{
  if (ijvar < 0 || ijvar >= _nVar * _nVar)
  {
    mesArg("Multivariate Index",ijvar,_nVar * _nVar);
    return false;
  }
  return true;
}

bool VarioC::_isDirectionValid(int idir) const
{
  if (idir < 0 || idir >= getDirectionNumber())
  {
    mesArg("Direction Index",idir,getDirectionNumber());
    return false;
  }
  return true;
}

bool VarioC::_isAddressValid(int idir, int i) const
{
  if (! _isDirectionValid(idir)) return false;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  if (i < 0 || i >= getDirSize(idir)) return false;
  return true;
}

int VarioC::deSerialize(const String& filename, bool verbose)
{
  int ndim, nvar, ndir, npas, opt_code, flag_calcul, flag_regular;
  double dpas, tolang, scale, tolcode, toldis;
  VectorDouble codir, grloc, vars;
  VectorInt grincr;

  // Open the Neutral File

  if (_fileOpen(filename, "VarioC", "r", verbose)) return 1;

  /* Create the VarioC structure */

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

  internalResize(nvar);
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

    DirParam dirparam = DirParam(ndim);
    for (int idim = 0; idim < ndim; idim++)
      grincr[idim] = (int) grloc[idim];
    dirparam.init(ndim, npas, dpas, toldis, tolang, getFlagAsym(), opt_code, 0,
                  TEST, TEST, tolcode, VectorDouble(), codir, grincr);
    _varioparam.addDirs(dirparam);

    /* Read the arrays of results (optional) */

    if (flag_calcul)
    {
      directionResize(idir);
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

  /* Close the file */
  _fileClose(verbose);

  return 0;
}

int VarioC::serialize(const String& filename, bool verbose) const
{
  double value;
  static int flag_calcul = 1;

  if (_fileOpen(filename, "VarioC", "w", verbose)) return 1;

  /* Write the VarioC structure */

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
    _recordWrite("%d", dirparam.getLagRegular());
    _recordWrite("#", "Regular lags");
    _recordWrite("%d", dirparam.getNPas());
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

void VarioC::patchCenter(int idir, int nech, double rho)
{
  if (! _varioparam.getDirParam(idir).getFlagAsym()) return;
  for (int ivar=0; ivar<_nVar; ivar++)
    for (int jvar=0; jvar<=ivar; jvar++)
    {
      // Get the central address
      int iad = _getDirAddress(idir, ivar, jvar, 0, false, 0);
      if (IFFFF(iad)) continue;
      setSw(idir, iad, (double) nech);
      setHh(idir, iad, 0.);
      if (ivar == jvar)
        setGg(idir, iad, 1.);
      else
        setGg(idir, iad, rho);
    }
}

int VarioC::fill(int idir,
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

VectorDouble VarioC::_getVariableInterval(int ivar) const
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

VectorDouble VarioC::_getDirectionInterval(int idir) const
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
int VarioC::getDirSize(int idir) const
{
  const DirParam dirparam = _varioparam.getDirParam(idir);
  return (dirparam.getLagTotalNumber() * _nVar * (_nVar + 1) / 2);
}

/**
 * Convert the Calculation Name into a Calculation Type (enum)
 * @return
 */

int VarioC::getCalculType() const
{
  return identifyCalculTypeC(_varioparam.getCalculName());
}

int VarioC::getFlagAsym() const
{
  return identifyFlagAsymC(_varioparam.getCalculName());
}

int identifyFlagAsymC(const String& calcul_name)
{
  int flagAsym = 0;

  switch (identifyCalculTypeC(calcul_name))
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

