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
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_f_private.h"

#include "Enum/EAnam.hpp"

#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Stats/Classical.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"

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
      ICloneable(),
      _nVar(0),
      _varioparam(),
      _means(means),
      _vars(vars),
      _calcul(),
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

Vario::Vario(const Vario& r)
    : AStringable(r),
      ASerializable(r),
      _nVar(r._nVar),
      _varioparam(r._varioparam),
      _means(r._means),
      _vars(r._vars),
      _calcul(r._calcul),
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
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _nVar = r._nVar;
    _varioparam = r._varioparam;
    _means = r._means;
    _vars  = r._vars;
    _calcul = r._calcul;
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

Vario* Vario::create(const VarioParam* varioparam,
                     Db* db,
                     const VectorDouble& means,
                     const VectorDouble& vars)
{
  return new Vario(varioparam, db, means, vars);
}

Vario* Vario::createFromNF(const String& neutralFilename, bool verbose)
{
  Vario* vario = nullptr;
  std::ifstream is;
  VarioParam* varioparam = new VarioParam();
  vario = new Vario(varioparam);
  bool success = false;
  if (vario->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  vario->deserialize(is, verbose);
  }
  if (! success)
  {
    delete vario;
    delete varioparam;
    vario = nullptr;
  }
  return vario;
}

Vario* Vario::computeFromDb(const VarioParam* varioparam,
                            Db* db,
                            const ECalcVario& calcul,
                            bool flag_gen,
                            bool flag_sample,
                            bool verr_mode,
                            Model *model,
                            bool verbose)
{
  Vario* vario = nullptr;
  vario = new Vario(varioparam, db);
  if (vario->compute(calcul, flag_gen, flag_sample, verr_mode, model,
                 verbose))
  {
    return nullptr;
  }
  return vario;
}

Vario* Vario::createRegularizeFromModel(const Model* model,
                                        const VarioParam* varioparam,
                                        const VectorDouble& ext,
                                        const VectorInt& ndisc,
                                        const VectorDouble& angles)
{
  Vario* vario = nullptr;
  vario = new Vario(varioparam, nullptr);
  if (vario->modelRegularize(model, ext, ndisc, angles))
  {
    messerr("Error when calculating the regularized variogram");
    return nullptr;
  }
  return vario;
}

Vario* Vario::createTransformZToY(const Vario* varioZ,
                                  const AAnam* anam,
                                  double cvv)
{
  Vario* varioY = varioZ->clone();
  if (varioY->transformZToY(anam, cvv))
  {
    messerr("Error when transforming Raw Variogram into Gaussian");
    return nullptr;
  }
  return varioY;
}

Vario* Vario::createTransformYToZ(const Vario* varioY,
                                  const AAnam* anam,
                                  const Model* model)
{
  Vario* varioZ = varioY->clone();
  if (varioZ->transformYToZ(anam, model))
  {
    messerr("Error when transforming Gaussian Variogram into Raw");
    return nullptr;
  }
  return varioZ;
}

Vario* Vario::createReduce(const Vario *varioIn,
                           const VectorInt &varcols,
                           const VectorInt &dircols,
                           bool asSymmetric)
{
  Vario* varioOut = varioIn->clone();
  varioOut->reduce(varcols, dircols, asSymmetric);
  return varioOut;
}

int Vario::compute(const ECalcVario &calcul,
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

  _calcul = calcul;
  _setFlagAsym();
  _setDPasFromGrid(isDefinedForGrid());
  if (internalVariableResize()) return 1;
  internalDirectionResize();

  if (_variogram_compute(_db, this, flag_gen, flag_sample, verr_mode, model,
                         verbose))
  {
    messerr("Error when calculating the Variogram");
    return 1;
  }
  return 0;
}

int Vario::computeByKey(const String& calcul_name,
                        bool flag_gen,
                        bool flag_sample,
                        bool verr_mode,
                        Model *model,
                        bool verbose)
{
  ECalcVario calcul = getCalculType(calcul_name);
  if (calcul == ECalcVario::UNDEFINED) return 1;
  return compute(calcul, flag_gen, flag_sample, verr_mode, model, verbose);
}

/**
 * Reduce this current variogram by keeping a subset of variables and/or directions
 *
 * @param varcols Vector of variable ranks (starting from 0)
 * @param dircols Vector of direction ranks (starting from 0)
 * @param asSymmetric Turn the result into as Symmetrical function (i.e. variogram)
 */
void Vario::reduce(const VectorInt& varcols,
                        const VectorInt& dircols,
                        bool asSymmetric)
{
  VectorInt selvars;
  VectorInt seldirs;
  Vario vario_in(*this); // Copy this as input variogram
  int nvar_in = vario_in.getVariableNumber();
  int ndir_in = vario_in._varioparam.getDirectionNumber();

  // Checking arguments
  if (varcols.empty())
    selvars = VH::sequence(nvar_in);
  else
  {
    selvars = varcols;
    for (int i = 0; i < static_cast<int>(varcols.size()); i++)
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
  _nVar = static_cast<int>(selvars.size());
  if (_nVar <= 0)
  {
    my_throw("The number of variables extracted cannot be zero");
  }

  if (seldirs.empty())
    seldirs = VH::sequence(ndir_in);
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
  int ndir = static_cast<int>(seldirs.size());
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
        setMean(ivar, vario_in.getMean(selvars[ivar]));
    }

    if (! vario_in.getVars().empty())
    {
      _vars.resize(_nVar * _nVar);
      for (int ivar = 0; ivar < _nVar; ivar++)
        for (int jvar = 0; jvar < _nVar; jvar++)
          setVar(ivar, jvar, vario_in.getVar(selvars[ivar], selvars[jvar]));
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
            _sw[idir][iadto] = vario_in.getSwByIndex(idir0,iadfrom);
            _gg[idir][iadto] = vario_in.getGgByIndex(idir0,iadfrom);
            _hh[idir][iadto] = vario_in.getHhByIndex(idir0,iadfrom);
            _utilize[idir][iadto] = vario_in.getUtilizeByIndex(idir0,iadfrom);
          }
          else
          {
            int iadf1 = vario_in.getDirAddress(idir,ivar,jvar,ipas,false,-1);
            int iadf2 = vario_in.getDirAddress(idir,ivar,jvar,ipas,false,1);
            _sw[idir][iadto] = (vario_in.getSwByIndex(idir0, iadf1)
                + vario_in.getSwByIndex(idir0, iadf2)) / 2.;
            _gg[idir][iadto] = (vario_in.getGgByIndex(idir0, iadf1)
                + vario_in.getGgByIndex(idir0, iadf2)) / 2.;
            _hh[idir][iadto] = (ABS(vario_in.getHhByIndex(idir0, iadf1))
                + ABS(vario_in.getHhByIndex(idir0, iadf2))) / 2.;
            _utilize[idir][iadto] = (vario_in.getUtilizeByIndex(idir0, iadf1)
                + vario_in.getUtilizeByIndex(idir0, iadf2)) / 2.;
            if (flagMakeSym)
            {
              double c0 = vario_in.getVar(ivar,jvar);
              _gg[idir][iadto] = c0 - _gg[idir][iadto];
            }
          }
        }
    }
  }
}

int Vario::computeIndic(const ECalcVario& calcul,
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
  int nclass = static_cast<int>(props.size());
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
  int iatt = _db->getUIDByLocator(ELoc::Z, 0);
  if (limits.toIndicatorByAttribute(_db, iatt))
  {
    messerr("Problem when translating Facies into Categories");
    return 1;
  }

  // Preparation

  _calcul = calcul;
  _nVar  = nclass;
  _means = props;
  _vars  = _varsFromProportions(props);
  _setFlagAsym();
  _setDPasFromGrid(isDefinedForGrid());
  if (internalVariableResize()) return 1;
  internalDirectionResize();

  // Calculate the variogram of indicators
  if (_variogram_compute(_db, this, flag_gen,
                         flag_sample, verr_mode, model, verbose))
  {
    messerr("Error when calculating the Variogram of Indicators");
    return 1;
  }

  // Delete the Indicators (created locally)
  _db->deleteColumnsByLocator(ELoc::Z);
  _db->setLocatorByUID(iatt, ELoc::Z);

  return 0;
}

int Vario::computeIndicByKey(const String& calcul_name,
                             bool flag_gen,
                             bool flag_sample,
                             bool verr_mode,
                             Model *model,
                             bool verbose,
                             int nfacmax)
{
  ECalcVario calcul = getCalculType(calcul_name);
  if (calcul == ECalcVario::UNDEFINED) return 1;
  return computeIndic(calcul, flag_gen, flag_sample, verr_mode,
                      model, verbose, nfacmax);
}

/*****************************************************************************/
/*!
 **  Transform the experimental variogram from raw to gaussian space
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        Point Hermite anamorphosis
 ** \param[in]  cvv         Block variance
 **
 *****************************************************************************/
int Vario::transformZToY(const AAnam *anam, double cvv)
{
  if (anam == nullptr)
  {
    messerr("This function needs an Anamorphosis");
    return 1;
  }
  AnamHermite* anamH = dynamic_cast<AnamHermite*>(anam->clone());
  if (anamH == nullptr)
  {
    messerr("This function needs a Hermite Anamorphosis");
    return 1;
  }
  if (getVariableNumber() != 1)
  {
    messerr("This function is restricted to Monovariate Variogram");
    return 1;
  }

  /* Loop on the directions of the variogram */

  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    /* Loop on the lags */

    for (int i = 0; i < getLagNumber(idir); i++)
    {
      // TODO. GG must be a variogram of Zv -> Cv(h)
      setGgByIndex(idir,i,1. - anamH->invertVariance(cvv-getGgByIndex(idir, i)));
    }
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate the experimental variogram of the Raw starting from the Model
 **  of the Gaussian variable
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        Point anamorphosis
 ** \param[in]  model       Model of the Punctual Gaussian
 **
 *****************************************************************************/
int Vario::transformYToZ(const AAnam *anam, const Model *model)
{
  CovCalcMode mode;

  /* Preliminary checks */

  int ndim = getDimensionNumber();
  if (anam == (AAnam*) NULL) return 1;
  if (model == nullptr) return 1;
  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("This function is restricted to Gaussian Anamorphosis");
    return 1;
  }
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam->clone());
  if (anam_hermite == nullptr)
  {
    messerr("This function needs a Hermite Anamorphosis");
    return 1;
  }
  if (anam_hermite->getRCoef() != 1.)
  {
    messerr("This function is restricted to Punctual Anamorphosis");
    return 1;
  }
  if (getVariableNumber() != 1)
  {
    messerr("This function is restricted to Monovariate Variogram");
    return 1;
  }
  if (model->getVariableNumber() != 1)
  {
    messerr("This function requires a Monovariate Model");
    return 1;
  }
  if (model->getDimensionNumber() != ndim)
  {
    messerr("Variogram and Model should share the same Space Dimension");
    return 1;
  }

  /* Calculate the theoretical variance of Z */

  double varz = anam_hermite->computeVariance(1.);

  /* Loop on the directions of the variogram */

  VectorDouble shift(ndim);
  SpacePoint p1(model->getContext().getSpace());
  SpacePoint p2(model->getContext().getSpace());
  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    /* Loop on the lags */

    for (int ipas = 0; ipas < getLagNumber(idir); ipas++)
    {
      for (int idim = 0; idim < ndim; idim++)
        shift[idim] = (ipas + 1) * getDPas(idir) * getCodir(idir, idim);

      p2 = p1;
      p2.move(shift);
      double chh = model->eval(0,0,p1,p2,mode);
      if (chh < 0.)
      {
        messerr("Gaussian covariance is negative in direction %d for lag %d",
                idir + 1, ipas + 1);
        messerr("Calculation is impossible");
        return 1;
      }

      double cov = anam_hermite->computeVariance(chh);
      setGg(idir, 0, 0, ipas, varz - cov);
      setHh(idir, 0, 0, ipas, (ipas + 1) * getDPas(idir));
      setSw(idir, 0, 0, ipas, 1.);
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the regularized model as an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  ext       Vector of Block extension
 ** \param[in]  ndisc     Vector of discretization counts
 ** \param[in]  angles    Vector of rotation angles (optional)
 ** \param[in]  mode      CovCalcMode structure
 **
 *****************************************************************************/
int Vario::modelRegularize(const Model* model,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles,
                           const CovCalcMode& mode)
{
  if (model == nullptr)
  {
    messerr("Model must be provided");
    return 1;
  }
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();

  /* Preliminary checks */

  setNVar(nvar);
  internalVariableResize();
  internalDirectionResize();

  /* Initialize the variance array */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = model->evalCvv(ext, ndisc, angles, ivar, jvar, mode);
      setVar(ivar, jvar, value);
    }

  /* Loop on the directions */

  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {

    /* Loop on the number of lags */

    for (int ipas = 0; ipas < getLagNumber(idir); ipas++)
    {

      // Calculate the shift vector

      double dist = ipas * getDPas(idir);
      VectorDouble shift(ndim);
      for (int idim = 0; idim < ndim; idim++)
        shift[idim] = dist * getCodir(idir, idim);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          double value = model->evalCvvShift(ext, ndisc, shift, angles, ivar, jvar, mode);
          int iad = getDirAddress(idir, ivar, jvar, ipas, false, 0);
          setGgByIndex(idir, iad, value);
          setHhByIndex(idir, iad, dist);
          setSwByIndex(idir, iad, 1);
        }
    }
  }
  return 0;
}

int Vario::attachDb(Db* db, const VectorDouble& vars, const VectorDouble& means)
{
  _db = db;
  if (db != nullptr)
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

int Vario::internalVariableResize() //TODO: to be called when nvar is modified...
{
  if (! _means.empty())
  {
    int nloc = static_cast<int>(_means.size());
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
    int nloc = static_cast<int>(_vars.size());
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

double Vario::getHmax(int ivar, int jvar, int idir) const
{
  VectorInt ivb = _getVariableInterval(ivar);
  VectorInt jvb = _getVariableInterval(jvar);
  VectorInt idb = _getDirectionInterval(idir);

  double hmax = 0.;
  for (int id = idb[0]; id < idb[1]; id++)
    for (int iv = ivb[0]; iv < ivb[1]; iv++)
      for (int jv = jvb[0]; jv < jvb[1]; jv++)
      {
        VectorDouble hh = getHhVec(id, iv, jv);
        double hhloc = VH::maximum(hh);
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
  VectorInt ivb = _getVariableInterval(ivar);
  VectorInt jvb = _getVariableInterval(jvar);
  VectorInt idb = _getDirectionInterval(idir);

  VectorDouble vec(2);
  vec[0] =  1.e30;
  vec[1] = -1.e30;
  for (int id = idb[0]; id < idb[1]; id++)
    for (int iv = ivb[0]; iv < ivb[1]; iv++)
      for (int jv = jvb[0]; jv < jvb[1]; jv++)
      {
        VectorDouble hh = getHhVec(id, iv, jv);
        double hhmin = VH::minimum(hh);
        double hhmax = VH::maximum(hh);
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
  VectorInt ivb = _getVariableInterval(ivar);
  VectorInt jvb = _getVariableInterval(jvar);
  VectorInt idb = _getDirectionInterval(idir);

  double gmax = 0.;

  for (int id = idb[0]; id < idb[1]; id++)
    for (int iv = ivb[0]; iv < ivb[1]; iv++)
      for (int jv = jvb[0]; jv < jvb[1]; jv++)
      {
        VectorDouble gg = getGgVec(id, iv, jv);
        double ggloc = VH::maximum(gg, flagAbs);
        if (ggloc > gmax) gmax = ggloc;
        if (flagSill)
        {
          double sill = ABS(getVar(iv, jv));
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
  VectorInt ivb = _getVariableInterval(ivar);
  VectorInt jvb = _getVariableInterval(jvar);
  VectorInt idb = _getDirectionInterval(idir);

  VectorDouble vec(2);
  vec[0] =  1.e30;
  vec[1] = -1.e30;

  for (int id = idb[0]; id < idb[1]; id++)
    for (int iv = ivb[0]; iv < ivb[1]; iv++)
      for (int jv = jvb[0]; jv < jvb[1]; jv++)
      {
        VectorDouble gg = getGgVec(id, iv, jv);
        double ggmin = VH::minimum(gg);
        double ggmax = VH::maximum(gg);
        if (ggmin < vec[0]) vec[0] = ggmin;
        if (ggmax > vec[1]) vec[1] = ggmax;
        if (flagSill)
        {
          double sill = getVar(iv, jv);
          if (sill < vec[0]) vec[0] = sill;
          if (sill > vec[1]) vec[1] = sill;
        }
      }
  return vec;
}

String Vario::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  // Print the calculation type

  switch (getCalcul().toEnum())
  {
    case ECalcVario::E_UNDEFINED:
      sstr << toTitle(0,"Undefined");
      break;

    case ECalcVario::E_VARIOGRAM:
      sstr << toTitle(0,"Variogram characteristics");
      break;

    case ECalcVario::E_MADOGRAM:
      sstr << toTitle(0,"Madogram characteristics");
      break;

    case ECalcVario::E_RODOGRAM:
      sstr << toTitle(0,"Rodogram characteristics");
      break;

    case ECalcVario::E_POISSON:
      sstr << toTitle(0,"Poisson variogram characteristics");
      break;

    case ECalcVario::E_COVARIANCE:
      sstr << toTitle(0,"Covariance characteristics");
      break;

    case ECalcVario::E_COVARIANCE_NC:
      sstr << toTitle(0,"Non-centered Covariance characteristics");
      break;

    case ECalcVario::E_COVARIOGRAM:
      sstr << toTitle(0,"Transitive Covariogram characteristics");
      break;

    case ECalcVario::E_GENERAL1:
      sstr << toTitle(0,"Generalized Variogram of order 1 characteristics");
      break;

    case ECalcVario::E_GENERAL2:
      sstr << toTitle(0,"Generalized Variogram of order 2 characteristics");
      break;

    case ECalcVario::E_GENERAL3:
      sstr << toTitle(0,"Generalized Variogram of order 3 characteristics");
      break;

    case ECalcVario::E_ORDER4:
      sstr << toTitle(0,"Order-4 Variogram");
      break;

    case ECalcVario::E_TRANS1:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/G1");
      break;

    case ECalcVario::E_TRANS2:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/G2");
      break;

    case ECalcVario::E_BINORMAL:
      sstr << toTitle(0,"Cross-to_simple Variogram ratio G12/sqrt(G1*G2)");
      break;

    default:
      break;
  }
  if (getCalcul() == ECalcVario::UNDEFINED) return sstr.str();
  sstr << "Number of variable(s)       = " << _nVar << std::endl;

  // Print the environment

  sstr << _varioparam.toStringMain(strfmt);

  // Print the variance matrix

  sstr << toMatrix("Variance-Covariance Matrix",VectorString(),VectorString(),
                    0,_nVar,_nVar,getVars());

  if (getCalcul() == ECalcVario::UNDEFINED) return sstr.str();

  /* Loop on the directions (only if the resulting arrays have been defined) */

  if (!_sw.empty())
  {
    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {
      sstr << toTitle(1,"Direction #%d",idir+1);
      sstr << getDirParam(idir).toString(strfmt);
      sstr << _toStringByDirection(strfmt,idir);
    }
  }
  return sstr.str();
}

String Vario::_toStringByDirection(const AStringFormat* /*strfmt*/, int idir) const
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
 * Convert the Calculation Name into a Calculation Type (ECalcVario)
 *
 * @return The corresponding ECalcVario enum
 */
ECalcVario Vario::getCalculType(const String& calcul_name) const
{
  ECalcVario calcul_type;

  if (calcul_name == "undefined")
    calcul_type = ECalcVario::UNDEFINED;
  else if (calcul_name == "vg")
    calcul_type = ECalcVario::VARIOGRAM;
  else if (calcul_name == "cov")
    calcul_type = ECalcVario::COVARIANCE;
  else if (calcul_name == "covnc")
    calcul_type = ECalcVario::COVARIANCE_NC;
  else if (calcul_name == "covg")
    calcul_type = ECalcVario::COVARIOGRAM;
  else if (calcul_name =="mado")
    calcul_type = ECalcVario::MADOGRAM;
  else if (calcul_name =="rodo")
    calcul_type = ECalcVario::RODOGRAM;
  else if (calcul_name =="poisson")
    calcul_type = ECalcVario::POISSON;
  else if (calcul_name =="general1")
    calcul_type = ECalcVario::GENERAL1;
  else if (calcul_name =="general2")
    calcul_type = ECalcVario::GENERAL2;
  else if (calcul_name =="general3")
    calcul_type = ECalcVario::GENERAL3;
  else if (calcul_name =="order4")
    calcul_type = ECalcVario::ORDER4;
  else if (calcul_name =="trans1")
    calcul_type = ECalcVario::TRANS1;
  else if (calcul_name =="trans2")
    calcul_type = ECalcVario::TRANS2;
  else if (calcul_name =="binormal")
    calcul_type = ECalcVario::BINORMAL;
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

    calcul_type = ECalcVario::UNDEFINED;
  }
  return calcul_type;
}

double Vario::getMean(int ivar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  return _means[ivar];
}

double Vario::getVar(int ivar, int jvar) const
{
  int iad = getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return TEST;
  return _vars[iad];
}

double Vario::getVarIndex(int ijvar) const
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
  if (! means.empty() && static_cast<int>(means.size()) == _nVar)
    _means = means;
}

void Vario::setMean(int ivar, double mean)
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
  if (! vars.empty() && static_cast<int>(vars.size()) == _nVar * _nVar)
    _vars = vars;
}

void Vario::setVarIndex(int ijvar, double value)
{
  if (_vars.empty()) _initVars();
  if (! _isBivariableValid(ijvar)) return;
  _vars[ijvar] = value;
}

void Vario::setVar(int ivar, int jvar, double value)
{
  if (_vars.empty()) _initVars();
  int iad = getVarAddress(ivar, jvar);
  if (IFFFF(iad)) return;
  _vars[iad] = value;
}

double Vario::getGgByIndex(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _gg[idir][i];
}

double Vario::getHhByIndex(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _hh[idir][i];
}

double Vario::getSwByIndex(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _sw[idir][i];
}

double Vario::getUtilizeByIndex(int idir, int i) const
{
  if (! _isAddressValid(idir, i)) return(TEST);
  return _utilize[idir][i];
}

void Vario::setGgByIndex(int idir, int i, double gg)
{
  if (! _isAddressValid(idir, i)) return;
  _gg[idir][i] = gg;
}

void Vario::setHhByIndex(int idir, int i, double hh)
{
  if (! _isAddressValid(idir, i)) return;
  _hh[idir][i] = hh;
}

void Vario::setSwByIndex(int idir, int i, double sw)
{
  if (! _isAddressValid(idir, i)) return;
  _sw[idir][i] = sw;
}

void Vario::setUtilizeByIndex(int idir, int i, double utilize)
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

void Vario::updateSwByIndex(int idir, int i, double sw)
{
  if (! _isAddressValid(idir, i)) return;
  _sw[idir][i] += sw;
}

void Vario::updateHhByIndex(int idir, int i, double hh)
{
  if (! _isAddressValid(idir, i)) return;
  _hh[idir][i] += hh;
}

void Vario::updateGgByIndex(int idir, int i, double gg)
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
    double c0 = getVar(ivar, jvar);
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

/**
 * Returns a triple array of values: 0-Weight; 1-Distance; 2-Variogram
 * @param idir Target Direction
 * @param ivar Target first variable
 * @param jvar Target second variable
 * @return
 */
VectorVectorDouble Vario::getVec(int idir, int ivar, int jvar) const
{
  VectorVectorDouble vec;
  if (!_isVariableValid(ivar)) return vec;
  if (!_isVariableValid(jvar)) return vec;
  if (!_isDirectionValid(idir)) return vec;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  int npas = dirparam.getLagNumber();
  vec.resize(3);
  for (int i = 0; i < 3; i++) vec[i].resize(npas);

  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
    vec[0][ipas] = _sw[idir][iad];
    vec[1][ipas] = _hh[idir][iad];
    vec[2][ipas] = _gg[idir][iad];
  }
  return vec;
}

/**
 * Returns the vector of variogram values for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @param asCov True if result should be provided as Covariance, False for Variogram
 * @param flagNorm If the result should be provided as a Normalized covariance / Variogram
 * @return The vector of 'gg' (not calculated lags are suppressed)
 */
VectorDouble Vario::getGgVec(int idir,
                             int ivar,
                             int jvar,
                             bool asCov,
                             bool flagNorm) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);

  VectorDouble gg;
  double c0 = 0.;
  if (asCov || flagNorm) c0 = getVar(ivar, jvar);
  int npas = dirparam.getLagNumber();

  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas-1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] <= 0.) continue;
      double val = _gg[idir][iad];
      if (asCov && !getFlagAsym()) val = c0 - val;
      if (flagNorm) val /= c0;
      gg.push_back(val);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (! IFFFF(iad) && _sw[idir][iad] > 0.)
    {
      double val = _gg[idir][iad];
      if (asCov && !getFlagAsym()) val = c0 - val;
      if (flagNorm) val /= c0;
      gg.push_back(val);
    }
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (! IFFFF(iad) && _sw[idir][iad] > 0.)
      {
        double val = _gg[idir][iad];
        if (asCov && !getFlagAsym()) val = c0 - val;
        if (flagNorm) val /= c0;
        gg.push_back(val);
      }
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (! IFFFF(iad) && _sw[idir][iad] > 0.)
      {
        double val = _gg[idir][iad];
        if (asCov && !getFlagAsym()) val = c0 - val;
        if (flagNorm) val /= c0;
        gg.push_back(val);
      }
    }
  }
  return gg;
}

/**
 * Returns the vector of distances for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @return The vector of 'hh' (not calculated lags are suppressed)
 */
VectorDouble Vario::getHhVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);

  VectorDouble hh;
  int npas = dirparam.getLagNumber();
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] > 0) hh.push_back(_hh[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    double value = _hh[idir][iad];
    if (_sw[idir][iad] <= 0) value = 0.;
    hh.push_back(value);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] > 0.) hh.push_back(_hh[idir][iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] > 0.) hh.push_back(_hh[idir][iad]);
    }
  }
  return hh;
}

/**
 * Returns the vector of weights for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @return The vector of 'sw' (not calculated lags are suppressed)
 */
VectorDouble Vario::getSwVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar))  return VectorDouble();
  if (!_isVariableValid(jvar))  return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);

  VectorDouble sw;
  int npas = dirparam.getLagNumber();
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (!IFFFF(iad)) sw.push_back(_sw[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (!IFFFF(iad)) sw.push_back(_sw[idir][iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (!IFFFF(iad)) sw.push_back(_sw[idir][iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] > 0.) sw.push_back(_sw[idir][iad]);
    }
  }
  return sw;
}

/**
 * Returns the vector of utilization flags for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @return The vector of 'utilize' (not calculated lags are suppressed)
 */
VectorDouble Vario::getUtilizeVec(int idir, int ivar, int jvar) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();
  const DirParam dirparam = _varioparam.getDirParam(idir);

  VectorDouble utilize;
  int npas = dirparam.getLagNumber();
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
     {
       iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
       if (!IFFFF(iad)) utilize.push_back(_utilize[idir][iad]);
     }
     iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
     if (!IFFFF(iad)) utilize.push_back(_utilize[idir][iad]);
     for (int ipas = 0; ipas < npas; ipas++)
     {
       iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
       if (!IFFFF(iad)) utilize.push_back(_utilize[idir][iad]);
     }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (IFFFF(iad)) continue;
      if (_sw[idir][iad] > 0.) utilize.push_back(_utilize[idir][iad]);
    }
  }
  return utilize;
}

const VectorDouble& Vario::getAllGg(int idir) const
{
  return _gg[idir];
}

const VectorDouble& Vario::getAllHh(int idir) const
{
  return _hh[idir];
}

const VectorDouble& Vario::getAllSw(int idir) const
{
  return _sw[idir];
}

const VectorDouble& Vario::getAllUtilize(int idir) const
{
  return _utilize[idir];
}

int Vario::getCenter(int ivar, int jvar, int idir) const
{
  if (! _isDirectionValid(idir)) return ITEST;
  if (! _isVariableValid(ivar))  return ITEST;
  if (! _isVariableValid(jvar))  return ITEST;
  if (! getFlagAsym()) return ITEST;
  int i = getDirAddress(idir, ivar, jvar, 0, false, 0);
  return i;
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
  if (!_isVariableValid(ivar))  return ITEST;
  if (!_isVariableValid(jvar))  return ITEST;

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
    if (! dirparam.isLagValid(ipas, getFlagAsym())) return ITEST;
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
      if (! dirparam.isLagValid(ipas, getFlagAsym())) return ITEST;
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

bool Vario::_deserialize(std::istream& is, bool /*verbose*/)
{
  int flag_calcul = 0;
  int ndim = 0;
  int nvar = 0;
  int ndir = 0;
  int npas = 0;
  int opt_code = 0;
  int flag_regular = 0;

  double dpas = 0.;
  double tolang = 0.;
  double scale = 0.;
  double tolcode = 0.;
  double toldis = 0.;

  VectorDouble vars;

  /* Create the Vario structure */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Number of Variables", nvar);
  ret = ret && _recordRead<int>(is, "Number of Variogram Directions", ndir);
  ret = ret && _recordRead<double>(is, "Scale", scale);

  /* Read the variances (optional) */

  ret = ret && _recordRead<int>(is, "Variogram calculation Option", flag_calcul);
  vars.resize(nvar * nvar);
  if (flag_calcul)
  {
    int ecr = 0;
    for (int ivar = 0; ret && ivar < nvar; ivar++)
      for (int jvar = 0; ret && jvar < nvar; jvar++, ecr++)
        ret = ret && _recordRead<double>(is, "Experimental Variance term", vars[ecr]);
  }
  if (! ret) return ret;

  /* Initialize the variogram structure */

  _nVar = nvar;
  internalDirectionResize(ndir,false);
  setVars(vars);
  setCalculName("vg");
  setScale(scale);
  int isDefinedForGrid = 0;

  /* Reading the variogram calculation directions */

  for (int idir = 0; ret && idir < ndir; idir++)
  {
    ret = ret && _recordRead<int>(is, "Regular Variogram Calculation", flag_regular);
    ret = ret && _recordRead<int>(is, "Number of Variogram Lags", npas);
    ret = ret && _recordRead<int>(is, "Variogram Code Option", opt_code);
    ret = ret && _recordRead<double>(is, "Tolerance on Code", tolcode);
    ret = ret && _recordRead<double>(is, "Lag Value", dpas);
    ret = ret && _recordRead<double>(is, "Tolerance on Distance", toldis);
    ret = ret && _recordRead<int>(is, "Grid Definition",isDefinedForGrid);

    VectorDouble codir;
    VectorInt grincr;
    if (! isDefinedForGrid)
    {

      // Direction definition

      ret = ret && _recordRead<double>(is, "Tolerance on Direction", tolang);
      ret = ret && _recordReadVec<double>(is, "Direction vector", codir, ndim);
    }
    else
    {
      // Grid definition

      ret = ret && _recordReadVec<int>(is, "Grid Increment", grincr, ndim);
    }
    if (! ret) return ret;

    SpaceRN space(ndim);
    DirParam dirparam = DirParam(npas, dpas, toldis, tolang, opt_code, 0,
                                 TEST, TEST, tolcode, VectorDouble(), codir, grincr,
                                 &space);
    _varioparam.addDir(dirparam);

    /* Read the arrays of results (optional) */

    if (flag_calcul)
    {
      _directionResize(idir);
      for (int i = 0; ret && i < getDirSize(idir); i++)
      {
        double sw = 0.;
        double hh = 0.;
        double gg = 0.;
        ret = ret && _recordRead<double>(is, "Experimental Variogram Weight", sw);
        setSwByIndex(idir, i, sw);
        ret = ret && _recordRead<double>(is, "Experimental Variogram Distance", hh);
        setHhByIndex(idir, i, hh);
        ret = ret && _recordRead<double>(is, "Experimental Variogram Value", gg);
        setGgByIndex(idir, i, gg);
      }
    }
  }
  return ret;
}

bool Vario::_serialize(std::ostream& os, bool /*verbose*/) const
{
  double value;
  static int flag_calcul = 1;

  /* Write the Vario structure */

  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", _varioparam.getDimensionNumber());
  ret = ret && _recordWrite<int>(os, "Number of variables", getVariableNumber());
  ret = ret && _recordWrite<int>(os, "Number of directions", getDirectionNumber());
  ret = ret && _recordWrite<double>(os, "Scale", _varioparam.getScale());
  ret = ret && _recordWrite<int>(os, "Calculation Flag", flag_calcul);

  /* Dumping the Variances */

  if (flag_calcul)
  {
    ret = ret && _commentWrite(os, "Variance");
    for (int ivar = 0; ret && ivar < getVariableNumber(); ivar++)
    {
      for (int jvar = 0; ret && jvar < getVariableNumber(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getVar(ivar,jvar));
      ret = ret && _commentWrite(os, "");
    }
  }

  /* Loop on the directions */

  for (int idir = 0; ret && idir < getDirectionNumber(); idir++)
  {
    const DirParam dirparam = _varioparam.getDirParam(idir);
    ret = ret && _commentWrite(os, "Direction characteristics");
    ret = ret && _recordWrite<int>(os, "Regular lags", dirparam.getFlagRegular());
    ret = ret && _recordWrite<int>(os, "Number of lags", dirparam.getLagNumber());
    ret = ret && _recordWrite<int>(os, "", dirparam.getOptionCode());
    ret = ret && _recordWrite<double>(os, "Code selection: Option - Tolerance", dirparam.getTolCode());
    ret = ret && _recordWrite<double>(os, "Lag value", dirparam.getDPas());
    ret = ret && _recordWrite<double>(os, "Tolerance on distance", dirparam.getTolDist());
    ret = ret && _recordWrite<int>(os, "Grid Definition",dirparam.isDefinedForGrid());

    if (! dirparam.isDefinedForGrid())
    {
      ret = ret && _recordWrite<double>(os, "Tolerance on angle", dirparam.getTolAngle());
      for (int idim = 0; idim < (int) dirparam.getNDim() && ret; idim++)
        ret = ret && _recordWrite<double>(os, "", dirparam.getCodir(idim));
      ret = ret && _commentWrite(os, "Direction coefficients");
    }

    if (dirparam.isDefinedForGrid())
    {
      for (int idim = 0; ret && idim < (int) dirparam.getNDim() && ret; idim++)
        ret = ret && _recordWrite(os, "", (double) dirparam.getGrincr(idim));
      ret = ret && _commentWrite(os, "Direction increments on grid");
    }

    if (!flag_calcul) continue;
    ret = ret && _commentWrite(os, "Variogram results (Weight, Distance, Variogram)");
    for (int i = 0; i < getDirSize(idir) && ret; i++)
    {
      value = FFFF(getSwByIndex(idir, i)) ? 0. : getSwByIndex(idir, i);
      ret = ret && _recordWrite<double>(os, "", value);
      value = FFFF(getHhByIndex(idir, i)) ? 0. : getHhByIndex(idir, i);
      ret = ret && _recordWrite<double>(os, "", value);
      value = FFFF(getGgByIndex(idir, i)) ? 0. : getGgByIndex(idir, i);
      ret = ret && _recordWrite<double>(os, "", value);
      ret = ret && _commentWrite(os, "");
    }
  }
  return ret;
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
      setSwByIndex(idir, iad, (double) nech);
      setHhByIndex(idir, iad, 0.);
      if (ivar == jvar)
        setGgByIndex(idir, iad, 1.);
      else
        setGgByIndex(idir, iad, rho);
    }
}

int Vario::fill(int idir,
                 const VectorDouble& sw,
                 const VectorDouble& gg,
                 const VectorDouble& hh)
{
  if (! _isDirectionValid(idir)) return 1;
  int size = getDirSize(idir);
  if (size != static_cast<int>(sw.size()) ||
      size != static_cast<int>(hh.size()) ||
      size != static_cast<int>(gg.size()))
  {
    messerr("The argument do not have correct dimension");
    return 1;
  }
  for (int i=0; i<size; i++)
  {
    setSwByIndex(idir, i, sw[i]);
    setHhByIndex(idir, i, hh[i]);
    setGgByIndex(idir, i, gg[i]);
  }
  return 0;
}

VectorInt Vario::_getVariableInterval(int ivar) const
{
  VectorInt bounds(2);
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

VectorInt Vario::_getDirectionInterval(int idir) const
{
  VectorInt bounds(2);
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
  switch (getCalcul().toEnum())
  {
    case ECalcVario::E_VARIOGRAM:
    case ECalcVario::E_MADOGRAM:
    case ECalcVario::E_RODOGRAM:
    case ECalcVario::E_POISSON:
    case ECalcVario::E_GENERAL1:
    case ECalcVario::E_GENERAL2:
    case ECalcVario::E_GENERAL3:
    case ECalcVario::E_ORDER4:
    case ECalcVario::E_TRANS1:
    case ECalcVario::E_TRANS2:
    case ECalcVario::E_BINORMAL:
      _flagAsym = false;
      break;

    case ECalcVario::E_COVARIANCE:
    case ECalcVario::E_COVARIANCE_NC:
    case ECalcVario::E_COVARIOGRAM:
      _flagAsym = true;
      break;
    default:
      messerr("Wrong Variogram Calculation enum!");
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
    _nVar = static_cast<int>(_means.size());
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
  _calcul = getCalculType(calcul_name);
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
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(_db);
    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {
      _varioparam.setDPas(idir, dbgrid);
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

VectorDouble Vario::_varsFromProportions(VectorDouble props)
{
  if (props.empty()) return VectorDouble();

  int nvar = static_cast<int>(props.size());
  VectorDouble vars = VectorDouble(nvar * nvar);
  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      vars[ecr++] = (ivar == jvar) ?
          props[ivar] * (1. - props[ivar]) : - props[ivar] * props[jvar];
    }
  return vars;
}

bool Vario::drawOnlyPositiveX(int ivar, int jvar) const
{
  if (ivar == jvar || ! getFlagAsym()) return true;
  return false;
}

bool Vario::drawOnlyPositiveY(int ivar, int jvar) const
{
  if (ivar == jvar && ! getFlagAsym()) return true;
  return false;
}
