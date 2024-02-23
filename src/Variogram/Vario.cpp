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
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_f_private.h"
#include "geoslib_f.h"

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
#include "Basic/OptDbg.hpp"
#include "Stats/Classical.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "Geometry/BiTargetCheckCode.hpp"
#include "Geometry/BiTargetCheckDate.hpp"
#include "Geometry/BiTargetCheckFaults.hpp"
#include "Geometry/BiTargetCheckGeometry.hpp"

static Model *MODEL;
static Vario *VARIO;
static int IECH1, IECH2, IDIRLOC;

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int VARWGT[4][5] = { { 1, -1, 0, 0, 0 },
                            { 1, -2, 1, 0, 0 },
                            { 1, -3, 3, -1, 0 },
                            { 1, -4, 6, -4, 1 } };


/**
 * Build a Vario object by calculating the experimental variogram
 * @param varioparam VarioParam structure
 */
Vario::Vario(const VarioParam& varioparam)
    : AStringable(),
      ASerializable(),
      ICloneable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _calcul(),
      _flagSample(false),
      _db(nullptr),
      _sw(),
      _gg(),
      _hh(),
      _utilize(),
      _biPtsPerDirection(0),
      _bipts(),
      _flagAsym(false)
{
  _varioparam = varioparam;
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
      _biPtsPerDirection(r._biPtsPerDirection),
      _flagAsym(r._flagAsym)
{
  for (int ipt = 0, npt = _getBiPtsNumber(); ipt < npt; ipt++)
    _bipts.push_back(r._bipts[ipt]);
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
    _biPtsPerDirection = r._biPtsPerDirection;
    _flagAsym = r._flagAsym;

    for (int ipt = 0, npt = _getBiPtsNumber(); ipt < npt; ipt++)
      _bipts.push_back(r._bipts[ipt]);
  }
  return *this;
}

Vario::~Vario()
{
}

Vario* Vario::create(const VarioParam& varioparam)
{
  return new Vario(varioparam);
}

Vario* Vario::createSkeleton(const VarioParam &varioparam,
                             int nvar,
                             const ECalcVario &calcul,
                             const VectorDouble &means,
                             const VectorDouble &vars)
{
  Vario* vario = new Vario(varioparam);
  vario->setNVar(nvar);
  vario->setMeans(means);
  vario->setVars(vars);

  if (vario->prepare(calcul, false)) return nullptr;
  return vario;
}

Vario* Vario::createFromNF(const String& neutralFilename, bool verbose)
{
  Vario* vario = nullptr;
  std::ifstream is;
  VarioParam* varioparam = new VarioParam();
  vario = new Vario(*varioparam);
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

Vario* Vario::computeFromDb(const VarioParam& varioparam,
                            Db* db,
                            const ECalcVario& calcul,
                            bool flag_gen,
                            bool flag_sample,
                            bool verr_mode,
                            Model *model,
                            bool verbose)
{
  Vario* vario = nullptr;
  vario = new Vario(varioparam);
  if (vario->compute(db, calcul, flag_gen, flag_sample, verr_mode, model, verbose))
  {
    return nullptr;
  }
  return vario;
}

Vario* Vario::createRegularizeFromModel(const Model& model,
                                        const VarioParam& varioparam,
                                        const VectorDouble& ext,
                                        const VectorInt& ndisc,
                                        const VectorDouble& angles,
                                        bool asCov)
{
  DECLARE_UNUSED(asCov);
  Vario* vario = nullptr;
  vario = new Vario(varioparam);
  if (vario->modelRegularize(model, ext, ndisc, angles))
  {
    messerr("Error when calculating the regularized variogram");
    return nullptr;
  }
  return vario;
}

Vario* Vario::createTransformZToY(const Vario& varioZ,
                                  const AAnam* anam)
{
  Vario* varioY = new Vario(varioZ);
  if (varioY->transformZToY(anam))
  {
    messerr("Error when transforming Raw Variogram into Gaussian");
    return nullptr;
  }
  return varioY;
}

Vario* Vario::createTransformYToZ(const Vario& varioY,
                                  const AAnam* anam)
{
  Vario* varioZ = new Vario(varioY);
  if (varioZ->transformYToZ(anam))
  {
    messerr("Error when transforming Gaussian Variogram into Raw");
    return nullptr;
  }
  return varioZ;
}

Vario* Vario::createReduce(const Vario& varioIn,
                           const VectorInt &varcols,
                           const VectorInt &dircols,
                           bool asSymmetric)
{
  Vario* varioOut = new Vario(varioIn);
  varioOut->resetReduce(varcols, dircols, asSymmetric);
  return varioOut;
}

void Vario::_clearBiTargetCheck()
{
  for (int ipt = 0, npt = _getBiPtsNumber(); ipt < npt; ipt++)
    delete _bipts[ipt];
  _bipts.clear();
  _biPtsPerDirection = 0;
}

void Vario::_addBiTargetCheck(ABiTargetCheck* abpc)
{
  //_bipts.push_back(dynamic_cast<ABiTargetCheck*>(abpc->clone()));
  _bipts.push_back(abpc);
}

void Vario::_setListBiTargetCheck()
{
  _clearBiTargetCheck();
  for (int idir = 0, ndir = getDirectionNumber(); idir < ndir; idir++)
  {
    const DirParam dirparam = getDirParam(idir);
    _biPtsPerDirection = 0;

    // Add constraints linked to the Geometry (not performed for a Grid)
    if (! (_db->isGrid() && isDefinedForGrid() ))
    {
      BiTargetCheckGeometry *bipts = BiTargetCheckGeometry::create(_db->getNDim(),
                                                                   dirparam.getCodirs(),
                                                                   dirparam.getTolAngle(),
                                                                   dirparam.getBench(),
                                                                   dirparam.getCylRad(),
                                                                   getFlagAsym());
      _addBiTargetCheck(bipts);
      _biPtsPerDirection++;
    }

    // Add constraints linked to the Faults
    if (_varioparam.hasFaults())
    {
      BiTargetCheckFaults* bipts = BiTargetCheckFaults::create(_varioparam.getFaults());
      _addBiTargetCheck(bipts);
      _biPtsPerDirection++;
    }

    // Add constraints linked to the Date
    if (_varioparam.hasDate())
    {
      int idate = getIdate(idir);
      BiTargetCheckDate *bipts = BiTargetCheckDate::create(_varioparam.getDate(idate, 0),
                                                           _varioparam.getDate(idate, 1));
      _addBiTargetCheck(bipts);
      _biPtsPerDirection++;
    }

    // Add constraints linked to the Code
    if (dirparam.getOptionCode() != 0)
    {
      BiTargetCheckCode* bipts = BiTargetCheckCode::create(dirparam.getOptionCode(), dirparam.getTolCode());
      _addBiTargetCheck(bipts);
      _biPtsPerDirection++;
    }
  }
}

int Vario::_getBiPtsRank(int idir, int rank) const
{
  return (idir * _biPtsPerDirection + rank);
}

int Vario::prepare(const ECalcVario &calcul, bool defineList)
{
  if (_nVar <= 0)
  {
    messerr("The number of variables must be positive");
    return 1;
  }
  if (getDirectionNumber() <= 0)
  {
    messerr("The 'varioParam' argument must have some Direction defined");
    return 1;
  }

  // Preparation
  _calcul = calcul;
  _setFlagAsym();
  _setDPasFromGrid(isDefinedForGrid());
  if (internalVariableResize()) return 1;
  internalDirectionResize();

  // Define the list of ABiTargetCheckers corresponding to the calculation constraints
  if (defineList) _setListBiTargetCheck();

  return 0;
}

int Vario::compute(Db* db,
                   const ECalcVario &calcul,
                   bool flag_gen,
                   bool flag_sample,
                   bool verr_mode,
                   Model *model,
                   bool verbose)
{
  _db = db;
  _nVar = _db->getLocNumber(ELoc::Z);

  if (prepare(calcul)) return 1;

  if (_variogram_compute(_db, flag_gen, flag_sample, verr_mode, model, verbose))
  {
    messerr("Error when calculating the Variogram");
    return 1;
  }
  return 0;
}

int Vario::computeIndic(Db *db,
                        const ECalcVario& calcul,
                        bool flag_gen,
                        bool flag_sample,
                        bool verr_mode,
                        Model *model,
                        bool verbose,
                        int nfacmax)
{
  _db = db;
  int nvar = _db->getLocNumber(ELoc::Z);
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

  _nVar  = nclass;
  _means = props;
  _vars  = _varsFromProportions(props);

  if (prepare(calcul)) return 1;

  // Calculate the variogram of indicators
  if (_variogram_compute(_db, flag_gen, flag_sample, verr_mode, model, verbose))
  {
    messerr("Error when calculating the Variogram of Indicators");
    return 1;
  }

  // Delete the Indicators (created locally)
  _db->deleteColumnsByLocator(ELoc::Z);
  _db->setLocatorByUID(iatt, ELoc::Z);

  return 0;
}

/**
 * Reduce the current variogram by keeping a subset of variables and/or directions
 *
 * @param varcols Vector of variable ranks (starting from 0)
 * @param dircols Vector of direction ranks (starting from 0)
 * @param asSymmetric Turn the result into as Symmetrical function (i.e. variogram)
 */
void Vario::resetReduce(const VectorInt &varcols,
                   const VectorInt &dircols,
                   bool asSymmetric)
{
  VectorInt selvars;
  VectorInt seldirs;
  Vario vario_in(*this); // Copy the current variogram as input variogram
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
          setVar(vario_in.getVar(selvars[ivar], selvars[jvar]), ivar, jvar);
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
    for (int ipas = 0, npas = getLagNumber(idir); ipas < npas; ipas++)
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

/*****************************************************************************/
/*!
 **  Transform the experimental variogram from raw to gaussian space
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        Point Hermite anamorphosis
 **
 *****************************************************************************/
int Vario::transformZToY(const AAnam *anam)
{
  if (anam == nullptr)
  {
    messerr("The function 'transformZToY' needs an Anamorphosis");
    return 1;
  }
  AnamHermite* anamH = dynamic_cast<AnamHermite*>(anam->clone());
  if (anamH == nullptr)
  {
    messerr("The function 'transformZToY' needs a Hermite Anamorphosis");
    return 1;
  }
  if (getVariableNumber() != 1)
  {
    messerr("The function 'transformZToY' is restricted to Monovariate Variogram");
    return 1;
  }

  /* Loop on the directions of the variogram */

  double cvv = anam->getVariance();
  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    /* Loop on the lags */

    for (int ipas = 0, npas = getLagNumber(idir); ipas < npas; ipas++)
    {
      // TODO. GG must be a variogram of Zv -> Cv(h)
      setGgByIndex(idir,ipas,1. - anamH->invertVariance(cvv-getGgByIndex(idir, ipas)));
    }
  }

  // Modify the variance array
  setVar(1., 0,  0);

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
 **
 *****************************************************************************/
int Vario::transformYToZ(const AAnam *anam)
{
  CovCalcMode mode;

  /* Preliminary checks */

  if (anam == (AAnam*) NULL) return 1;
  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("The function 'transformYToZ' is restricted to Gaussian Anamorphosis");
    return 1;
  }
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam->clone());
  if (anam_hermite == nullptr)
  {
    messerr("The function 'transformYToZ' needs a Hermite Anamorphosis");
    return 1;
  }
  if (getVariableNumber() != 1)
  {
    messerr("The function 'transformYToZ' is restricted to Monovariate Variogram");
    return 1;
  }

  /* Loop on the directions of the variogram */

  double c0 = anam_hermite->computeVariance(1.);
  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    /* Loop on the lags */

    for (int ipas = 0, npas = getLagNumber(idir); ipas < npas; ipas++)
    {
      double chh = 1. - getGg(idir, 0, 0, ipas, false);
      double var = anam_hermite->computeVariance(chh);
      setGg(idir, 0, 0, ipas, c0 - var);
    }
  }

  setVar(c0, 0, 0);
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
 ** \param[in]  asCov     When true; produces a covariance
 **
 *****************************************************************************/
int Vario::modelRegularize(const Model& model,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles,
                           const CovCalcMode* mode,
                           bool asCov)
{
  int ndim = model.getDimensionNumber();
  int nvar = model.getVariableNumber();

  /* Preliminary checks */

  setNVar(nvar);
  internalVariableResize();
  internalDirectionResize();
  if (asCov)
    setCalcul(ECalcVario::COVARIANCE);
  else
    setCalcul(ECalcVario::VARIOGRAM);

  /* Initialize the variance array */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = model.evalCvv(ext, ndisc, angles, ivar, jvar, mode);
      setVar(value, ivar, jvar);
    }

  /* Loop on the directions */

  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {

    /* Loop on the number of lags */

    for (int ipas = 0, npas = getLagNumber(idir); ipas < npas; ipas++)
    {

      // Calculate the shift vector

      double dist = ipas * getDPas(idir);
      VectorDouble shift(ndim);
      for (int idim = 0; idim < ndim; idim++)
        shift[idim] = dist * getCodir(idir, idim);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          double value = model.evalCvvShift(ext, ndisc, shift, angles, ivar,
                                            jvar, mode);
          if (! asCov) value = getVar(ivar, jvar) - value;
          int iad = getDirAddress(idir, ivar, jvar, ipas, false, 0);
          setGgByIndex(idir, iad, value);
          setHhByIndex(idir, iad, dist);
          setSwByIndex(idir, iad, 1);
        }
    }
  }
  return 0;
}

void Vario::setDb(Db* db)
{
  _db = db;
  if (db != nullptr)
    _nVar = _db->getLocNumber(ELoc::Z);
}

int Vario::internalVariableResize() //TODO: to be called when nvar is modified...
{
  if (! _means.empty())
  {
    int nloc = static_cast<int>(_means.size());
    if (nloc != _nVar) _initMeans();
  }
  else
  {
    _initMeans();
  }

  if (! _vars.empty())
  {
    int nloc = static_cast<int>(_vars.size());
    if (nloc != _nVar * _nVar) _initVars();
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
        int rank = (!getFlagAsym()) ? i : i - getLagNumber(idir);
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
const ECalcVario Vario::getCalculType(const String& calcul_name)
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
  _means.resize(_nVar, 0.);
}

void Vario::setMeans(const VectorDouble& means)
{
  if (_means.empty()) _initMeans();
  if (! means.empty() && static_cast<int>(means.size()) == _nVar)
    _means = means;
}

void Vario::setMean(double mean, int ivar)
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

void Vario::setVar(double value, int ivar, int jvar)
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
                    bool asCov,
                    bool flagNorm) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ipas,true,0);
  if (IFFFF(iad)) return TEST;
  double val = _gg[idir][iad];
  double c0 = getVar(ivar, jvar);

  if (_flagAsym)
  {
    if (! asCov) val = c0 - val;
  }
  else
  {
    if (asCov) val = c0 - val;
  }

  if (flagNorm)
  {
    val /= c0;
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

  int npas = getLagNumber(idir);
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
 * @param compress When true, suppress lags where 'sw' <= 0
 * @return The vector of 'gg'
 */
VectorDouble Vario::getGgVec(int idir,
                             int ivar,
                             int jvar,
                             bool asCov,
                             bool flagNorm,
                             bool compress) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble gg;
  double c0 = 0.;
  if (asCov || flagNorm) c0 = getVar(ivar, jvar);
  int npas = getLagNumber(idir);

  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas-1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
      {
        double val = _gg[idir][iad];
        if (asCov && !getFlagAsym()) val = c0 - val;
        if (flagNorm) val /= c0;
        gg.push_back(val);
      }
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
    {
      double val = _gg[idir][iad];
      if (asCov && !getFlagAsym()) val = c0 - val;
      if (flagNorm) val /= c0;
      gg.push_back(val);
    }
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
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
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
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

void Vario::setGgVec(int idir, int ivar, int jvar, const VectorDouble& gg)
{
  if (!_isVariableValid(ivar))  return;
  if (!_isVariableValid(jvar))  return;
  if (!_isDirectionValid(idir)) return;
  int npas = getLagNumber(idir);
  if (npas != (int) gg.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      setGg(idir, ivar, jvar, ipas, gg[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setGg(idir, ivar, jvar, 0, gg[iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      setGg(idir, ivar, jvar, ipas, gg[iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      setGg(idir, ivar, jvar, ipas, gg[iad]);
    }
  }
}

/**
 * Returns the vector of distances for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @param compress When true, suppress lags where 'sw' <= 0
 * @return The vector of 'hh'
 */
VectorDouble Vario::getHhVec(int idir, int ivar, int jvar, bool compress) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble hh;
  int npas = getLagNumber(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        hh.push_back(_hh[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
      hh.push_back(_hh[idir][iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        hh.push_back(_hh[idir][iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        hh.push_back(_hh[idir][iad]);
    }
  }
  return hh;
}

void Vario::setHhVec(int idir, int ivar, int jvar, const VectorDouble& hh)
{
  if (!_isVariableValid(ivar))  return;
  if (!_isVariableValid(jvar))  return;
  if (!_isDirectionValid(idir)) return;

  int npas = getLagNumber(idir);
  if (npas != (int) hh.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      setHh(idir, ivar, jvar, ipas, hh[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setHh(idir, ivar, jvar, 0, hh[iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      setHh(idir, ivar, jvar, ipas, hh[iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      setHh(idir, ivar, jvar, ipas, hh[iad]);
    }
  }
}

/**
 * Returns the vector of weights for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @param compress When true, suppress lags where 'sw' <= 0
 * @return The vector of 'sw'
 */
VectorDouble Vario::getSwVec(int idir, int ivar, int jvar, bool compress) const
{
  if (!_isVariableValid(ivar))  return VectorDouble();
  if (!_isVariableValid(jvar))  return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();

  VectorDouble sw;
  int npas = getLagNumber(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        sw.push_back(_sw[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
      sw.push_back(_sw[idir][iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        sw.push_back(_sw[idir][iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        sw.push_back(_sw[idir][iad]);
    }
  }
  return sw;
}

void Vario::setSwVec(int idir, int ivar, int jvar, const VectorDouble& sw)
{
  if (!_isVariableValid(ivar))  return;
  if (!_isVariableValid(jvar))  return;
  if (!_isDirectionValid(idir)) return;
  int npas = getLagNumber(idir);
  if (npas != (int) sw.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
      setSw(idir, ivar, jvar, ipas, sw[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setSw(idir, ivar, jvar, 0, sw[iad]);
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
      setSw(idir, ivar, jvar, ipas, sw[iad]);
    }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      setSw(idir, ivar, jvar, ipas, sw[iad]);
    }
  }
}

/**
 * Returns the vector of utilization flags for a given pair of variables in a given direction
 * @param idir Direction
 * @param ivar First variable
 * @param jvar Second variable
 * @param compress When true, suppress lags where 'sw' <= 0
 * @return The vector of 'utilize'
 */
VectorDouble Vario::getUtilizeVec(int idir, int ivar, int jvar, bool compress) const
{
  if (!_isVariableValid(ivar)) return VectorDouble();
  if (!_isVariableValid(jvar)) return VectorDouble();
  if (!_isDirectionValid(idir)) return VectorDouble();
  VectorDouble utilize;
  int npas = getLagNumber(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ipas = npas - 1; ipas >= 0; ipas--)
     {
       iad = getDirAddress(idir, ivar, jvar, ipas, false, -1);
       if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
         utilize.push_back(_utilize[idir][iad]);
     }
     iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
     if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
       utilize.push_back(_utilize[idir][iad]);
     for (int ipas = 0; ipas < npas; ipas++)
     {
       iad = getDirAddress(idir, ivar, jvar, ipas, false, 1);
       if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        utilize.push_back(_utilize[idir][iad]);
     }
  }
  else
  {
    for (int ipas = 0; ipas < npas; ipas++)
    {
      iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
      if (! (IFFFF(iad) || (compress && _sw[idir][iad] <= 0.)))
        utilize.push_back(_utilize[idir][iad]);
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
  int i = getDirAddress(idir, ivar, jvar, 0, false, 0);
  return i;
}

int Vario::getNext(int ivar, int jvar, int idir, int shift) const
{
  if (!_isVariableValid(ivar))  return ITEST;
  if (!_isVariableValid(jvar))  return ITEST;
  if (!_isDirectionValid(idir)) return ITEST;

  VectorDouble sw;
  int npas = getLagNumber(idir);
  int count;
  if (_flagAsym) return ITEST;
  int iad = getDirSize(idir) - 1;
  count = 0;
  for (int ipas = 0; ipas < npas && count < shift; ipas++)
  {
    iad = getDirAddress(idir, ivar, jvar, ipas, true, 0);
    if (IFFFF(iad)) continue;
    if (_sw[idir][iad] != 0. && _hh[idir][iad] != 0.) count++;
  }
  return iad;
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
      int npas = getLagNumber(idir);
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
      ret = ret && _recordReadVec<double>(is, "Direction vector", codir, ndim);
    }
    if (! ret) return ret;

    SpaceRN space(ndim);
    DirParam dirparam = DirParam(npas, dpas, toldis, tolang, opt_code, 0,
                                 TEST, TEST, tolcode, VectorDouble(), codir, TEST,
                                 &space);
    if (isDefinedForGrid)
      dirparam.setGrincr(grincr);
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
      for (int idim = 0; idim < (int) dirparam.getNDim() && ret; idim++)
        ret = ret && _recordWrite<double>(os, "", dirparam.getCodir(idim));
      ret = ret && _commentWrite(os, "Direction coefficients");
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
    _nVar = db->getLocNumber(ELoc::Z);
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
  int npas = getLagNumber(idir);
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
  if (_db != nullptr && _db->isGrid() && flag_grid)
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

VectorDouble Vario::getGgs(int idir, int ivar, int jvar, const VectorInt& ipas) const
{
  VectorDouble values;
  if (ipas.empty()) return values;
  if (! _isDirectionValid(idir)) return values;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  for (int i = 0; i < (int) ipas.size(); i++)
  {
    if (ipas[i] >= 0 && ipas[i] < getDirSize(idir)) values.push_back(getGg(idir,ivar,jvar,ipas[i]));
  }
  return values;
}

VectorDouble Vario::setGgs(int idir, int ivar, int jvar, const VectorInt& ipas, const VectorDouble& values)
{
  if (ipas.empty()) return values;
  if (values.empty()) return values;
  if (! _isDirectionValid(idir)) return values;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  for (int i = 0; i < (int) ipas.size(); i++)
  {
    if (ipas[i] >= 0 && ipas[i] < getDirSize(idir) && i < (int) values.size())
      setGg(idir,ivar,jvar,ipas[i], values[i]);
  }
  return values;
}

VectorDouble Vario::getCodirs(int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return getDirParam(idir).getCodirs();
}

double Vario::getCodir(int idir, int idim) const
{
  if (! _isDirectionValid(idir)) return TEST;
  return getDirParam(idir).getCodir(idim);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  flag_gen     1 for calculation of generalized variogram
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 ** \param[in]  model        Model structure (triggers the KU option)
 ** \param[in]  verbose      Verbose flag
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
int Vario::_variogram_compute(Db *db,
                              int flag_gen,
                              int flag_sample,
                              int verr_mode,
                              Model *model,
                              int verbose)
{
  int error = 0;

  if (isDefinedForGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
    if (dbgrid == nullptr)
    {
      messerr("'Vario' is defined for Grid but 'db' is not organized as a grid");
      return 1;
    }
    if (flag_gen)
      error = _variogen_grid_calcul(dbgrid);
    else
      error = _variogrid_calcul(dbgrid);
  }
  else
  {
    if (flag_gen)
      error = _variogen_line_calcul(db);
    else
      error = _variogram_general(db, model, flag_sample, verr_mode, verbose);
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  model        Model structure (triggers the KU option)
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 ** \param[in]  verbose      Verbose flag
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
int Vario::_variogram_general(Db *db,
                              Model *model,
                              int flag_sample,
                              int verr_mode,
                              int verbose)
{
  int idir, error, flag_verr, flag_ku, nbfl;
  Vario_Order *vorder;
  VectorInt rindex;

  /* Initializations */

  error = 1;
  MODEL = model;
  flag_verr = flag_ku = nbfl = 0;
  if (db == nullptr) return (1);
  vorder = (Vario_Order*) NULL;
  manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    goto label_end;
  }
  if (_get_generalized_variogram_order() > 0)
  {
    messerr("Calculation does not allow generalized variogram definition");
    goto label_end;
  }

  /* Particular case of Transitive Covariogram */
  /* It is only coded in the by_sample case and uses the regression technique */

  if (getCalcul() == ECalcVario::COVARIOGRAM) flag_sample = 1;

  /* Auxiliary check for Variance Measurement Error */

  if (db->getLocNumber(ELoc::V) > 0 && verr_mode > 0)
  {
    vorder = vario_order_manage(1, 1, 0, vorder);
    flag_verr = 1;
  }

  // Auxiliary check for Drift removal. This is triggered only if the drift
  // contains at least one drift function different from Universality condition

  if (model != nullptr && model->isDriftDifferentDefined(VectorInt()))
  {
    if (vorder == (Vario_Order*) NULL)
      vorder = vario_order_manage(1, 1, 0, vorder);
    flag_ku = 1;
    manage_drift_removal(1, db, model);
  }

  /* Complementary checks */

  if (flag_verr && flag_ku)
  {
    messerr("These two options are incompatible");
    messerr("- Correction for the Variance of Error Measurements");
    messerr("- Correction for bias when removing the Drift");
    goto label_end;
  }
  if (flag_verr || flag_ku)
  {
    if (flag_sample)
    {
      messerr("The special Variogram option is incompatible with flag.sample");
      goto label_end;
    }
    if (!db->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Evaluate the drift coefficients */

  if (flag_ku)
  {
    if (estimate_drift_coefficients(db, verbose)) goto label_end;
  }

  /* Update the global statistics */

  _variogram_stats(db);

  /* Loop on the directions to evaluate */

  rindex = db->getSortArray();
  for (idir = 0; idir < getDirectionNumber(); idir++)
  {
    if (!flag_sample)
    {
      if (_variogram_calcul1(db, idir, rindex.data(), vorder)) goto label_end;
    }
    else
    {
      if (_variogram_calcul2(db, idir, rindex.data()))
        goto label_end;
    }

    if (vorder != (Vario_Order*) NULL)
      _variogram_calcul_internal(db, idir, vorder);
  }

  /* Posterior calculations when presence of Variance of Measurement errors */

  if (flag_verr)
  {
    for (idir = 0; idir < getDirectionNumber(); idir++)
    {
      if (_update_variogram_verr(db, idir, vorder, verr_mode)) goto label_end;
    }
  }

  /* Posterior update when filtering the bias attached to drift removal */

  if (flag_ku)
  {
    if (_update_variogram_ku(db, vorder, verbose)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: vorder = vario_order_manage(-1, 1, 0, vorder);
  manage_drift_removal(-1, db, model);
  return (error);
}

/****************************************************************************/
/*!
 **  Update the Variogram of Residuals when Drift has been removed
 **
 ** \return  Error returned code
 **
 ** \param[in]  db         Db description
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  verbose    Verbose flag
 **
 ** \remark The number of iterations used for the debiasing procedure
 ** \remark in presence of a drift can be defined using:
 ** \remark set_keypair("KU_Niter",newval)
 **
 *****************************************************************************/
int Vario::_update_variogram_ku(Db *db, Vario_Order *vorder, int verbose)
{
  Option_VarioFit optvar;
  Option_AutoFit mauto;
  double newval;
  int ifirst, ilast;
  Constraints constraints;

  /* Initializations */

  int ndim = MODEL->getDimensionNumber();
  int nbfl = MODEL->getDriftNumber();
  int niter_ku = (int) get_keypone("KU_Niter", 0);

  /* Core allocation */

  VectorDouble d1(ndim,0.);

  /* Loop on the iterations */

  for (int iter = 0; iter < niter_ku; iter++)
  {

    /* Perform the Automatic structure recognition */

    if (MODEL != nullptr)
    {
      if (model_auto_fit(this, MODEL, verbose, mauto, constraints, optvar))
        return 1;
    }

    /* Calculate the global bias correction terms */

    calculate_bias_global(db, d1);

    /* Optional printout */

    if (verbose)
    {
      message("Drift removal at iteration #d/%d\n", iter + 1, niter_ku);
      print_matrix("Drift Coefficients Matrix", 0, 1, nbfl, nbfl, NULL, get_DRFXGX().getValues().data());
    }

    /* Loop on the directions */

    for (int idir = 0; idir < getDirectionNumber(); idir++)
    {

      /* Loop on the lags */

      for (int ipas = 0; ipas < getLagNumber(idir); ipas++)
      {
        vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);
        if (ifirst > ilast) continue;

        /* Calculate the local bias correction terms */

        newval = _calculate_bias_local(db, vorder, ifirst, ilast);

        /* Patch the new value */

        setGg(idir, 0, 0, ipas, newval);
      }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the local bias terms
 **
 ** \param[in]  db        Db description
 ** \param[in]  vorder    Vario_Order structure
 ** \param[in]  ifirst    Rank of the first lag
 ** \param[in]  ilast     Rank of the last lag
 **
 *****************************************************************************/
double Vario::_calculate_bias_local(Db *db,
                                    Vario_Order *vorder,
                                    int ifirst,
                                    int ilast)
{
  int iech, jech, iiech, jjech;
  double dist, diff, tot0, tot1, tot2, totnum, result, v1, v2;

  /* Initializations */

  int nbfl = MODEL->getDriftNumber();

  /* Calculate the first corrected term */

  tot0 = tot1 = tot2 = totnum = 0.;
  for (int ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
    v1 = _get_IVAR(db, iech, 0);
    v2 = _get_IVAR(db, jech, 0);
    if (FFFF(v1) || FFFF(v2)) continue;

    iiech = _get_relative_sample_rank(db, iech);
    jjech = _get_relative_sample_rank(db, jech);

    diff = v1 - v2;
    tot0 += diff * diff;
    tot1 += get_bias_value(db, nbfl, iiech, jjech);
    tot2 += (get_DRFDIAG(iiech) + get_DRFDIAG(jjech)) / 2.;
    totnum += 1.;
  }

  tot0 /= 2.;
  result = tot0 - tot1 + tot2;
  if (totnum > 0.) result /= totnum;

  return (result);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram along lines
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 **
 *****************************************************************************/
int Vario::_variogen_line_calcul(Db *db)
{
  int idir, error, norder;

  /* Initializations */

  error = 1;
  if (db == nullptr) return (1);
  norder = _get_generalized_variogram_order();
  manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    return (1);
  }
  if (getVariableNumber() != 1)
  {
    messerr("The generalized variogram requires a single variable");
    return (1);
  }
  if (_get_generalized_variogram_order() == 0)
  {
    messerr("Calculation requires a generalized variogram definition");
    return (1);
  }
  if (! db->isGrid())
  {
    messerr("Calculation facility is dedicated to line architecture");
    return (1);
  }
  if (!db->hasLocVariable(ELoc::C))
  {
    messerr("Calculation facility requires the definition of a CODE");
    return (1);
  }

  /* Update the global statistics */

  _variogram_stats(db);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < getDirectionNumber(); idir++)
    _variogen_line(db, idir, norder);

  /* Set the error return code */

  error = 0;

  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the generalized variogram along lines
 **
 ** \param[in]  db      Db description
 ** \param[in]  idir    Rank of the Direction
 ** \param[in]  norder  Order of the generalized variogram
 **
 *****************************************************************************/
void Vario::_variogen_line(Db *db, int idir, int norder)
{
  int jech, keep, nvar;
  double value, zz;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  int nech = db->getSampleNumber();
  int npas = getLagNumber(idir);
  double dist0 = 0.;
  double dist = 0.;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);

  /* Loop on the first point */

  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    db->getSampleAsST(iech, T1);

    for (int ipas = 1; ipas < npas; ipas++)
    {
      value = _get_IVAR(db, iech, 0);
      if (FFFF(value)) break;
      dist0 = 0.;

      for (int iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
      {
        keep = 0;
        jech = iech + iwgt * ipas;
        if (jech < 0 || jech > nech) break;
        if (hasSel && !db->isActive(jech)) break;
        db->getSampleAsST(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! keepPair(idir, T1, T2, &dist)) continue;
        if (iwgt == 1) dist0 = dist;

        /* Evaluate the variogram */

        zz = _get_IVAR(db, jech, 0);
        if (FFFF(zz)) break;
        keep = 1;
        value += zz * VARWGT[norder][iwgt];
      }
      if (keep)
      {
        VARIO = this;
        value = value * value / NORWGT[norder];
        nvar = getVariableNumber();
        variogram_set(getCalcul(), nvar, ipas, 0, 0, 0, 1., dist0, value);
      }
    }
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  return;
}

/****************************************************************************/
/*!
 **  Patch the value of C(0) for covariances
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
void Vario::_variogram_patch_c00(Db *db, int idir)
{
  double z1, z2, s12w, s12wzz, ww, scale, value, m1, m2, sumw;

  /* Initializations */

  if (!getFlagAsym()) return;

  /* Calculate the C00 term */

  for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      int i = getDirAddress(idir, ivar, jvar, 0, false, 0);
      setHhByIndex(idir, i, 0.);

      scale = 1.;
      m1 = m2 = s12w = s12wzz = sumw = 0.;

      /* Calculate the statistics for each variable */

      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = _get_IVAR(db, iech, ivar);
        z2 = _get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;
        m1 += ww * z1;
        m2 += ww * z2;
        sumw += ww;
        value = z1 * z2;
        if (getCalcul() == ECalcVario::COVARIOGRAM)
        {
          scale = ww;
        }
        else
        {
          scale = ww * ww;
          s12w += scale;
        }
        s12wzz += scale * value;
        if (OptDbg::query(EDbg::VARIOGRAM))
          _print_debug(iech, iech, ivar, jvar, i, scale, value);
      }

      if (sumw > 0 && (getCalcul() == ECalcVario::COVARIANCE
          || getCalcul() == ECalcVario::COVARIANCE_NC))
      {
        m1 /= sumw;
        m2 /= sumw;
      }

      /* Final centering and normation */

      setSwByIndex(idir, i, sumw);
      if (getCalcul() == ECalcVario::COVARIOGRAM)
        setGgByIndex(idir, i, s12wzz);
      else if (getCalcul() == ECalcVario::COVARIANCE_NC)
        setGgByIndex(idir, i, s12wzz / s12w);
      else
        setGgByIndex(idir, i, s12wzz / s12w - m1 * m2);
    }
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram on grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 **
 *****************************************************************************/
int Vario::_variogen_grid_calcul(DbGrid *db)
{
  int idir, error, norder;

  /* Initializations */

  error = 1;
  if (db == nullptr) return (1);
  norder = _get_generalized_variogram_order();
  manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    return (1);
  }
  if (getVariableNumber() != 1)
  {
    messerr("The generalized variogram requires a single variable");
    return (1);
  }
  if (_get_generalized_variogram_order() == 0)
  {
    messerr("This calculation requires a generalized variogram definition");
    return (1);
  }
  if (! db->isGrid())
  {
    messerr("This calculation facility is dedicated to grid architecture");
    return (1);
  }

  /* Update the global statistics */

  _variogram_stats(db);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < getDirectionNumber(); idir++)
  {
    error = _variogen_grid(db, idir, norder);
    if (error) break;
  }

  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram on grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 **
 *****************************************************************************/
int Vario::_variogrid_calcul(DbGrid *db)
{
  int idir, error, iadd_new, iatt_old, iech;
  double maille;

  /* Initializations */

  error = 1;
  iadd_new = iatt_old = -1;
  if (db == nullptr) return (1);
  manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    goto label_end;
  }
  if (! db->isGrid())
  {
    messerr("This calculation facility is dedicated to grid architecture");
    goto label_end;
  }
  if (_get_generalized_variogram_order() > 0)
  {
    messerr("This calculation does not allow generalized variogram definition");
    goto label_end;
  }

  /* In the case of Covariogram, add the weight set to the scale */

  if (getCalcul() == ECalcVario::COVARIOGRAM)
  {
    iatt_old = db_attribute_identify(db, ELoc::W, 0);
    iadd_new = db->addColumnsByConstant(1, 0.);
    if (iadd_new < 0) goto label_end;
    db->setLocatorByUID(iadd_new, ELoc::W);
    maille = db_grid_maille(db);
    for (iech = 0; iech < db->getSampleNumber(); iech++)
      db->setLocVariable(ELoc::W, iech, 0, maille);
  }

  /* Update the global statistics */

  _variogram_stats(db);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < getDirectionNumber(); idir++)
  {
    error = _variogram_grid(db, idir);
    if (error) break;
  }

  /* Delete the additional weight variable (optional) */

  if (getCalcul() == ECalcVario::COVARIOGRAM)
  {
    if (iadd_new > 0) db->deleteColumnByUID(iadd_new);
    if (iatt_old > 0) db->setLocatorByUID(iatt_old, ELoc::W);
  }

  /* Set the error return code */

  error = 0;

  label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Return the Order of the generalized variogram
 **
 ** \return Order of the Generalized covariance
 **
 *****************************************************************************/
int Vario::_get_generalized_variogram_order()
{
  int norder = 0;
  if (getCalcul() == ECalcVario::GENERAL1) norder = 1;
  if (getCalcul() == ECalcVario::GENERAL2) norder = 2;
  if (getCalcul() == ECalcVario::GENERAL3) norder = 3;
  return (norder);
}

/****************************************************************************/
/*!
 **  Calculates the statistics for the variogram calculations
 **
 ** \param[in]  db      Db descriptor
 **
 *****************************************************************************/
void Vario::_variogram_stats(Db *db)
{
  double z1, z2, ww;

  /* Initializations */

  for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
  {
    setMean(0., ivar);
    for (int jvar = 0; jvar < db->getLocNumber(ELoc::Z); jvar++)
      setVar(0., ivar, jvar);
  }

  /* Loop on the variables */

  for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
  {
    double s1w = 0.;
    double s1z = 0.;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      ww = db->getWeight(iech);
      if (FFFF(ww) || ww < 0.) continue;
      z1 = _get_IVAR(db, iech, ivar);
      s1w += ww;
      s1z += z1;
    }

    if (s1w <= 0.) continue;
    setMean(s1z / s1w, ivar);
  }

  for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      double s12w = 0.;
      double s12wz1 = 0.;
      double s12wz2 = 0.;
      double s12wzz = 0.;

      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;

        z1 = _get_IVAR(db, iech, ivar);
        z2 = _get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;

        s12w += ww;
        s12wz1 += ww * z1;
        s12wz2 += ww * z2;
        s12wzz += ww * z1 * z2;
      }
      if (s12w <= 0.) continue;

      if (getCalcul() == ECalcVario::COVARIOGRAM)
      {
        setVar(s12wzz, ivar, jvar);
        setVar(s12wzz, jvar, ivar);
      }
      else
      {
        setVar(s12wzz / s12w - (s12wz1 / s12w) * (s12wz2 / s12w), ivar, jvar);
        setVar(s12wzz / s12w - (s12wz1 / s12w) * (s12wz2 / s12w), jvar, ivar);
      }
    }

  // Modification when the ultimate variogram is a transformed one

  if (getCalcul() == ECalcVario::TRANS1)
  {
    for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -getVar(ivar, jvar) / getVar(jvar, jvar);
        setVar(value, ivar, jvar);
        setVar(value, jvar, ivar);
      }
  }
  else if (getCalcul() == ECalcVario::TRANS2)
  {
    for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -getVar(ivar, jvar) / getVar(ivar, ivar);
        setVar(value, ivar, jvar);
        setVar(value, jvar, ivar);
      }
  }
  else if (getCalcul() == ECalcVario::BINORMAL)
  {
    for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
      for (int jvar = 0; jvar < db->getLocNumber(ELoc::Z); jvar++)
        if (ivar != jvar)
          setVar(getVar(ivar, jvar) / sqrt(getVar(ivar, ivar) * getVar(jvar, jvar)), ivar, jvar);
  }
  return;
}

/****************************************************************************/
/*!
 **  Update the Variogram when Variance of measurement error is available
 **
 ** \return  Error returned code
 **
 ** \param[in]  db        Db description
 ** \param[in]  idir      Rank of the direction
 ** \param[in]  vorder    Vario_Order structure
 ** \param[in]  verr_mode Mode of variogram correction (1, 2 or 3)
 **
 *****************************************************************************/
int Vario::_update_variogram_verr(Db *db,
                                  int idir,
                                  Vario_Order *vorder,
                                  int verr_mode)
{
  int ifirst, ilast, iech, jech, number, nfois;
  double dist, value, g_old, diff, sumt, sumb, wgt, sval, gval;
  static double tol = EPSILON5;
  static int maxiter = 100;

  /* Initializations */

  int npas = getLagNumber(idir);

  /* Loop on the lags */

  for (int ipas = 0; ipas < npas; ipas++)
  {
    vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);
    if (ifirst > ilast) continue;

    /* Dispatch according to the method */

    switch (verr_mode)
    {

      case 1: /* Simple bias correction */
        number = 0;
        value = 0.;
        for (int ipair = ifirst; ipair < ilast; ipair++)
        {
          vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
          value += _s(db, iech, jech);
          number++;
        }
        value = (number > 0) ? value / number : 0.;
        setGg(idir, 0, 0, ipas, MAX(0, getGg(idir, 0, 0, ipas) - value));
        break;

      case 2:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          nfois++;
          g_old = getGg(idir, 0, 0, ipas);
          sumt = sumb = 0.;
          for (int ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = _s(db, iech, jech);
            gval = _g(db, iech, jech);
            value = sval + getGg(idir, 0, 0, ipas);
            wgt = 1. / (value * value);
            sumt += wgt * (gval - sval);
            sumb += wgt;
          }
          setGg(idir, 0, 0, ipas, sumt / sumb);
          diff = ABS(getGg(idir, 0, 0, ipas) - g_old);
        }
        setGg(idir, 0, 0, ipas, MAX(0, getGg(idir, 0, 0, ipas)));
        if (nfois == maxiter && OptDbg::query(EDbg::CONVERGE))
          message("Convergence not reached for lag %d\n", ipas + 1);
        break;

      case 3:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          g_old = getGg(idir, 0, 0, ipas);
          sumt = sumb = 0.;
          for (int ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = _s(db, iech, jech);
            gval = _g(db, iech, jech);
            value = sval + getGgByIndex(idir, ipas);
            wgt = getGg(idir, 0, 0, ipas) / value;
            sumt += wgt * gval;
            sumb += 1.;
          }
          setGg(idir, 0, 0, ipas, sumt / sumb);
          diff = ABS(getGg(idir, 0, 0, ipas) - g_old);
        }
        if (nfois == maxiter && OptDbg::query(EDbg::CONVERGE))
          message("Convergence not reached for lag %d\n", ipas + 1);
        break;

      default:
        messerr("The method (%d) for updating the Variogram", verr_mode);
        messerr("calculated with Variance of Measurement Error is unknown");
        return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Local function defined as the half-sum of the Variance of Measurement
 **  Error on the two points constituting a pair
 **
 ** \return  Returned value
 **
 ** \param[in]  db     Db description
 ** \param[in]  iech   Rank of the first sample
 ** \param[in]  jech   Rank of the second sample
 **
 *****************************************************************************/
double Vario::_s(Db *db, int iech, int jech)
{
  double value = 0.5 * (db->getLocVariable(ELoc::V, iech, 0) +
                        db->getLocVariable(ELoc::V, jech, 0));
  return (value);
}

/****************************************************************************/
/*!
 **  Local function defined as the half-sum squared difference of values
 **  between the two points constituting a pair
 **
 ** \return  Returned value
 **
 ** \param[in]  db     Db description
 ** \param[in]  iech   Rank of the first sample
 ** \param[in]  jech   Rank of the second sample
 **
 *****************************************************************************/
double Vario::_g(Db *db, int iech, int jech)
{
  double value = _get_IVAR(db, iech, 0) - _get_IVAR(db, jech, 0);
  value = value * value / 2.;
  return (value);
}

/****************************************************************************/
/*!
 **  Return the relative sample rank corresponding to an abolsute number
 **
 ** \return The corresponding relative sample rank
 **
 ** \param[in]  db     Db structure
 ** \param[in]  iech0  Absolute rank of the first sample
 **
 *****************************************************************************/
int Vario::_get_relative_sample_rank(Db *db, int iech0)
{
  int iiech = 0;
  for (int iech = 0, nech = db->getSampleNumber(); iech < nech; iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    if (iech == iech0) return (iiech);
    iiech++;
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Perform the Variogram evaluation when the Geometry has been established
 **
 ** \param[in]  db     Db description
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  vorder Vario_Order structure
 **
 *****************************************************************************/
void Vario::_variogram_calcul_internal(Db *db, int idir, Vario_Order *vorder)
{
  int iech, jech, ifirst, ilast;
  double dist;

  /* Initializations */

  int npas = getLagNumber(idir);

  /* Loop on the lags */

  for (int ipas = 0; ipas < npas; ipas++)
  {
    vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);

    /* Loop on the pairs contributing to this lag */

    for (int ipair = ifirst; ipair < ilast; ipair++)
    {
      vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);

      /* Evaluate the variogram */

      VARIO = this;
      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      variogram_evaluate(db, getCalcul(), getVariableNumber(), iech, jech, ipas,
                         dist, 1, variogram_set);
    }
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  return;
}

/****************************************************************************/
/*!
 **  Evaluate the variogram using the traditional method
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  rindex Array of sorted samples
 ** \param[in]  vorder Vario_Order structure
 **
 *****************************************************************************/
int Vario::_variogram_calcul1(Db *db,
                              int idir,
                              int *rindex,
                              Vario_Order *vorder)
{
  int iech, jech, ipas, npair, ideb;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  DirParam dirparam = getDirParam(idir);
  int nech = db->getSampleNumber();
  double maxdist = getMaximumDistance(idir);
  const VarioParam& varioparam = getVarioParam();

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam.isDateUsed(db);
  double dist = 0.;

  /* Loop on the first point */

  for (int iiech = 0; iiech < nech - 1; iiech++)
  {
    iech = rindex[iiech];
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsST(iech, T1);

    ideb = (hasDate) ? 0 : iiech + 1;
    for (int jjech = ideb; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (db->getDistance1D(iech, jech) > maxdist) break;
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      ipas = dirparam.getLagRank(dist);
      if (IFFFF(ipas)) continue;

      /* Case of internal storage */

      if (vorder != (Vario_Order*) NULL)
      {
        vario_order_add(vorder, iech, jech, NULL, NULL, ipas, idir, dist);
      }
      else
      {

        /* Evaluate the variogram */

        VARIO = this;
        IDIRLOC = idir;
        IECH1 = iech;
        IECH2 = jech;
        variogram_evaluate(db, getCalcul(), getVariableNumber(), iech, jech,
                           ipas, dist, 1, variogram_set);
      }
    }
  }

  /* Internal storage */

  if (vorder != (Vario_Order*) NULL)
  {
    vorder = vario_order_final(vorder, &npair);
  }
  else
  {

    /* Scale the variogram calculations */

    _variogram_scale(idir);

    /* Center the covariance function */

    _covariance_center(db, idir);

    /* Patch the central value */

    _variogram_patch_c00(db, idir);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the variogram by sample
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  idir   Rank of the direction
 ** \param[in]  rindex Array of sorted samples
 **
 *****************************************************************************/
int Vario::_variogram_calcul2(Db *db, int idir, int *rindex)
{
  int iech, jech, i, ipas, ideb;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Initializations */

  const VarioParam &varioparam = getVarioParam();
  const DirParam &dirparam = getDirParam(idir);
  int nech = db->getSampleNumber();
  int size = getDirSize(idir);
  double maxdist = getMaximumDistance(idir);

  /* Core allocation */

  VectorDouble gg_sum(size, 0);
  VectorDouble hh_sum(size, 0);
  VectorDouble sw_sum(size, 0);

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam.isDateUsed(db);
  double w1 = 1.;
  double dist = 0.;

  /* Loop on the first sample */

  for (int iiech = 0; iiech < nech; iiech++)
  {
    iech = rindex[iiech];
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight)
    {
      w1 = db->getWeight(iech);
      if (FFFF(w1)) continue;
    }
    db->getSampleAsST(iech, T1);

    /* Looking for the second sample */

    ideb = (hasDate) ? 0 : iiech + 1;
    for (int jjech = ideb; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (db->getDistance1D(iech, jech) > maxdist) break;
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      ipas = dirparam.getLagRank(dist);
      if (IFFFF(ipas)) continue;

      /* Evaluate the variogram */

      VARIO = this;
      IECH1 = iech;
      IECH2 = jech;
      variogram_evaluate(db, getCalcul(), getVariableNumber(), iech, jech, ipas,
                         dist, 1, variogram_set);
    }

    /* Cumulate to the global variogram */

    for (i = 0; i < size; i++)
    {
      if (getSwByIndex(idir, i) <= 0) continue;
      sw_sum[i] += w1;
      gg_sum[i] += w1 * getGgByIndex(idir, i) / getSwByIndex(idir, i);
      hh_sum[i] += w1 * getHhByIndex(idir, i) / getSwByIndex(idir, i);
    }
  }

  /* Copy the cumulated variogram into the Vario structure */

  for (i = 0; i < size; i++)
  {
    setGgByIndex(idir, i, gg_sum[i]);
    setHhByIndex(idir, i, hh_sum[i]);
    setSwByIndex(idir, i, sw_sum[i]);
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the variogram on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db    Db description
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
int Vario::_variogram_grid(DbGrid *db, int idir)
{
  int *indg1, *indg2, jech;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Initializations */

  int error = 1;
  int nech = db->getSampleNumber();
  int npas = getLagNumber(idir);
  const DirParam &dirparam = getDirParam(idir);
  indg1 = indg2 = nullptr;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double dist = 0.;

  /* Core allocation */

  indg1 = db_indg_alloc(db);
  if (indg1 == nullptr) goto label_end;
  indg2 = db_indg_alloc(db);
  if (indg2 == nullptr) goto label_end;

  /* Loop on the first point */

  for (int iech = 0; iech < nech; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsST(iech, T1);
    db_index_sample_to_grid(db, iech, indg1);

    for (int ipas = 1; ipas < npas; ipas++)
    {
      for (int idim = 0; idim < db->getNDim(); idim++)
        indg2[idim] = indg1[idim] + (int) (ipas * getGrincr(idir, idim));
      jech = db_index_grid_to_sample(db, indg2);
      if (jech < 0) continue;

      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Evaluate the variogram */

      VARIO = this;
      dist = ipas * dirparam.getDPas();
      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      variogram_evaluate(db, getCalcul(), getVariableNumber(), iech, jech, ipas,
                         dist, 1, variogram_set);
    }
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  /* Set the error return status */

  error = 0;

  label_end:

  /* Core deallocation */

  indg1 = db_indg_free(indg1);
  indg2 = db_indg_free(indg2);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the generalized variogram on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db description
 ** \param[in]  idir    Rank of the direction
 ** \param[in]  norder  Order of the generalized variogram
 **
 *****************************************************************************/
int Vario::_variogen_grid(DbGrid *db, int idir, int norder)
{
  int *indg1, *indg2;
  int jech, nech, idim, error, npas, keep, nvar;
  double zz, value;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  indg1 = indg2 = nullptr;
  npas = getLagNumber(idir);
  const DirParam &dirparam = getDirParam(idir);

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  double dist = 0.;

  /* Core allocation */

  indg1 = db_indg_alloc(db);
  if (indg1 == nullptr) goto label_end;
  indg2 = db_indg_alloc(db);
  if (indg2 == nullptr) goto label_end;

  /* Loop on the first point */

  for (int iech = 0; iech < nech; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    db->getSampleAsST(iech, T1);
    db_index_sample_to_grid(db, iech, indg1);

    for (int ipas = 1; ipas < npas; ipas++)
    {
      value = _get_IVAR(db, iech, 0);
      if (FFFF(value)) break;
      dist = ipas * getDPas(idir);

      for (int iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
      {
        keep = 0;
        for (idim = 0; idim < db->getNDim(); idim++)
          indg2[idim] = indg1[idim] + (int) (ipas * iwgt * dirparam.getGrincr(idim));

        jech = db_index_grid_to_sample(db, indg2);
        if (jech < 0) continue;
        if (hasSel && !db->isActive(jech)) continue;
        db->getSampleAsST(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! keepPair(idir, T1, T2, &dist)) continue;

        /* Evaluate the variogram */

        zz = _get_IVAR(db, jech, 0);
        if (FFFF(zz)) break;
        keep = 1;
        value += zz * VARWGT[norder][iwgt];
      }
      if (keep)
      {
        VARIO = this;
        value = value * value / NORWGT[norder];
        nvar = getVariableNumber();
        variogram_set(getCalcul(), nvar, ipas, 0, 0, 0, 1., dist, value);
      }
    }
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  /* Set the error return status */

  error = 0;

  label_end:

  /* Core deallocation */

  indg1 = db_indg_free(indg1);
  indg2 = db_indg_free(indg2);
  return (error);
}

/****************************************************************************/
/*!
 **  Internal function for setting a variogram value
 **
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  ipas        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  orient      Orientation
 ** \param[in]  ww          Weight
 ** \param[in]  dist        Distance
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
void variogram_set(const ECalcVario &calcul_type,
                   int /*nvar*/,
                   int ipas,
                   int ivar,
                   int jvar,
                   int orient,
                   double ww,
                   double dist,
                   double value)
{
  int i = VARIO->getDirAddress(IDIRLOC, ivar, jvar, ipas, false, orient);
  VARIO->updateGgByIndex(IDIRLOC, i, ww * value);
  if (calcul_type == ECalcVario::POISSON)
    VARIO->updateGgByIndex(IDIRLOC, i, -VARIO->getMean(ivar) / 2.);
  VARIO->updateHhByIndex(IDIRLOC, i, ww * dist);
  VARIO->updateSwByIndex(IDIRLOC, i, ww);
  return;
}

/****************************************************************************/
/*!
 **  Printout function for Debug case
 **
 ** \param[in]  iech1 Rank of the first sample
 ** \param[in]  iech2 Rank of the second sample
 ** \param[in]  ivar  Rank of the first variable
 ** \param[in]  jvar  Rank of the second variable
 ** \param[in]  ilag  Rank of the Lag
 ** \param[in]  scale Weighting factor
 ** \param[in]  value Variogram value
 **
 *****************************************************************************/
void Vario::_print_debug(int iech1,
                         int iech2,
                         int ivar,
                         int jvar,
                         int ilag,
                         double scale,
                         double value)
{
  message("Samples: %d/%d - Variables: %d/%d - Weight: %lf - Lag: %d - Variogram: %lf\n",
          iech1 + 1, iech2 + 1, ivar + 1, jvar + 1, scale, ilag, value);
}

/****************************************************************************/
/*!
 **  Center the covariance calculations
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  idir  Rank of the direction
 **
 *****************************************************************************/
void Vario::_covariance_center(Db *db, int idir)
{
  double m1, m2, sumw, z1, z2, ww;
  if (!getFlagAsym()) return;

  /* Scale the experimental variogram quantities */

  for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      /* Calculate the mean for each variable */

      m1 = m2 = sumw = 0.;
      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = _get_IVAR(db, iech, ivar);
        z2 = _get_IVAR(db, iech, jvar);
        if (FFFF(z1) || FFFF(z2)) continue;
        m1 += ww * z1;
        m2 += ww * z2;
        sumw += ww;
      }

      if (sumw > 0 && (getCalcul() == ECalcVario::COVARIANCE
          || getCalcul() == ECalcVario::COVARIANCE_NC))
      {
        m1 /= sumw;
        m2 /= sumw;
      }

      /* Perform the Centering */

      if (!(getCalcul() == ECalcVario::COVARIOGRAM || getCalcul() == ECalcVario::COVARIANCE_NC))
        for (int i = 0; i < getLagTotalNumber(idir); i++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          if (getSwByIndex(idir, j) > 0)
            setGgByIndex(idir, j, getGgByIndex(idir, j) - m1 * m2);
        }
    }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the geometry for a given direction
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 **
 ** \param[out] vorder Vario_Order structure
 ** \param[out] npair  Number of pairs
 **
 *****************************************************************************/
int Vario::geometryCompute(Db *db, Vario_Order *vorder, int *npair)
{
  double maxdist;
  int iiech, iech, jjech, jech, nech, ipas, idir, ideb;
  VectorInt rindex;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Initializations */

  if (db == nullptr) return (1);
  const VarioParam &varioparam = getVarioParam();
  double dist = 0.;

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    return 1;
  }
  if (_get_generalized_variogram_order() > 0)
  {
    messerr("This calculation does not allow generalized variogram definition");
    return 1;
  }

  /* Sort the data */
  rindex = db->getSortArray();

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam.isDateUsed(db);
  nech = db->getSampleNumber();

  /* Loop on the directions */

  for (idir = 0; idir < getDirectionNumber(); idir++)
  {
    const DirParam &dirparam = getDirParam(idir);
    maxdist = getMaximumDistance(idir);

    /* Loop on the first point */

    for (iiech = 0; iiech < nech - 1; iiech++)
    {
      iech = rindex[iiech];
      if (hasSel && !db->isActive(iech)) continue;
      if (hasWeight && FFFF(db->getWeight(iech))) continue;
      db->getSampleAsST(iech, T1);

      ideb = (hasDate) ? 0 : iiech + 1;
      for (jjech = ideb; jjech < nech; jjech++)
      {
        jech = rindex[jjech];
        if (db->getDistance1D(iech, jech) > maxdist) break;
        if (hasSel && !db->isActive(jech)) continue;
        if (hasWeight && FFFF(db->getWeight(jech))) continue;
        db->getSampleAsST(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! keepPair(idir, T1, T2, &dist)) continue;

        /* Get the rank of the lag */

        ipas = dirparam.getLagRank(dist);
        if (IFFFF(ipas)) continue;

        /* Case of internal storage */

        vario_order_add(vorder, iech, jech, NULL, NULL, ipas, idir, dist);
      }
    }
  }

  /* Sort the geometry */

  vorder = vario_order_final(vorder, npair);

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the variogram extension for a pair of variables
 **
 ** \param[in]  ivar      Rank of the first variable
 ** \param[in]  jvar      Rank of the second variable
 ** \param[in]  idir0     Rank of the direction (-1 for all)
 ** \param[in]  flag_norm 1 if the variogram must be normalized by variance
 ** \param[in]  flag_vars 1 if the global statistics must be taken into account
 ** \param[in]  distmin   Minimum along the distance axis
 ** \param[in]  distmax   Maximum along the distance axis
 ** \param[in]  varmin    Minimum along the variogram (or covariance) axis
 ** \param[in]  varmax    Maximum along the variogram (or covariance) axis
 **
 ** \param[out]  flag_hneg 1 if the distance scale can be negative
 ** \param[out]  flag_gneg 1 if the variogram scale can be negative
 ** \param[out]  c0        Value of the variogram at the origin
 ** \param[out]  hmin      Minimum distance
 ** \param[out]  hmax      Maximum distance
 ** \param[out]  gmin      Minimum variogram value
 ** \param[out]  gmax      Maximum variogram value
 **
 *****************************************************************************/
void Vario::getExtension(int ivar,
                         int jvar,
                         int idir0,
                         int flag_norm,
                         int flag_vars,
                         double distmin,
                         double distmax,
                         double varmin,
                         double varmax,
                         int *flag_hneg,
                         int *flag_gneg,
                         double *c0,
                         double *hmin,
                         double *hmax,
                         double *gmin,
                         double *gmax)
{
  double hh, gg;
  int i, j, idir, jdir, ndir;
  double tol = 0.1;

  /* Initializations */

  (*hmin) = 0.;
  (*gmin) = 0.;
  (*hmax) = -1.e30;
  (*gmax) = -1.e30;
  if (getFlagAsym())
  {
    *c0 = getGgByIndex(0, getDirAddress(0, ivar, jvar, 0, false, 0));
  }
  else
  {
    *c0 = getVar(ivar, jvar);
  }
  if (_get_generalized_variogram_order() > 0) (*c0) = TEST;
  if (FFFF(*c0) && flag_norm)
  {
    messerr("The Normalization option is discarded for this variogram");
    messerr("probably as it corresponds to a generalized variogram");
    flag_norm = 0;
  }
  (*flag_hneg) = (ivar != jvar && getFlagAsym());
  (*flag_gneg) = (ivar != jvar || getFlagAsym());

  /* Loop on the directions */

  ndir = getDirectionNumber();
  if (idir0 >= 0) ndir = 1;
  for (jdir = 0; jdir < ndir; jdir++)
  {
    idir = (idir0 >= 0) ? idir0 : jdir;
    for (i = 0; i < getLagTotalNumber(idir); i++)
    {
      j = getDirAddress(idir, ivar, jvar, i, true, 0);
      if (getSwByIndex(idir, j) <= 0) continue;
      hh = getHhByIndex(idir, j);
      gg = getGgByIndex(idir, j);
      if (FFFF(hh) || FFFF(gg)) continue;
      if (flag_norm) gg /= (*c0);
      if (!FFFF(distmin) && hh < distmin) continue;
      if (!FFFF(distmax) && hh > distmax) continue;
      if (hh < (*hmin)) (*hmin) = hh;
      if (hh > (*hmax)) (*hmax) = hh;
      if (gg < (*gmin)) (*gmin) = gg;
      if (gg > (*gmax)) (*gmax) = gg;
    }
  }

  if (flag_norm) (*c0) = 1.;
  if (!FFFF(*c0) && flag_vars)
  {
    if ((*c0) < (*gmin)) (*gmin) = (*c0);
    if ((*c0) > (*gmax)) (*gmax) = (*c0);
  }

  /* Expand the vertical graphic scales */

  (*gmax) *= (1. + tol);
  if ((*gmin) < 0) (*gmin) *= (1. + tol);

  /* Correction due to variogram type */

  if (*flag_hneg)
  {
    (*hmax) = MAX(ABS(*hmin), ABS(*hmax));
    (*hmin) = -(*hmax);
  }

  if (*flag_gneg)
  {
    (*gmax) = MAX(ABS(*gmin), ABS(*gmax));
    if (ivar != jvar) (*gmin) = -(*gmax);
  }
  else
  {
    (*gmin) = 0.;
  }

  /* Truncation to limits provided by the user */

  if (!FFFF(distmax)) (*hmax) = distmax;
  if (!FFFF(distmin)) (*hmin) = distmin;
  if (!FFFF(varmax)) (*gmax) = varmax;
  if (!FFFF(varmin)) (*gmin) = varmin;

  return;
}

/****************************************************************************/
/*!
 **  Evaluate the covariance for directional variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  idir   Rank of the Direction
 ** \param[in]  ncomp  Number of components
 ** \param[in]  rindex Array of sorted samples
 **
 *****************************************************************************/
int Vario::_variovect_calcul(Db *db,
                             int idir,
                             int ncomp,
                             int *rindex)
{
  int iech, jech, ipas, i, icomp;
  double w1, w2, zi1, zi2, zj1, zj2, v12, v21, di1, di2, dj1, dj2;
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  const DirParam &dirparam = getDirParam(idir);
  int nech = db->getSampleNumber();
  int nvar = getVariableNumber();
  double maxdist = getMaximumDistance(idir);

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double dist = 0.;

  /* Loop on the first point */

  for (int iiech = 0; iiech < nech - 1; iiech++)
  {
    iech = rindex[iiech];
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsST(iech, T1);

    for (int jjech = iiech + 1; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (db->getDistance1D(iech, jech) > maxdist) break;
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      ipas = dirparam.getLagRank(dist);
      if (IFFFF(ipas)) continue;

      w1 = db->getWeight(iech);
      w2 = db->getWeight(jech);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {

          /* Evaluate the variogram */

          v12 = v21 = di1 = di2 = dj1 = dj2 = 0.;
          for (icomp = 0; icomp < ncomp; icomp++)
          {
            zi1 = _get_IVAR(db, iech, ivar * ncomp + icomp);
            zi2 = _get_IVAR(db, iech, jvar * ncomp + icomp);
            zj1 = _get_IVAR(db, jech, ivar * ncomp + icomp);
            zj2 = _get_IVAR(db, jech, jvar * ncomp + icomp);
            if (FFFF(zi1) || FFFF(zi2) || FFFF(zj1) || FFFF(zj2))
            {
              v12 = v21 = TEST;
              break;
            }
            else
            {
              v12 += zi1 * zj2;
              v21 += zi2 * zj1;
              di1 += zi1 * zi1;
              di2 += zi2 * zi2;
              dj1 += zj1 * zj1;
              dj2 += zj2 * zj2;
            }
          }
          if (FFFF(v12) || FFFF(v21)) continue;
          if (ABS(di1) < EPSILON8 || ABS(di2) < EPSILON8) continue;
          if (ABS(dj1) < EPSILON8 || ABS(dj2) < EPSILON8) continue;
          di1 = sqrt(di1);
          di2 = sqrt(di2);
          dj1 = sqrt(dj1);
          dj2 = sqrt(dj2);
          v12 = ABS(v12) / (di1 * dj2);
          v21 = ABS(v21) / (di2 * dj1);

          i = getDirAddress(idir, ivar, jvar, ipas, false, 1);
          setGgByIndex(idir, i, getGgByIndex(idir, i) + w1 * w2 * v12);
          setHhByIndex(idir, i, getHhByIndex(idir, i) + w1 * w2 * dist);
          setSwByIndex(idir, i, getSwByIndex(idir, i) + w1 * w2);

          i = getDirAddress(idir, ivar, jvar, ipas, false, -1);
          setGgByIndex(idir, i, getGgByIndex(idir, i) + w1 * w2 * v21);
          setHhByIndex(idir, i, getHhByIndex(idir, i) + w1 * w2 * dist);
          setSwByIndex(idir, i, getSwByIndex(idir, i) + w1 * w2);
        }
    }
  }

  /* Scale the variogram calculations */

  _variogram_scale(idir);

  /* Center the covariance function */

  _covariance_center(db, idir);

  /* Patch the central value */

  _variogram_patch_c00(db, idir);

  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental covariance for directional variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  ncomp  Number of components
 **
 *****************************************************************************/
int Vario::variovectCompute(Db *db, int ncomp)
{
  if (db == nullptr) return (1);
  manage_drift_removal(0, NULL, NULL);

  /* Preliminary checks */

  if (db->getNDim() != getDimensionNumber() || db->getLocNumber(ELoc::Z) != getVariableNumber() * ncomp)
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getLocNumber(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getDimensionNumber(),
            getVariableNumber());
    messerr("Number of components = %d", ncomp);
    return (1);
  }

  /* Update the global statistics */

  _variovect_stats(db, ncomp);

  /* Loop on the directions to evaluate */

  VectorInt rindex = db->getSortArray();
  for (int idir = 0; idir < getDirectionNumber(); idir++)
  {
    if (_variovect_calcul(db, idir, ncomp, rindex.data())) return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculates the statistics for the covariance calculations
 **  for directional variables
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  ncomp  Number of components
 **
 *****************************************************************************/
void Vario::_variovect_stats(Db *db, int ncomp)
{
  double vi, vj, vij, s12ww, s12wzz, zi, zj, ww;

  /* Loop on the variables */

  int nb_neg = 0;
  for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      s12ww = s12wzz = 0.;
      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww)) continue;
        if (ww < 0.)
        {
          nb_neg++;
          continue;
        }

        /* Loop on the components */

        vi = vj = vij = 0;
        for (int icomp = 0; icomp < ncomp; icomp++)
        {
          zi = _get_IVAR(db, iech, ivar * ncomp + icomp);
          zj = _get_IVAR(db, iech, jvar * ncomp + icomp);
          if (!FFFF(zi) && !FFFF(zj))
          {
            vi += zi * zi;
            vj += zj * zj;
            vij += zi * zj;
          }
          else
          {
            vij = TEST;
            break;
          }
        }
        if (ABS(vi * vj) < 1.e-10 || FFFF(vij)) continue;
        vij = ABS(vij) / sqrt(vi * vj);

        s12ww += ww * ww;
        s12wzz += ww * ww * vij;
      }
      setVar((s12ww > 0) ? s12wzz / s12ww : 0., ivar, jvar);
      setVar((s12ww > 0) ? s12wzz / s12ww : 0., jvar, ivar);
    }

  if (nb_neg > 0)
    message("There were %d negative weights. They have been set to zero\n",
            nb_neg);
  return;
}

/****************************************************************************/
/*!
 **  Scale the variogram calculations
 **
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
void Vario::_variogram_scale(int idir)
{
  int nvar = getVariableNumber();

  /* Scale the experimental variogram quantities */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      for (int i = 0; i < getLagTotalNumber(idir); i++, ecr++)
      {
        int j = getDirAddress(idir, ivar, jvar, i, true, 0);
        if (getSwByIndex(idir, j) <= 0)
        {
          setHhByIndex(idir, j, TEST);
          setGgByIndex(idir, j, TEST);
        }
        else
        {
          setHhByIndex(idir, j, getHhByIndex(idir, j) / getSwByIndex(idir, j));
          if (getFlagAsym() && i < getLagNumber(idir))
            setHhByIndex(idir, j, -ABS(getHhByIndex(idir, j)));
          if (getCalcul() != ECalcVario::COVARIOGRAM)
            setGgByIndex(idir, j, getGgByIndex(idir, j) / getSwByIndex(idir, j));
        }
      }
    }

  // Process the variogram transformations

  if (getCalcul() == ECalcVario::TRANS1)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < getLagTotalNumber(idir); i++, ecr++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          int j0 = getDirAddress(idir, jvar, jvar, i, true, 0);
          setGgByIndex(idir, j,-getGgByIndex(idir, j) / getGgByIndex(idir, j0));
        }
      }
  }
  else if (getCalcul() == ECalcVario::TRANS2)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < getLagTotalNumber(idir); i++, ecr++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          int j0 = getDirAddress(idir, ivar, ivar, i, true, 0);
          setGgByIndex(idir, j, -getGgByIndex(idir, j) / getGgByIndex(idir, j0));
        }
      }
  }
  else if (getCalcul() == ECalcVario::BINORMAL)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        for (int i = 0; i < getLagTotalNumber(idir); i++, ecr++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          int j1 = getDirAddress(idir, ivar, ivar, i, true, 0);
          int j2 = getDirAddress(idir, jvar, jvar, i, true, 0);
          setGgByIndex(idir, j,
              getGgByIndex(idir, j) / sqrt(getGgByIndex(idir, j1) * getGgByIndex(idir, j2)));
        }
      }
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the data value (possibly after removing the global trend)
 **
 ** \return    The data value (or the residual)
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 **
 ** \remark  The trend removal only applies on the first variable
 ** \remark  Therefore, if applied on any variable rank other than 0,
 ** \remark  TEST is returned
 **
 *****************************************************************************/
double Vario::_get_IVAR(const Db *db, int iech, int ivar)
{
  double zz = db->getLocVariable(ELoc::Z, iech, ivar);
  if (FFFF(zz)) return (TEST);
  if (MODEL == nullptr) return (zz);
  if (ivar != 0) return (TEST);
  double drfval = model_drift_evaluate(0, MODEL, db, iech, 0, get_BETA().data());
  if (FFFF(drfval)) return (TEST);
  return (zz - drfval);
}

bool Vario::keepPair(int idir, SpaceTarget &T1, SpaceTarget &T2, double *dist)
{
  for (int ipt = 0, npt = getBiPtsNumberPerDirection(); ipt < npt; ipt++)
  {
    const ABiTargetCheck* bipts = getBipts(idir, ipt);
    if (! bipts->isOK(T1, T2)) return false;
    const BiTargetCheckGeometry* bigeom = dynamic_cast<const BiTargetCheckGeometry*>(bipts);
    if (bigeom != nullptr) *dist = bigeom->getDist();
  }
  return true;
}

/****************************************************************************
 **
 ** FUNCTION: st_identify_calcul_type
 **
 ** PURPOSE:  Identify the type of variogram calculation
 **
 ** IN_ARGS:  calcul_name  : Type of the variogram
 **
 ** REMARKS:  In case the calculation type is not identified,
 ** REMARKS:  the routine returns ECalcVario::UNDEFINED
 ** REMARKS:  The error message is produced internally
 **
 *****************************************************************************/
ECalcVario identifyVarioTypeByName(const String &calcul_name)

{
  ECalcVario calcul_type;

  if (!strcmp(calcul_name.c_str(), "vg"))
    calcul_type = ECalcVario::VARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "cov"))
    calcul_type = ECalcVario::COVARIANCE;
  else if (!strcmp(calcul_name.c_str(), "covnc"))
    calcul_type = ECalcVario::COVARIANCE_NC;
  else if (!strcmp(calcul_name.c_str(), "covg"))
    calcul_type = ECalcVario::COVARIOGRAM;
  else if (!strcmp(calcul_name.c_str(), "mado"))
    calcul_type = ECalcVario::MADOGRAM;
  else if (!strcmp(calcul_name.c_str(), "rodo"))
    calcul_type = ECalcVario::RODOGRAM;
  else if (!strcmp(calcul_name.c_str(), "poisson"))
    calcul_type = ECalcVario::POISSON;
  else if (!strcmp(calcul_name.c_str(), "general1"))
    calcul_type = ECalcVario::GENERAL1;
  else if (!strcmp(calcul_name.c_str(), "general2"))
    calcul_type = ECalcVario::GENERAL2;
  else if (!strcmp(calcul_name.c_str(), "general3"))
    calcul_type = ECalcVario::GENERAL3;
  else if (!strcmp(calcul_name.c_str(), "order4"))
    calcul_type = ECalcVario::ORDER4;
  else if (!strcmp(calcul_name.c_str(), "trans1"))
    calcul_type = ECalcVario::TRANS1;
  else if (!strcmp(calcul_name.c_str(), "trans2"))
    calcul_type = ECalcVario::TRANS2;
  else if (!strcmp(calcul_name.c_str(), "binormal"))
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
  return (calcul_type);
}

