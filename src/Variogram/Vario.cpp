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
#include "geoslib_define.h"
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
#include "Space/SpaceRN.hpp"
#include "Geometry/BiTargetCheckCode.hpp"
#include "Geometry/BiTargetCheckDate.hpp"
#include "Geometry/BiTargetCheckFaults.hpp"
#include "Geometry/BiTargetCheckGeometry.hpp"
#include "Polynomials/Hermite.hpp"

static int IECH1, IECH2, IDIRLOC;

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int VARWGT[4][5] = {{1, -1, 0, 0, 0}, {1, -2, 1, 0, 0}, {1, -3, 3, -1, 0}, {1, -4, 6, -4, 1}};
#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)    wt[IJDIR(ijvar, ipadir)]

/**
 * Build a Vario object by calculating the experimental variogram
 * @param varioparam VarioParam structure
 */
Vario::Vario(const VarioParam& varioparam)
    : AVario(),
      ASerializable(),
      _nVar(0),
      _varioparam(),
      _means(),
      _vars(),
      _flagSample(false),
      _db(nullptr),
      _sw(),
      _gg(),
      _hh(),
      _utilize(),
      _biPtsPerDirection(0),
      _bipts(),
      _flagAsym(false),
      _verbose(false),
      _flag_UK(false),
      _niter_UK(0),
      _variableNames(),
      _model(),
      _BETA(),
      _DRFDIAG(),
      _DRFXA(),
      _DRFGX(),
      _DRFTAB(),
      _DRFXGX()
{
  _varioparam = varioparam;
}

Vario::Vario(const Vario& r)
  : AVario(r)
  , ASerializable(r)
  , _nVar(r._nVar)
  , _varioparam(r._varioparam)
  , _means(r._means)
  , _vars(r._vars)
  , _flagSample(r._flagSample)
  , _db(r._db)
  , _sw(r._sw)
  , _gg(r._gg)
  , _hh(r._hh)
  , _utilize(r._utilize)
  , _biPtsPerDirection(r._biPtsPerDirection)
  , _flagAsym(r._flagAsym)
  , _verbose(r._verbose)
  , _flag_UK(r._flag_UK)
  , _niter_UK(r._niter_UK)
  , _variableNames(r._variableNames)
  , _model(r._model)
  , _BETA(r._BETA)
  , _DRFDIAG(r._DRFDIAG)
  , _DRFXA(r._DRFXA)
  , _DRFGX(r._DRFGX)
  , _DRFTAB(r._DRFTAB)
  , _DRFXGX(r._DRFXGX)
{
  for (int ipt = 0, npt = _getNBiPts(); ipt < npt; ipt++)
    _bipts.push_back(r._bipts[ipt]);
}

Vario& Vario::operator=(const Vario& r)
{
  if (this != &r)
  {
    AVario::operator=(r);
    ASerializable::operator=(r);
    _nVar = r._nVar;
    _varioparam = r._varioparam;
    _means = r._means;
    _vars  = r._vars;
    _flagSample = r._flagSample;
    _db = r._db;
    _sw = r._sw;
    _gg = r._gg;
    _hh = r._hh;
    _utilize = r._utilize;
    _biPtsPerDirection = r._biPtsPerDirection;
    _flagAsym = r._flagAsym;
    _verbose = r._verbose;
    _flag_UK = r._flag_UK;
    _niter_UK = r._niter_UK;
    _variableNames = r._variableNames;
    _model = r._model;
    _BETA = r._BETA;
    _DRFDIAG = r._DRFDIAG;
    _DRFXA = r._DRFXA;
    _DRFGX = r._DRFGX;
    _DRFTAB = r._DRFTAB;
    _DRFXGX = r._DRFXGX;

    for (int ipt = 0, npt = _getNBiPts(); ipt < npt; ipt++)
      _bipts.push_back(r._bipts[ipt]);
  }
  return *this;
}

Vario::~Vario()
{
  for (int ipt = 0, npt = _getNBiPts(); ipt < npt; ipt++)
    delete _bipts[ipt];
}

Vario* Vario::create(const VarioParam& varioparam)
{
  return new Vario(varioparam);
}

Vario* Vario::createFromNF(const String& neutralFilename, bool verbose)
{
  std::ifstream is;
  VarioParam varioparam = VarioParam();
  Vario* vario = new Vario(varioparam);
  bool success = false;
  if (vario->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  vario->deserialize(is, verbose);
  }
  if (! success)
  {
    delete vario;
    vario = nullptr;
  }
  return vario;
}

Vario* Vario::computeFromDb(const VarioParam& varioparam,
                            Db* db,
                            const ECalcVario& calcul,
                            bool flag_sample,
                            bool verr_mode,
                            Model *model,
                            int niter_UK,
                            bool verbose)
{
  Vario* vario = nullptr;
  vario = new Vario(varioparam);
  if (vario->compute(db, calcul, flag_sample, verr_mode, model, niter_UK, verbose))
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
  if (vario->regularizeFromModel(model, ext, ndisc, angles))
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
    delete varioZ;
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
  for (int ipt = 0, npt = _getNBiPts(); ipt < npt; ipt++)
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
  for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
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
  if (getNDir() <= 0)
  {
    messerr("The 'varioParam' argument must have some Direction defined");
    return 1;
  }

  // Preparation
  setCalcul(calcul);
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
                   bool flag_sample,
                   bool verr_mode,
                   const Model *model,
                   int niter_UK,
                   bool verbose)
{
  _db   = db;
  _nVar = _db->getNLoc(ELoc::Z);
  if (_nVar <= 0)
  {
    messerr(
      "You need some Variable defined (Z locator) to calculate variogram");
    return 1;
  }

  if (prepare(calcul)) return 1;

  if (_compute(_db, flag_sample, verr_mode, model, niter_UK, verbose))
  {
    messerr("Error when calculating the Variogram");
    return 1;
  }
  return 0;
}

int Vario::computeIndic(Db *db,
                        const ECalcVario& calcul,
                        bool flag_sample,
                        bool verr_mode,
                        const Model *model,
                        int niter_UK,
                        bool verbose,
                        int nfacmax)
{
  _db = db;
  int nvar = _db->getNLoc(ELoc::Z);
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
  if (_compute(_db, flag_sample, verr_mode, model, niter_UK, verbose))
  {
    messerr("Error when calculating the Variogram of Indicators");
    return 1;
  }

  // Delete the Indicators (created locally)
  _db->deleteColumnsByLocator(ELoc::Z);
  _db->setLocatorByUID(iatt, ELoc::Z, 0);

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
  int nvar_in = vario_in.getNVar();
  int ndir_in = vario_in._varioparam.getNDir();

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
    setCalculByName("vg");
  }

  // Reset Mean and variance arrays (only if variable number has been modified)
  if (_nVar != nvar_in)
  {
    if (! vario_in.getMeans().empty())
    {
      _means.resize(_nVar);
      for (int ivar = 0; ivar < _nVar; ivar++)
        setMean(vario_in.getMean(selvars[ivar]), ivar);
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
    _vars  = vario_in.getVars();
  }

  // Add the directions

  internalDirectionResize(ndir_in,true);

  for (int idir0 = 0; idir0 < ndir; idir0++)
  {
    int idir = seldirs[idir0];
    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {
      for (int ivar0 = 0; ivar0 < _nVar; ivar0++)
        for (int jvar0 = 0; jvar0 < _nVar; jvar0++)
        {
          int ivar = selvars[ivar0];
          int jvar = selvars[jvar0];

          int iadto = getDirAddress(idir0,ivar0,jvar0,ilag);

          if (! flagMakeSym)
          {
            int iadfrom = vario_in.getDirAddress(idir,ivar,jvar,ilag,false,0);
            _sw[idir][iadto] = vario_in.getSwByIndex(idir0,iadfrom);
            _gg[idir][iadto] = vario_in.getGgByIndex(idir0,iadfrom);
            _hh[idir][iadto] = vario_in.getHhByIndex(idir0,iadfrom);
            _utilize[idir][iadto] = vario_in.getUtilizeByIndex(idir0,iadfrom);
          }
          else
          {
            int iadf1 = vario_in.getDirAddress(idir,ivar,jvar,ilag,false,-1);
            int iadf2 = vario_in.getDirAddress(idir,ivar,jvar,ilag,false,1);
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
  if (getNVar() != 1)
  {
    messerr("The function 'transformZToY' is restricted to Monovariate Variogram");
    return 1;
  }

  /* Loop on the directions of the variogram */

  double cvv = anam->getVariance();
  for (int idir = 0; idir < getNDir(); idir++)
  {
    /* Loop on the lags */

    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {
      // TODO. GG must be a variogram of Zv -> Cv(h)
      setGgByIndex(idir,ilag,1. - anamH->invertVariance(cvv-getGgByIndex(idir, ilag)));
    }
  }

  // Modify the variance array
  setVar(1., 0,  0);

  delete anamH;
  
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
  if (getNVar() != 1)
  {
    messerr("The function 'transformYToZ' is restricted to Monovariate Variogram");
    return 1;
  }

  /* Loop on the directions of the variogram */

  double c0 = anam_hermite->computeVariance(1.);
  for (int idir = 0; idir < getNDir(); idir++)
  {
    /* Loop on the lags */

    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {
      double chh = 1. - getGg(idir, 0, 0, ilag, false);
      double var = anam_hermite->computeVariance(chh);
      setGg(idir, 0, 0, ilag, c0 - var);
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
int Vario::regularizeFromModel(const Model &model,
                               const VectorDouble &ext,
                               const VectorInt &ndisc,
                               const VectorDouble &angles,
                               const CovCalcMode *mode,
                               bool asCov)
{
  int ndim = model.getNDim();
  int nvar = model.getNVar();

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

  for (int idir = 0; idir < getNDir(); idir++)
  {

    /* Loop on the number of lags */

    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {

      // Calculate the shift vector

      double dist = ilag * getDPas(idir);
      VectorDouble shift(ndim);
      for (int idim = 0; idim < ndim; idim++)
        shift[idim] = dist * getCodir(idir, idim);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          double value = model.evalCvvShift(ext, ndisc, shift, angles, ivar, jvar, mode);
          if (! asCov) value = getVar(ivar, jvar) - value;
          int iad = getDirAddress(idir, ivar, jvar, ilag, false, 0);
          setGgByIndex(idir, iad, value);
          setHhByIndex(idir, iad, dist);
          setSwByIndex(idir, iad, 1);
        }
    }
  }
  return 0;
}

MatrixSquareGeneral Vario::_evalAverageDbIncr(Model *model,
                                              const Db &db,
                                              const VectorDouble &incr,
                                              const CovCalcMode *mode) const
{
  int nvar = getNVar();
  int nech = db.getNSample(true);
  int ndim = getNDim();
  int norme = nech * nech;

  MatrixSquareGeneral mat(nvar);
  VectorDouble dd(ndim);
  MatrixSquareGeneral covtab(nvar);

  for (int iech = 0; iech < nech; iech++)
  {
    if (! db.isActive(iech)) continue;
    for (int jech = 0; jech < nech; jech++)
    {
      if (! db.isActive(jech)) continue;

      // Calculate the distance between the two samples
      db.getDistanceVecInPlace(iech, jech, dd);
      if (! incr.empty()) VH::addInPlace(dd, incr);

      // Evaluate the covariance matrix between two samples
      model->evaluateMatInPlace(nullptr, dd, covtab, false, 1., mode);
    }
  }

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.setValue(ivar, jvar, covtab.getValue(ivar, jvar) / norme);
  return mat;
}

/****************************************************************************/
/*!
 **  Calculate the regularized model as an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  db        Db discretization structure
 ** \param[in]  mode      CovCalcMode structure
 **
 *****************************************************************************/
int Vario::regularizeFromDbGrid(Model* model,
                                const Db& db,
                                const CovCalcMode* mode)
{
  int nvar = model->getNVar();
  setNVar(nvar);
  internalVariableResize();
  internalDirectionResize();

  /* Calculate the Cvv (for a zero-shift) */

  MatrixSquareGeneral c00tab = _evalAverageDbIncr(model, db, VectorDouble(), mode);

  /* Initialize the variance array */

  setVars(c00tab.getValues());

  /* Loop on the directions */

  for (int idir = 0; idir < getNDir(); idir++)
  {

    /* Loop on the number of lags */

    for (int ilag = 0; ilag < getNLag(idir); ilag++)
    {
      double dist = ilag * getDPas(idir);
      VectorDouble dd = getCodirs(idir);
      VH::multiplyConstant(dd, dist);

      MatrixSquareGeneral covtab = _evalAverageDbIncr(model, db, dd, mode);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          int iad = getDirAddress(idir, ivar, jvar, ilag, false, 0);
          setGgByIndex(idir, iad, c00tab.getValue(ivar,jvar) - covtab.getValue(ivar,jvar));
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
    _nVar = _db->getNLoc(ELoc::Z);
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
  if (ndir <= 0) ndir = getNDir();
  _sw.resize(ndir);
  _gg.resize(ndir);
  _hh.resize(ndir);
  _utilize.resize(ndir);

  if (flagDirs)
    for (int idir = 0; idir < getNDir(); idir++)
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

  sstr << _elemString(strfmt);
  if (getCalcul() == ECalcVario::UNDEFINED) return sstr.str();

  sstr << "Number of variable(s)       = " << _nVar << std::endl;

  // Print the environment

  sstr << _varioparam.toStringMain(strfmt);

  // Print the variable names

  if (! _variableNames.empty())
    sstr << "Variable(s)                 = " << _variableNames.toString() << std::endl;

  // Print the variance matrix

  sstr << toMatrix("Variance-Covariance Matrix",VectorString(),VectorString(),
                    0,_nVar,_nVar,getVars());

  /* Loop on the directions (only if the resulting arrays have been defined) */

  if (!_sw.empty())
  {
    for (int idir = 0; idir < getNDir(); idir++)
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

      for (int i = 0; i < getNLagTotal(idir); i++)
      {
        int j = getDirAddress(idir, ivar, jvar, i, true, 0);
        if (_sw[idir][j] <= 0) continue;
        int rank = (!getFlagAsym()) ? i : i - getNLag(idir);
        sstr << toInt(rank);
        sstr << toDouble(_sw[idir][j]);
        sstr << toDouble(_hh[idir][j]);
        sstr << toDouble(_gg[idir][j]);
        sstr << std::endl;
      }
    }
  return sstr.str();
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

MatrixSquareSymmetric Vario::getVarMatrix() const
{
  MatrixSquareSymmetric mat(_nVar);
  for (int ivar = 0; ivar < _nVar; ivar++)
    for (int jvar = 0; jvar < _nVar; jvar++)
      mat.setValue(ivar, jvar, getVar(ivar, jvar));
  return mat;
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

void Vario::setGgByIndex(int idir, int i, double gg, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _gg[idir][i] = gg;
}

void Vario::setHhByIndex(int idir, int i, double hh, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _hh[idir][i] = hh;
}

void Vario::setSwByIndex(int idir, int i, double sw, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _sw[idir][i] = sw;
}

void Vario::setUtilizeByIndex(int idir, int i, double utilize, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _utilize[idir][i] = utilize;
}

void Vario::setSw(int idir, int ivar, int jvar, int ilag, double sw, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
  }
  int iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
  if (IFFFF(iad)) return;
  _sw[idir][iad] = sw;
}

void Vario::setHh(int idir, int ivar, int jvar, int ilag, double hh, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
  }
  int iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
  if (IFFFF(iad)) return;
  _hh[idir][iad] = hh;
}

void Vario::setGg(int idir, int ivar, int jvar, int ilag, double gg, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
  }
  int iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
  if (IFFFF(iad)) return;
  _gg[idir][iad] = gg;
}

void Vario::setUtilize(int idir, int ivar, int jvar, int ilag, double utilize, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
  }
  int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
  if (IFFFF(iad)) return;
  _utilize[idir][iad] = utilize;
}

void Vario::updateSwByIndex(int idir, int i, double sw, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _sw[idir][i] += sw;
}

void Vario::updateHhByIndex(int idir, int i, double hh, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _hh[idir][i] += hh;
}

void Vario::updateGgByIndex(int idir, int i, double gg, bool flagCheck)
{
  if (! _isAddressValid(idir, i, flagCheck)) return;
  _gg[idir][i] += gg;
}

double Vario::getGg(int idir,
                    int ivar,
                    int jvar,
                    int ilag,
                    bool asCov,
                    bool flagNorm) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
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

double Vario::getHh(int idir, int ivar, int jvar, int ilag) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
  if (IFFFF(iad)) return TEST;
  return _hh[idir][iad];
}

double Vario::getSw(int idir, int ivar, int jvar, int ilag) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
  if (IFFFF(iad)) return TEST;
  return _sw[idir][iad];
}

double Vario::getUtilize(int idir, int ivar, int jvar, int ilag) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
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

  int nlag = getNLag(idir);
  vec.resize(3);
  for (int i = 0; i < 3; i++) vec[i].resize(nlag);

  for (int ilag = 0 ; ilag < nlag; ilag++)
  {
    int iad = getDirAddress(idir,ivar,jvar,ilag,true,0);
    vec[0][ilag] = _sw[idir][iad];
    vec[1][ilag] = _hh[idir][iad];
    vec[2][ilag] = _gg[idir][iad];
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
  int nlag = getNLag(idir);

  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag-1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
      {
        double val = _gg[idir][iad];
        if (asCov && !getFlagAsym()) val = c0 - val;
        if (flagNorm) val /= c0;
        gg.push_back(val);
      }
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
    {
      double val = _gg[idir][iad];
      if (asCov && !getFlagAsym()) val = c0 - val;
      if (flagNorm) val /= c0;
      gg.push_back(val);
    }
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
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
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
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
  int nlag = getNLag(idir);
  if (nlag != (int) gg.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      setGg(idir, ivar, jvar, ilag, gg[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setGg(idir, ivar, jvar, 0, gg[iad]);
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      setGg(idir, ivar, jvar, ilag, gg[iad]);
    }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      setGg(idir, ivar, jvar, ilag, gg[iad]);
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
  int nlag = getNLag(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
        hh.push_back(_hh[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
      hh.push_back(_hh[idir][iad]);
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
        hh.push_back(_hh[idir][iad]);
    }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
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

  int nlag = getNLag(idir);
  if (nlag != (int) hh.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      setHh(idir, ivar, jvar, ilag, hh[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setHh(idir, ivar, jvar, 0, hh[iad]);
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      setHh(idir, ivar, jvar, ilag, hh[iad]);
    }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      setHh(idir, ivar, jvar, ilag, hh[iad]);
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
  int nlag = getNLag(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
        sw.push_back(_sw[idir][iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
      sw.push_back(_sw[idir][iad]);
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
        sw.push_back(_sw[idir][iad]);
    }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
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
  int nlag = getNLag(idir);
  if (nlag != (int) sw.size()) return;

  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
      setSw(idir, ivar, jvar, ilag, sw[iad]);
    }
    iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
    setSw(idir, ivar, jvar, 0, sw[iad]);
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
      setSw(idir, ivar, jvar, ilag, sw[iad]);
    }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      setSw(idir, ivar, jvar, ilag, sw[iad]);
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
  int nlag = getNLag(idir);
  int iad;
  if (_flagAsym)
  {
    for (int ilag = nlag - 1; ilag >= 0; ilag--)
     {
       iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
       if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
         utilize.push_back(_utilize[idir][iad]);
     }
     iad = getDirAddress(idir, ivar, jvar, 0, false, 0);
     if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
       utilize.push_back(_utilize[idir][iad]);
     for (int ilag = 0; ilag < nlag; ilag++)
     {
       iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
       if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
        utilize.push_back(_utilize[idir][iad]);
     }
  }
  else
  {
    for (int ilag = 0; ilag < nlag; ilag++)
    {
      iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
      if (!IFFFF(iad) && (!compress || _sw[idir][iad] > 0.))
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
  int nlag = getNLag(idir);
  int count;
  if (_flagAsym) return ITEST;
  int iad = getDirSize(idir) - 1;
  count = 0;
  for (int ilag = 0; ilag < nlag && count < shift; ilag++)
  {
    iad = getDirAddress(idir, ivar, jvar, ilag, true, 0);
    if (IFFFF(iad)) continue;
    if (! isZero(_sw[idir][iad]) && ! isZero(_hh[idir][iad])) count++;
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
                         int ilag,
                         bool flag_abs,
                         int sens,
                         bool flagCheck) const
{
  if (!_isDirectionValid(idir, flagCheck)) return ITEST;
  if (!_isVariableValid(ivar, flagCheck))  return ITEST;
  if (!_isVariableValid(jvar, flagCheck)) return ITEST;

  int rank;

  /* Get the order of the variables */

  if (ivar > jvar)
    rank = ivar * (ivar + 1) / 2 + jvar;
  else
    rank = jvar * (jvar + 1) / 2 + ivar;

  /* Get the position in the array */

  int iad = 0;
  if (flagCheck)
  {
    const DirParam dirparam = _varioparam.getDirParam(idir);
    if (!dirparam.isLagValid(ilag, getFlagAsym(), flagCheck)) return ITEST;
  }

  if (! getFlagAsym())
  {
    iad = ilag;
  }
  else
  {
    if (flag_abs)
    {
      iad = ilag;
    }
    else
    {
      int nlag = getNLag(idir);
      switch (sens)
      {
        case 1:
          iad = nlag + ilag + 1;
          break;

        case -1:
          iad = nlag - ilag - 1;
          break;

        case 0:
          iad = nlag;
          break;

        default:
          break;
      }
    }
  }
  iad += rank * getNLagTotal(idir);
  return (iad);
}

bool Vario::_isVariableValid(int ivar, bool flagCheck) const
{
  if (!flagCheck) return true;
  return checkArg("Variable Index", ivar, _nVar);
}

bool Vario::_isBivariableValid(int ijvar, bool flagCheck) const
{
  if (!flagCheck) return true;
  return checkArg("Multivariate Index", ijvar, _nVar * _nVar);
}

bool Vario::_isDirectionValid(int idir, bool flagCheck) const
{
  if (!flagCheck) return true;
  return checkArg("Direction Index", idir, getNDir());
}

bool Vario::_isAddressValid(int idir, int i, bool flagCheck) const
{
  if (! flagCheck) return true;
  if (! _isDirectionValid(idir)) return false;
  const DirParam dirparam = _varioparam.getDirParam(idir);
  return (i >= 0 && i < getDirSize(idir));
}

bool Vario::_deserialize(std::istream& is, bool /*verbose*/)
{
  int flag_calcul = 0;
  int ndim = 0;
  int nvar = 0;
  int ndir = 0;
  int nlag = 0;
  int opt_code = 0;
  int flag_regular = 0;

  double dlag = 0.;
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

  // Reading the calculation flag
  ret = ret && _recordRead<int>(is, "Calculation Flag", flag_calcul);

  // Reading the variable names
  _variableNames.resize(nvar, "Unknown");
  if (flag_calcul == 2)
  {
    // Reading the variable names
    for (int ivar = 0; ivar < nvar; ivar++)
      ret = ret && _recordRead<String>(is, "Variable Name", _variableNames[ivar]);
  }

  /* Read the variances (optional) */

  vars.resize(nvar * nvar);
  if (flag_calcul)
  {
    int ecr = 0;
    for (int ivar = 0; ret && ivar < nvar; ivar++)
      for (int jvar = 0; ret && jvar < nvar; jvar++, ecr++)
        ret = ret && _recordRead<double>(is, "Variance", vars[ecr]);
  }
  if (! ret) return ret;

  /* Initialize the variogram structure */

  _nVar = nvar;
  internalDirectionResize(ndir,false);
  setVars(vars);
  setCalculByName("vg"); // TODO: read this information from NF file and treat accordingly
  setScale(scale);
  int isDefinedForGrid = 0;

  /* Reading the variogram calculation directions */

  for (int idir = 0; ret && idir < ndir; idir++)
  {
    ret = ret && _recordRead<int>(is, "Regular Variogram Calculation", flag_regular);
    ret = ret && _recordRead<int>(is, "Number of Variogram Lags", nlag);
    ret = ret && _recordRead<int>(is, "Variogram Code Option", opt_code);
    ret = ret && _recordRead<double>(is, "Tolerance on Code", tolcode);
    ret = ret && _recordRead<double>(is, "Lag Value", dlag);
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

    auto space  = SpaceRN::create(ndim);
    DirParam dirparam = DirParam(nlag, dlag, toldis, tolang, opt_code, 0,
                                 TEST, TEST, tolcode, VectorDouble(), codir, TEST,
                                 space);
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

// TODO: add the type of calculation which has been performed
bool Vario::_serialize(std::ostream& os, bool /*verbose*/) const
{
  double value;
  static int flag_calcul = 2;

  /* Write the Vario structure */

  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", _varioparam.getNDim());
  ret = ret && _recordWrite<int>(os, "Number of variables", getNVar());
  ret = ret && _recordWrite<int>(os, "Number of directions", getNDir());
  ret = ret && _recordWrite<double>(os, "Scale", _varioparam.getScale());
  ret = ret && _recordWrite<int>(os, "Calculation Flag", flag_calcul);

  // Dump the variable names

  ret = ret && _commentWrite(os, "Variable Names");
  for (int ivar = 0; ivar < getNVar(); ivar++)
  {
    if (ivar < (int) _variableNames.size())
      ret = ret && _recordWrite<String>(os, "", _variableNames[ivar]);
    else
      ret = ret && _recordWrite<String>(os, "", "Unknown");
  }
  ret = ret && _commentWrite(os, "");

  /* Dumping the Variances */

  if (flag_calcul)
  {
    ret = ret && _commentWrite(os, "Variance");
    for (int ivar = 0; ret && ivar < getNVar(); ivar++)
    {
      for (int jvar = 0; ret && jvar < getNVar(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getVar(ivar,jvar));
      ret = ret && _commentWrite(os, "");
    }
  }

  /* Loop on the directions */

  for (int idir = 0; ret && idir < getNDir(); idir++)
  {
    const DirParam dirparam = _varioparam.getDirParam(idir);
    ret = ret && _commentWrite(os, "Direction characteristics");
    ret = ret && _recordWrite<int>(os, "Regular lags", dirparam.getFlagRegular());
    ret = ret && _recordWrite<int>(os, "Number of lags", dirparam.getNLag());
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
  if (ivar < 0 || ivar >= getNVar())
  {
    bounds[0] = 0;
    bounds[1] = getNVar();
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
  if (idir < 0 || idir >= getNDir())
  {
    bounds[0] = 0;
    bounds[1] = getNDir();
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
  return (getNLagTotal(idir) * _nVar * (_nVar + 1) / 2);
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
    _nVar = db->getNLoc(ELoc::Z);
    return 0;
  }
  if (!_means.empty())
  {
    _nVar = static_cast<int>(_means.size());
    return 0;
  }

  messerr("Cannot determine the Number of Variables from arguments");
  return 1;
}

int Vario::getNLagTotal(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  int nlag = getNLag(idir);
  return ((_flagAsym) ? 2 * nlag + 1 : nlag);
}

void Vario::setCalculByName(const String& calcul_name)
{
  AVario::setCalculByName(calcul_name);
  _setFlagAsym();
}

double Vario::getMaximumDistance() const
{
  double distmax = 0.;
  for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
  {
    double dist = getMaximumDistance(idir);
    if (dist > distmax) distmax = dist;
  }
  return distmax;
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
    for (int idir = 0; idir < getNDir(); idir++)
    {
      _varioparam.setDPas(idir, dbgrid);
    }
  }
  else
  {
    for (int idir = 0; idir < getNDir(); idir++)
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
  return (ivar == jvar || ! getFlagAsym());
}

bool Vario::drawOnlyPositiveY(int ivar, int jvar) const
{
  return (ivar == jvar && ! getFlagAsym());
}

VectorDouble Vario::getGgs(int idir, int ivar, int jvar, const VectorInt& ilag) const
{
  VectorDouble values;
  if (ilag.empty()) return values;
  if (! _isDirectionValid(idir)) return values;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  for (int i = 0; i < (int) ilag.size(); i++)
  {
    if (ilag[i] >= 0 && ilag[i] < getDirSize(idir)) values.push_back(getGg(idir,ivar,jvar,ilag[i]));
  }
  return values;
}

VectorDouble Vario::setGgs(int idir, int ivar, int jvar, const VectorInt& ilag, const VectorDouble& values)
{
  if (ilag.empty()) return values;
  if (values.empty()) return values;
  if (! _isDirectionValid(idir)) return values;
  const DirParam dirparam = _varioparam.getDirParam(idir);

  for (int i = 0; i < (int) ilag.size(); i++)
  {
    if (ilag[i] >= 0 && ilag[i] < getDirSize(idir) && i < (int) values.size())
      setGg(idir,ivar,jvar,ilag[i], values[i]);
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
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 ** \param[in]  model        Model structure (triggers the KU option)
 ** \param[in]  niter_UK     Number of iteration for bias correction
 ** \param[in]  verbose      Verbose flag
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
int Vario::_compute(Db *db,
                    int flag_sample,
                    int verr_mode,
                    const Model *model,
                    int niter_UK,
                    bool verbose)
{
  if (model != nullptr) _model = model->clone();
  _verbose = verbose;
  int norder = _get_generalized_variogram_order();

  // Preliminary checks
  if (! _isCompatible(db)) return 1;

  if (_model != nullptr && _model->isDriftDifferentDefined(VectorInt()))
  {
    _flag_UK = true;
    _driftManage(db);
    _niter_UK = niter_UK;
    if (_niter_UK != 0)
    {
      if (isDefinedForGrid())
      {
        messerr("Drift Bias correction is not coded in the case of Grid");
        return 1;
      }
      // If the Model does not contain any covariance, add some defaulted basic structures, i.e.:
      // - a nugget effect
      // - an exponential covariance (with initial range set to 1 a,d sill to 1)
      // - a spherical covariance (with initial range set to 2, and sill to 1)
      int ncov = _model->getNCov();
      if (ncov <= 0)
      {
        _model->addCovFromParam(ECov::NUGGET);
        _model->addCovFromParam(ECov::EXPONENTIAL,  1.,  1.);
        _model->addCovFromParam(ECov::SPHERICAL, 2., 1.);
      }
    }
  }

  // Save the variable names
  int nvar = getNVar();
  _variableNames.resize(nvar, "Unknown");
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (ivar < db->getNLoc(ELoc::Z))
      setVariableName(ivar, db->getNameByLocator(ELoc::Z, ivar));
  }

  // Dispatch
  if (isDefinedForGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
    if (dbgrid == nullptr)
    {
      messerr("'Vario' is defined for Grid but 'db' is not organized as a grid");
      return 1;
    }
    if (norder > 0) return _calculateGenOnGrid(dbgrid, norder);
    return _calculateOnGrid(dbgrid);
  }
  if (norder > 0) return _calculateGenOnLine(db, norder);
  return _calculateGeneral(db, flag_sample, verr_mode);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  flag_sample  calculate the variogram per sample
 ** \param[in]  verr_mode    Mode of variogram correction (1, 2 or 3)
 **
 ** \note The number of iterations in KU is given can be updated using the
 ** \note keypair technique with name "KU_Niter".
 **
 *****************************************************************************/
int Vario::_calculateGeneral(Db *db,
                             int flag_sample,
                             int verr_mode)
{
  bool flag_verr = false;
  Vario_Order* vorder = (Vario_Order*) NULL;

  /* Particular case of Transitive Covariogram */
  /* It is only coded in the by_sample case and uses the regression technique */

  if (getCalcul() == ECalcVario::COVARIOGRAM) flag_sample = 1;

  /* Auxiliary check for Variance Measurement Error */

  if (db->getNLoc(ELoc::V) > 0 && verr_mode > 0)
  {
    vorder = vario_order_manage(1, 1, 0, vorder);
    flag_verr = true;
  }

  // Auxiliary check for Drift removal. This is triggered only if the drift
  // contains at least one drift function different from Universality condition

  if (_flag_UK)
  {
    if (vorder == (Vario_Order*) NULL)
      vorder = vario_order_manage(1, 1, 0, vorder);
  }

  /* Complementary checks */

  if (flag_verr && _flag_UK)
  {
    messerr("These two options are incompatible");
    messerr("- Correction for the Variance of Error Measurements");
    messerr("- Correction for bias when removing the Drift");
    return 1;
  }
  if (flag_verr || _flag_UK)
  {
    if (flag_sample)
    {
      messerr("The special Variogram option is incompatible with flag.sample");
      return 1;
    }
    if (!db->isNVarComparedTo(1)) return 1;
  }

  /* Evaluate the drift coefficients */

  if (_flag_UK)
  {
    if (_driftEstimateCoefficients(db)) return 1;
  }

  /* Update the global statistics */

  _getStatistics(db);

  /* Loop on the directions to evaluate */

  VectorInt rindex = db->getSortArray();
  for (int idir = 0; idir < getNDir(); idir++)
  {
    if (!flag_sample)
    {
      if (_calculateGeneralSolution1(db, idir, rindex.data(), vorder)) return 1;
    }
    else
    {
      if (_calculateGeneralSolution2(db, idir, rindex.data())) return 1;
    }

    if (vorder != (Vario_Order*) NULL)
      _calculateFromGeometry(db, idir, vorder);
  }

  /* Posterior calculations when presence of Variance of Measurement errors */

  if (flag_verr)
  {
    for (int idir = 0; idir < getNDir(); idir++)
    {
      if (_updateVerr(db, idir, vorder, verr_mode)) return 1;
    }
  }

  /* Posterior update when filtering the bias attached to drift removal */

  if (_flag_UK && _niter_UK != 0)
  {
    if (_updateUK(db, vorder)) return 1;
  }

  /* Set the error return code */

  vario_order_manage(-1, 1, 0, vorder);
  return 0;
}

/****************************************************************************/
/*!
 **  Update the Variogram of Residuals when Drift has been removed
 **
 ** \return  Error returned code
 **
 ** \param[in]  db         Db description
 ** \param[in]  vorder     Vario_Order structure
 **
 *****************************************************************************/
int Vario::_updateUK(Db *db, Vario_Order *vorder)
{
  Option_VarioFit optvar;
  Option_AutoFit mauto;
  int ifirst, ilast;
  Constraints constraints;

  // Do not allow reducing the number of covariances in the Model during iterations
  optvar.setFlagNoreduce(true);

  // Particular case when 'niter_UK' is negative, we want to exhibit the bias alone:
  // - the Model is supposed to be provided in input
  // - only one iteration is performed

  if (_niter_UK < 0)
  {
    // Calculate the global bias correction terms
    _calculateBiasGlobal(db);

    // Loop on the directions
    for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
    {

      // Loop on the lags
      for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
      {
        vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);
        if (ifirst > ilast) continue;
        _calculateBiasLocal(db, idir, ilag, vorder, ifirst, ilast);
      }
    }
    return 0;
  }

  // Loop on the iterations (when 'niter_UK' is positive
  for (int iter = 0; iter < _niter_UK; iter++)
  {

    // Perform the Automatic structure recognition
    if (model_auto_fit(this, _model, false, mauto, constraints, optvar)) return 1;

    // Calculate the global bias correction terms
    _calculateBiasGlobal(db);

    // Optional printout
    if (_verbose)
    {
      message("Drift removal at iteration #%d/%d\n", iter + 1, _niter_UK);
      _model->display();
      print_matrix("Drift Coefficients Matrix", 0, 1, _DRFXGX.getNRows(),
                   _DRFXGX.getNCols(), NULL, _DRFXGX.getValues().data());
    }

    // Loop on the directions
    for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
    {

      // Loop on the lags
      for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
      {
        vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);
        if (ifirst > ilast) continue;
        _calculateBiasLocal(db, idir, ilag, vorder, ifirst, ilast);
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
 ** \param[in]  idir      Rank of the current direction
 ** \param[in]  ilag      Rank of the current lag
 ** \param[in]  vorder    Vario_Order structure
 ** \param[in]  ifirst    Rank of the first lag
 ** \param[in]  ilast     Rank of the last lag
 **
 ** \remarks: When '_niter_UK' < 0: the bias is stored instead of the variogram
 **
 *****************************************************************************/
void Vario::_calculateBiasLocal(Db *db,
                                int idir,
                                int ilag,
                                Vario_Order *vorder,
                                int ifirst,
                                int ilast)
{
  int iech, jech;
  double dist;

  /* Calculate the first corrected term */

  double tot0 = 0.;
  double tot1 = 0.;
  double tot2 = 0.;
  double totnum = 0.;
  for (int ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
    double v1 = _getIVAR(db, iech, 0);
    double v2 = _getIVAR(db, jech, 0);
    if (FFFF(v1) || FFFF(v2)) continue;

    int iiech = _getRelativeSampleRank(db, iech);
    int jjech = _getRelativeSampleRank(db, jech);

    double diff = v1 - v2;
    tot0 += diff * diff;
    tot1 += _getBias(iiech, jjech);
    tot2 += (_DRFDIAG[iiech] + _DRFDIAG[jjech]) / 2.;
    totnum += 1.;
  }
  tot0 /= 2.;

  if (totnum > 0.)
  {
    double oldval = tot0 / totnum;
    double newval = (tot2 - tot1) / totnum;
    if (_niter_UK > 0) newval += oldval;
    setGg(idir, 0, 0, ilag, newval);
  }
}

/****************************************************************************/
/*!
 **  Calculate the global bias terms
 **
 ** \param[in]  db        Db description
 **
 *****************************************************************************/
void Vario::_calculateBiasGlobal(Db *db)
{
  double covtab, value;

  /* Initializations */

  int nbfl = _model->getNDrift();
  int ndim = _model->getNDim();
  int nech = db->getNSampleActiveAndDefined(0);
  VectorDouble d1(ndim,0.);

  /* Calculate the c00 term */

  double c00 = _model->evaluateOneGeneric(nullptr, d1);

  /* Calculate the term: G %*% X */

  int iiech = 0;
  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    for (int il = 0; il < nbfl; il++)
    {
      value = 0;
      int jjech = 0;
      for (int jech = 0; jech < db->getNSample(); jech++)
      {
        if (!db->isActiveAndDefined(jech, 0)) continue;
        for (int idim = 0; idim < ndim; idim++)
          d1[idim] = db->getDistance1D(iech, jech, idim);
        covtab = _model->evaluateOneGeneric(nullptr, d1);
        value += (c00 - covtab) * _DRFTAB.getValue(jjech, il);
        jjech++;
      }
      _DRFGX.setValue(iiech, il, value);
    }
    iiech++;
  }

  /* Calculate the term: t(X) %*% G %*% X */

  for (int il = 0; il < nbfl; il++)
    for (int jl = 0; jl < nbfl; jl++)
    {
      value = 0;
      for (int iech = 0; iech < nech; iech++)
        value += _DRFGX.getValue(iech, il) * _DRFTAB.getValue(iech, jl);
      _DRFXGX.setValue(il,jl,value);
    }

  /* Calculate the term: diag(bias) */

  iiech = 0;
  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    _DRFDIAG[iiech] = _getBias(iiech, iiech);
    iiech++;
  }
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram along lines
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  norder       Order of the generalized variogram
 **
 *****************************************************************************/
int Vario::_calculateGenOnLine(Db *db, int norder)
{
  /* Preliminary checks */

  if (getNVar() != 1)
  {
    messerr("The generalized variogram requires a single variable");
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

  _getStatistics(db);

  /* Loop on the directions to evaluate */

  for (int idir = 0; idir < getNDir(); idir++)
    _calculateOnLineSolution(db, idir, norder);

  return 0;
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
void Vario::_calculateOnLineSolution(Db *db, int idir, int norder)
{
  SpaceTarget T1(getSpace(),false);
  SpaceTarget T2(getSpace(),false);
  int jech, keep;
  double value, zz;

  int nech     = db->getNSample();
  int nlag     = getNLag(idir);
  int nvar     = getNVar();
  double dist0 = 0.;
  double dist  = 0.;
  bool hasSel  = db->hasLocVariable(ELoc::SEL);

  /* Loop on the first point */

  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    db->getSampleAsSTInPlace(iech, T1);

    for (int ilag = 1; ilag < nlag; ilag++)
    {
      value = _getIVAR(db, iech, 0);
      if (FFFF(value)) break;
      dist0 = 0.;

      for (int iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
      {
        keep = 0;
        jech = iech + iwgt * ilag;
        if (jech < 0 || jech > nech) break;
        if (hasSel && !db->isActive(jech)) break;
        db->getSampleAsSTInPlace(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! keepPair(idir, T1, T2, &dist)) continue;
        if (iwgt == 1) dist0 = dist;

        /* Evaluate the variogram */

        zz = _getIVAR(db, jech, 0);
        if (FFFF(zz)) break;
        keep = 1;
        value += zz * VARWGT[norder][iwgt];
      }
      if (keep)
      {
        value = value * value / NORWGT[norder];
        _setResult(iech, iech, nvar, ilag, 0, 0, 0, 1., dist0, value);
      }
    }
  }

  /* Scale the variogram calculations */

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);
}

/****************************************************************************/
/*!
 **  Patch the value of C(0) for covariances
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
void Vario::_patchC00(Db *db, int idir)
{
  double z1, z2, s12w, s12wzz, ww, scale, value, m1, m2, sumw;

  /* Initializations */

  if (!getFlagAsym()) return;

  /* Calculate the C00 term */

  for (int ivar = 0; ivar < db->getNLoc(ELoc::Z); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      int i = getDirAddress(idir, ivar, jvar, 0, false, 0);
      setHhByIndex(idir, i, 0.);

      m1 = m2 = s12w = s12wzz = sumw = 0.;

      /* Calculate the statistics for each variable */

      for (int iech = 0; iech < db->getNSample(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = _getIVAR(db, iech, ivar);
        z2 = _getIVAR(db, iech, jvar);
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
          _printDebug(iech, iech, ivar, jvar, i, scale, value);
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
}

/****************************************************************************/
/*!
 **  Evaluate the experimental generalized variogram on grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  norder     Order of the Generalized variogram
 **
 *****************************************************************************/
int Vario::_calculateGenOnGrid(DbGrid *db, int norder)
{
  /* Preliminary checks */

  if (getNVar() != 1)
  {
    messerr("The generalized variogram requires a single variable");
    return (1);
  }

  /* Update the global statistics */

  _getStatistics(db);

  /* Loop on the directions to evaluate */

  for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
  {
    if (_calculateGenOnGridSolution(db, idir, norder)) return 1;
  }
  return 0;
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
int Vario::_calculateOnGrid(DbGrid *db)
{
  int iadd_new = -1;
  int iatt_old = -1;
  double maille = 0.;

  /* In the case of Covariogram, add the weight set to the scale */

  if (getCalcul() == ECalcVario::COVARIOGRAM)
  {
    iatt_old = db->getUIDByLocator(ELoc::W, 0);
    iadd_new = db->addColumnsByConstant(1, 0.);
    if (iadd_new < 0) return 1;
    db->setLocatorByUID(iadd_new, ELoc::W, 0);
    maille = db->getCellSize();
    for (int iech = 0; iech < db->getNSample(); iech++)
      db->setLocVariable(ELoc::W, iech, 0, maille);
  }

  /* Evaluate the drift coefficients */

  if (_flag_UK)
  {
    if (_driftEstimateCoefficients(db)) return 1;
  }

  /* Update the global statistics */

  _getStatistics(db);

  /* Loop on the directions to evaluate */

  for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
  {
    if (_calculateOnGridSolution(db, idir)) return 1;
  }

  /* Delete the additional weight variable (optional) */

  if (getCalcul() == ECalcVario::COVARIOGRAM)
  {
    if (iadd_new > 0) db->deleteColumnByUID(iadd_new);
    if (iatt_old > 0) db->setLocatorByUID(iatt_old, ELoc::W, 0);
  }
  return 0;
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
void Vario::_getStatistics(Db *db)
{
  double z1, z2, ww;
  int nvar = db->getNLoc(ELoc::Z);

  /* Initializations */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    setMean(0., ivar);
    for (int jvar = 0; jvar < nvar; jvar++)
      setVar(0., ivar, jvar);
  }

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double s1w = 0.;
    double s1z = 0.;
    for (int iech = 0; iech < nvar; iech++)
    {
      if (!db->isActive(iech)) continue;
      ww = db->getWeight(iech);
      if (FFFF(ww) || ww < 0.) continue;
      z1 = _getIVAR(db, iech, ivar);
      s1w += ww;
      s1z += z1;
    }

    if (s1w <= 0.) continue;
    setMean(s1z / s1w, ivar);
  }

  for (int ivar = 0; ivar < db->getNLoc(ELoc::Z); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      double s12w = 0.;
      double s12wz1 = 0.;
      double s12wz2 = 0.;
      double s12wzz = 0.;

      for (int iech = 0; iech < db->getNSample(); iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;

        z1 = _getIVAR(db, iech, ivar);
        z2 = _getIVAR(db, iech, jvar);
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
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -getVar(ivar, jvar) / getVar(jvar, jvar);
        setVar(value, ivar, jvar);
        setVar(value, jvar, ivar);
      }
  }
  else if (getCalcul() == ECalcVario::TRANS2)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < ivar; jvar++)
      {
        double value = -getVar(ivar, jvar) / getVar(ivar, ivar);
        setVar(value, ivar, jvar);
        setVar(value, jvar, ivar);
      }
  }
  else if (getCalcul() == ECalcVario::BINORMAL)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        if (ivar != jvar)
          setVar(getVar(ivar, jvar) / sqrt(getVar(ivar, ivar) * getVar(jvar, jvar)), ivar, jvar);
  }
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
int Vario::_updateVerr(Db *db, int idir, Vario_Order *vorder, int verr_mode)
{
  int ifirst, ilast, iech, jech, number, nfois;
  double dist, value, g_old, diff, sumt, sumb, wgt, sval, gval;
  static double tol = EPSILON5;
  static int maxiter = 100;

  /* Initializations */

  int nlag = getNLag(idir);

  /* Loop on the lags */

  for (int ilag = 0; ilag < nlag; ilag++)
  {
    vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);
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
        setGg(idir, 0, 0, ilag, MAX(0, getGg(idir, 0, 0, ilag) - value));
        break;

      case 2:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          nfois++;
          g_old = getGg(idir, 0, 0, ilag);
          sumt = sumb = 0.;
          for (int ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = _s(db, iech, jech);
            gval = _g(db, iech, jech);
            value = sval + getGg(idir, 0, 0, ilag);
            wgt = 1. / (value * value);
            sumt += wgt * (gval - sval);
            sumb += wgt;
          }
          setGg(idir, 0, 0, ilag, sumt / sumb);
          diff = ABS(getGg(idir, 0, 0, ilag) - g_old);
        }
        setGg(idir, 0, 0, ilag, MAX(0, getGg(idir, 0, 0, ilag)));
        if (nfois == maxiter && OptDbg::query(EDbg::CONVERGE))
          message("Convergence not reached for lag %d\n", ilag + 1);
        break;

      case 3:
        nfois = 0;
        diff = 1.e30;
        while (diff > tol && nfois < maxiter)
        {
          g_old = getGg(idir, 0, 0, ilag);
          sumt = sumb = 0.;
          for (int ipair = ifirst; ipair < ilast; ipair++)
          {
            vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);
            sval = _s(db, iech, jech);
            gval = _g(db, iech, jech);
            value = sval + getGgByIndex(idir, ilag);
            wgt = getGg(idir, 0, 0, ilag) / value;
            sumt += wgt * gval;
            sumb += 1.;
          }
          setGg(idir, 0, 0, ilag, sumt / sumb);
          diff = ABS(getGg(idir, 0, 0, ilag) - g_old);
        }
        if (nfois == maxiter && OptDbg::query(EDbg::CONVERGE))
          message("Convergence not reached for lag %d\n", ilag + 1);
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
  return (0.5 * (db->getLocVariable(ELoc::V, iech, 0) +
                 db->getLocVariable(ELoc::V, jech, 0)));
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
double Vario::_g(Db *db, int iech, int jech) const
{
  double value = _getIVAR(db, iech, 0) - _getIVAR(db, jech, 0);
  return (value * value / 2.);
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
int Vario::_getRelativeSampleRank(Db *db, int iech0)
{
  int iiech = 0;
  for (int iech = 0, nech = db->getNSample(); iech < nech; iech++)
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
void Vario::_calculateFromGeometry(Db *db, int idir, Vario_Order *vorder)
{
  int iech, jech, ifirst, ilast;
  double dist;

  /* Initializations */

  int nlag = getNLag(idir);
  int nvar = getNVar();

  /* Loop on the lags */

  for (int ilag = 0; ilag < nlag; ilag++)
  {
    vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);

    /* Loop on the pairs contributing to this lag */

    for (int ipair = ifirst; ipair < ilast; ipair++)
    {
      vario_order_get_indices(vorder, ipair, &iech, &jech, &dist);

      /* Evaluate the variogram */

      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      (this->*_evaluate)(db, nvar, iech, jech, ilag, dist, true);
    }
  }

  /* Scale the variogram calculations */

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);
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
int Vario::_calculateGeneralSolution1(Db *db,
                                      int idir,
                                      const int *rindex,
                                      Vario_Order *vorder)
{
  SpaceTarget T1(getSpace(),false);
  SpaceTarget T2(getSpace(),false);
  int iech, jech, ilag, npair, ideb;

  DirParam dirparam = getDirParam(idir);
  int nech          = db->getNSample();
  int nvar = getNVar();
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
    db->getSampleAsSTInPlace(iech, T1);

    ideb = (hasDate) ? 0 : iiech + 1;
    for (int jjech = ideb; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (db->getDistance1D(iech, jech) > maxdist) break;
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsSTInPlace(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      ilag = dirparam.getLagRank(dist);
      if (IFFFF(ilag)) continue;

      /* Case of internal storage */

      if (vorder != (Vario_Order*) NULL)
      {
        vario_order_add(vorder, iech, jech, NULL, NULL, ilag, idir, dist);
      }
      else
      {

        /* Evaluate the variogram */

        IDIRLOC = idir;
        IECH1 = iech;
        IECH2 = jech;
        (this->* _evaluate)(db, nvar, iech, jech, ilag, dist, true);
      }
    }
  }

  /* Internal storage */

  if (vorder != (Vario_Order*) NULL)
  {
    vario_order_final(vorder, &npair);
  }
  else
  {

    /* Scale the variogram calculations */

    _rescale(idir);

    /* Center the covariance function */

    _centerCovariance(db, idir);

    /* Patch the central value */

    _patchC00(db, idir);
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
int Vario::_calculateGeneralSolution2(Db *db, int idir, const int *rindex)
{
 SpaceTarget T1(getSpace(),false);
 SpaceTarget T2(getSpace(),false);
 int iech, jech, i, ilag, ideb;

 /* Initializations */

 const VarioParam& varioparam = getVarioParam();
 const DirParam& dirparam     = getDirParam(idir);
 int nech                     = db->getNSample();
 int size                     = getDirSize(idir);
 int nvar                     = getNVar();
 double maxdist               = getMaximumDistance(idir);

 /* Core allocation */

 VectorDouble gg_sum(size, 0);
 VectorDouble hh_sum(size, 0);
 VectorDouble sw_sum(size, 0);

 // Local variables to speed up calculations
 bool hasSel    = db->hasLocVariable(ELoc::SEL);
 bool hasWeight = db->hasLocVariable(ELoc::W);
 bool hasDate   = varioparam.isDateUsed(db);
 double w1      = 1.;
 double dist    = 0.;

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
   db->getSampleAsSTInPlace(iech, T1);

   /* Looking for the second sample */

   ideb = (hasDate) ? 0 : iiech + 1;
   for (int jjech = ideb; jjech < nech; jjech++)
   {
     jech = rindex[jjech];
     if (db->getDistance1D(iech, jech) > maxdist) break;
     if (hasSel && !db->isActive(jech)) continue;
     if (hasWeight && FFFF(db->getWeight(jech))) continue;
     db->getSampleAsSTInPlace(jech, T2);

     // Reject the point as soon as one BiTargetChecker is not correct
     if (!keepPair(idir, T1, T2, &dist)) continue;

     /* Get the rank of the lag */

     ilag = dirparam.getLagRank(dist);
     if (IFFFF(ilag)) continue;

     /* Evaluate the variogram */

     IECH1 = iech;
     IECH2 = jech;
     (this->*_evaluate)(db, nvar, iech, jech, ilag, dist, true);
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

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);

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
int Vario::_calculateOnGridSolution(DbGrid *db, int idir)
{
  SpaceTarget T1(getSpace(), false);
  SpaceTarget T2(getSpace(), false);

  /* Initializations */

  int nech = db->getNSample();
  int nlag = getNLag(idir);
  const DirParam &dirparam = getDirParam(idir);

  // Local variables to speed up calculations
  bool hasSel    = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double dist    = 0.;
  int nvar       = getNVar();
  int ndim       = db->getNDim();

  /* Core allocation */

  VectorInt indg1(ndim);
  VectorInt indg2(ndim);

  /* Loop on the first point */

  for (int iech = 0; iech < nech; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsSTInPlace(iech, T1);
    db->rankToIndice(iech, indg1);

    for (int ilag = 1; ilag < nlag; ilag++)
    {
      for (int idim = 0; idim < db->getNDim(); idim++)
        indg2[idim] = indg1[idim] + (int) (ilag * getGrincr(idir, idim));
      int jech = db->indiceToRank(indg2);
      if (jech < 0) continue;

      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsSTInPlace(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Evaluate the variogram */

      dist = ilag * dirparam.getDPas();
      IDIRLOC = idir;
      IECH1 = iech;
      IECH2 = jech;
      (this->*_evaluate)(db, nvar, iech, jech, ilag, dist, true);
    }
  }

  /* Scale the variogram calculations */

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);

  return 0;
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
int Vario::_calculateGenOnGridSolution(DbGrid *db, int idir, int norder)
{
 SpaceTarget T1(getSpace(),false);
 SpaceTarget T2(getSpace(),false);
 int keep;
 double zz, value;

 /* Initializations */

 int nech                 = db->getNSample();
 int nlag                 = getNLag(idir);
 int ndim                 = db->getNDim();
 int nvar                 = getNVar();
 const DirParam& dirparam = getDirParam(idir);

 // Local variables to speed up calculations

 bool hasSel = db->hasLocVariable(ELoc::SEL);
 double dist = 0.;

 /* Core allocation */

 VectorInt indg1(ndim);
 VectorInt indg2(ndim);

 /* Loop on the first point */

 for (int iech = 0; iech < nech; iech++)
 {
   if (hasSel && !db->isActive(iech)) continue;
   db->getSampleAsSTInPlace(iech, T1);
   db->rankToIndice(iech, indg1);

   for (int ilag = 1; ilag < nlag; ilag++)
   {
     value = _getIVAR(db, iech, 0);
     if (FFFF(value)) break;
     dist = ilag * getDPas(idir);

     for (int iwgt = keep = 1; iwgt < NWGT[norder] && keep; iwgt++)
     {
       keep = 0;
       for (int idim = 0; idim < db->getNDim(); idim++)
         indg2[idim] = indg1[idim] + (int)(ilag * iwgt * dirparam.getGrincr(idim));

       int jech = db->indiceToRank(indg2);
       if (jech < 0) continue;
       if (hasSel && !db->isActive(jech)) continue;
       db->getSampleAsSTInPlace(jech, T2);

       // Reject the point as soon as one BiTargetChecker is not correct
       if (!keepPair(idir, T1, T2, &dist)) continue;

       /* Evaluate the variogram */

       zz = _getIVAR(db, jech, 0);
       if (FFFF(zz)) break;
       keep = 1;
       value += zz * VARWGT[norder][iwgt];
     }
     if (keep)
     {
       value = value * value / NORWGT[norder];

       _setResult(iech, iech, nvar, ilag, 0, 0, 0, 1., dist, value);
     }
   }
  }

  /* Scale the variogram calculations */

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);

  return 0;
}

/****************************************************************************/
/*!
 **  Internal function for setting a variogram value
 **
 ** \param[in]  iech1       Rank of the first sample
 ** \param[in]  iech2       Rank of the second sample
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ilag        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  orient      Orientation
 ** \param[in]  ww          Weight
 ** \param[in]  dist        Distance
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
void Vario::_setResult(int iech1,
                       int iech2,
                       int nvar,
                       int ilag,
                       int ivar,
                       int jvar,
                       int orient,
                       double ww,
                       double dist,
                       double value)
{
  DECLARE_UNUSED(iech1);
  DECLARE_UNUSED(iech2);
  DECLARE_UNUSED(nvar);
  int i = getDirAddress(IDIRLOC, ivar, jvar, ilag, false, orient, false);
  updateGgByIndex(IDIRLOC, i, ww * value, false);
  if (getCalcul() == ECalcVario::POISSON)
    updateGgByIndex(IDIRLOC, i, -getMean(ivar) / 2., false);
  updateHhByIndex(IDIRLOC, i, ww * dist, false);
  updateSwByIndex(IDIRLOC, i, ww, false);
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
void Vario::_printDebug(int iech1,
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
void Vario::_centerCovariance(Db *db, int idir)
{
  double m1, m2, sumw, z1, z2, ww;
  if (!getFlagAsym()) return;

  /* Scale the experimental variogram quantities */

  for (int ivar = 0, nvar = getNVar(); ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      /* Calculate the mean for each variable */

      m1 = m2 = sumw = 0.;
      for (int iech = 0, nech = db->getNSample(); iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        ww = db->getWeight(iech);
        if (FFFF(ww) || ww < 0.) continue;
        z1 = _getIVAR(db, iech, ivar);
        z2 = _getIVAR(db, iech, jvar);
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
        for (int i = 0, nlagtot = getNLagTotal(idir); i < nlagtot; i++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          if (getSwByIndex(idir, j) > 0)
            setGgByIndex(idir, j, getGgByIndex(idir, j) - m1 * m2);
        }
    }
}

bool Vario::_isCompatible(const Db *db) const
{
  if (db->getNDim() != getNDim() ||
      db->getNLoc(ELoc::Z) != getNVar())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d NVAR=%d", db->getNDim(),
            db->getNLoc(ELoc::Z));
    messerr("Variogram: NDIM=%d NVAR=%d", getNDim(),
            getNVar());
    return false;
  }
  return true;
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
int Vario::computeGeometry(Db *db, Vario_Order *vorder, int *npair)
{
   SpaceTarget T1(getSpace(),false);
   SpaceTarget T2(getSpace(),false);
   int iech, jech, ideb;

   /* Initializations */

   if (db == nullptr) return 1;
   const VarioParam& varioparam = getVarioParam();

   /* Preliminary checks */

   if (!_isCompatible(db)) return 1;
   if (_get_generalized_variogram_order() > 0)
   {
     messerr("This calculation does not allow generalized variogram definition");
     return 1;
  }

  /* Sort the data */
  VectorInt rindex = db->getSortArray();

  // Local variables to speed up calculations
  bool hasSel    = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate   = varioparam.isDateUsed(db);
  int nech       = db->getNSample();
  int ndir       = getNDir();
  double dist    = 0.;

  /* Loop on the directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    const DirParam &dirparam = getDirParam(idir);
    double maxdist = getMaximumDistance(idir);

    /* Loop on the first point */

    for (int iiech = 0; iiech < nech - 1; iiech++)
    {
      iech = rindex[iiech];
      if (hasSel && !db->isActive(iech)) continue;
      if (hasWeight && FFFF(db->getWeight(iech))) continue;
      db->getSampleAsSTInPlace(iech, T1);

      ideb = (hasDate) ? 0 : iiech + 1;
      for (int jjech = ideb; jjech < nech; jjech++)
      {
        jech = rindex[jjech];
        if (db->getDistance1D(iech, jech) > maxdist) break;
        if (hasSel && !db->isActive(jech)) continue;
        if (hasWeight && FFFF(db->getWeight(jech))) continue;
        db->getSampleAsSTInPlace(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! keepPair(idir, T1, T2, &dist)) continue;

        /* Get the rank of the lag */

        int ilag = dirparam.getLagRank(dist);
        if (IFFFF(ilag)) continue;

        /* Case of internal storage */

        vario_order_add(vorder, iech, jech, NULL, NULL, ilag, idir, dist);
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

  ndir = getNDir();
  if (idir0 >= 0) ndir = 1;
  for (jdir = 0; jdir < ndir; jdir++)
  {
    idir = (idir0 >= 0) ? idir0 : jdir;
    for (i = 0; i < getNLagTotal(idir); i++)
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
int Vario::_calculateVarioVectSolution(Db *db, int idir, int ncomp, const int *rindex)
{
  SpaceTarget T1(getSpace(),false);
  SpaceTarget T2(getSpace(),false);
  int iech, jech, ilag, i, icomp;
  double w1, w2, zi1, zi2, zj1, zj2, v12, v21, di1, di2, dj1, dj2;

  const DirParam &dirparam = getDirParam(idir);
  int nech = db->getNSample();
  int nvar = getNVar();
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
    db->getSampleAsSTInPlace(iech, T1);

    for (int jjech = iiech + 1; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (db->getDistance1D(iech, jech) > maxdist) break;
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsSTInPlace(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! keepPair(idir, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      ilag = dirparam.getLagRank(dist);
      if (IFFFF(ilag)) continue;

      w1 = db->getWeight(iech);
      w2 = db->getWeight(jech);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {

          /* Evaluate the variogram */

          v12 = v21 = di1 = di2 = dj1 = dj2 = 0.;
          for (icomp = 0; icomp < ncomp; icomp++)
          {
            zi1 = _getIVAR(db, iech, ivar * ncomp + icomp);
            zi2 = _getIVAR(db, iech, jvar * ncomp + icomp);
            zj1 = _getIVAR(db, jech, ivar * ncomp + icomp);
            zj2 = _getIVAR(db, jech, jvar * ncomp + icomp);
            if (FFFF(zi1) || FFFF(zi2) || FFFF(zj1) || FFFF(zj2))
            {
              v12 = v21 = TEST;
              break;
            }
            v12 += zi1 * zj2;
            v21 += zi2 * zj1;
            di1 += zi1 * zi1;
            di2 += zi2 * zi2;
            dj1 += zj1 * zj1;
            dj2 += zj2 * zj2;
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

          i = getDirAddress(idir, ivar, jvar, ilag, false, 1);
          setGgByIndex(idir, i, getGgByIndex(idir, i) + w1 * w2 * v12);
          setHhByIndex(idir, i, getHhByIndex(idir, i) + w1 * w2 * dist);
          setSwByIndex(idir, i, getSwByIndex(idir, i) + w1 * w2);

          i = getDirAddress(idir, ivar, jvar, ilag, false, -1);
          setGgByIndex(idir, i, getGgByIndex(idir, i) + w1 * w2 * v21);
          setHhByIndex(idir, i, getHhByIndex(idir, i) + w1 * w2 * dist);
          setSwByIndex(idir, i, getSwByIndex(idir, i) + w1 * w2);
        }
    }
  }

  /* Scale the variogram calculations */

  _rescale(idir);

  /* Center the covariance function */

  _centerCovariance(db, idir);

  /* Patch the central value */

  _patchC00(db, idir);

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
int Vario::computeVarioVect(Db *db, int ncomp)
{
  if (db == nullptr) return (1);

  /* Preliminary checks */

  if (! _isCompatible(db)) return 1;

  /* Update the global statistics */

  _getVarioVectStatistics(db, ncomp);

  /* Loop on the directions to evaluate */

  VectorInt rindex = db->getSortArray();
  for (int idir = 0; idir < getNDir(); idir++)
  {
    if (_calculateVarioVectSolution(db, idir, ncomp, rindex.data())) return 1;
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
void Vario::_getVarioVectStatistics(Db *db, int ncomp)
{
  double vi, vj, vij, s12ww, s12wzz, zi, zj, ww;

  /* Loop on the variables */

  int nb_neg = 0;
  for (int ivar = 0; ivar < getNVar(); ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {

      /* Loop on the samples */

      s12ww = s12wzz = 0.;
      for (int iech = 0; iech < db->getNSample(); iech++)
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
          zi = _getIVAR(db, iech, ivar * ncomp + icomp);
          zj = _getIVAR(db, iech, jvar * ncomp + icomp);
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
}

/****************************************************************************/
/*!
 **  Scale the variogram calculations
 **
 ** \param[in]  idir  Rank of the Direction
 **
 *****************************************************************************/
void Vario::_rescale(int idir)
{
  int nvar = getNVar();

  /* Scale the experimental variogram quantities */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      for (int i = 0; i < getNLagTotal(idir); i++, ecr++)
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
          if (getFlagAsym() && i < getNLag(idir))
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
        for (int i = 0; i < getNLagTotal(idir); i++, ecr++)
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
        for (int i = 0; i < getNLagTotal(idir); i++, ecr++)
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
        for (int i = 0; i < getNLagTotal(idir); i++, ecr++)
        {
          int j = getDirAddress(idir, ivar, jvar, i, true, 0);
          int j1 = getDirAddress(idir, ivar, ivar, i, true, 0);
          int j2 = getDirAddress(idir, jvar, jvar, i, true, 0);
          setGgByIndex(idir, j,
              getGgByIndex(idir, j) / sqrt(getGgByIndex(idir, j1) * getGgByIndex(idir, j2)));
        }
      }
  }
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
double Vario::_getIVAR(const Db *db, int iech, int ivar) const
{
  double zz = db->getZVariable(iech, ivar);
  if (FFFF(zz)) return (TEST);
  if (_BETA.empty()) return (zz);
  if (ivar != 0) return (TEST);
  if (_model == nullptr) return TEST;
  double drfval = _model->evalDriftVarCoef(db, iech, 0, _BETA);
  if (FFFF(drfval)) return (TEST);
  return (zz - drfval);
}

bool Vario::keepPair(int idir, SpaceTarget &T1, SpaceTarget &T2, double *dist) const
{
  for (int ipt = 0, npt = getNBiPtsPerDir(); ipt < npt; ipt++)
  {
    const ABiTargetCheck* bipts = getBipts(idir, ipt);
    if (! bipts->isOK(T1, T2)) return false;
    const BiTargetCheckGeometry* bigeom = dynamic_cast<const BiTargetCheckGeometry*>(bipts);
    if (bigeom != nullptr) *dist = bigeom->getDist();
  }
  return true;
}

/****************************************************************************/
/*!
 **  Ask for the rank of the 'vardir' structure, given direction and date
 **
 ** \return  Absolute rank (or -1 for error)
 **
 ** \param[in]  idir   Rank for the direction (starting from 0)
 ** \param[in]  idate  Rank for the Date (starting from 0)
 **
 ** \remark  An error occurs if 'idir' is negative or larger than 'ndir'
 ** \remark  or if 'idate' is negative or larger than 'ndate'
 **
 *****************************************************************************/
int Vario::getRankFromDirAndDate(int idir, int idate) const
{
  int rank = idir;
  int ndir = getNDir();
  int ndate = getNDate();
  if (idir < 0 || idir >= ndir) return (-1);
  if (ndate > 0)
  {
    if (idate < 0 || idate >= ndate) return (-1);
    rank = rank * idir + idate;
  }
  return (rank);
}

/****************************************************************************/
/*!
 **  Manage the drift removal option
 **
 ** \param[in]  db      Db structure
 **
 *****************************************************************************/
void Vario::_driftManage(Db *db)
{
  if (_model == nullptr) return;

  int nbfl = _model->getNDrift();
  int nech = db->getNSampleActiveAndDefined(0);

  _BETA.resize(nbfl,0.);
  _DRFDIAG.resize(nech, 0.);
  _DRFTAB.resetFromValue(nech, nbfl, 0.);
  _DRFXA.resetFromValue(nech, nbfl, 0.);
  _DRFGX.resetFromValue(nech, nbfl, 0.);
  _DRFXGX.resetFromValue(nbfl, nbfl, 0.);
}

/****************************************************************************/
/*!
 **  Estimate the coefficients of the global drift
 **
 ** \return Error return code
 **
 ** \param[in]  db      Db descriptor
 **
 *****************************************************************************/
int Vario::_driftEstimateCoefficients(Db *db)
{
  if (_model == nullptr) return 1;
  int iiech;
  int nbfl = _model->getNDrift();
  VectorDouble b(nbfl, 0.);
  MatrixSquareGeneral matdrf(nbfl);

  /* Calculate: t(X) %*% X */

  for (int iech = iiech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    VectorDouble drfloc = _model->evalDriftBySample(db, iech, ECalcMember::LHS);
    double zval = db->getZVariable( iech, 0);

    for (int il = 0; il < nbfl; il++)
    {
      if (FFFF(drfloc[il]))
      {
        messerr("Drift cannot be calculated: term (%d) is undefined at sample (%d)",
                il + 1, iech + 1);
        return 1;
      }
      _DRFTAB.setValue(iiech, il, drfloc[il]);
      b[il] += drfloc[il] * zval;
      for (int jl = 0; jl < nbfl; jl++)
        matdrf.setValue(il,jl, matdrf.getValue(il,jl) + drfloc[il] * drfloc[jl]);
    }
    iiech++;
  }

  /* Calculate: matdrf = (t(X) %*% X)-1 */

  if (matdrf.invert()) return 1;

  /* Calculate: _BETA = (t(X) %*% X)-1 %*% t(X) %*% Y */

  matdrf.prodMatVecInPlace(b, _BETA);

  /* Optional printout */

  if (_verbose)
  {
    message("Drift removal initial step\n");
    print_matrix("Drift Coefficients Matrix", 0, 1, nbfl, nbfl, NULL, matdrf.getValues().data());
  }

  /* Pre-process the vector X %*% (t(X) %*% X)-1 */

  _DRFXA.prodMatMatInPlace(&_DRFTAB, &matdrf);

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the a bias term between samples iech and jech
 **
 ** \param[in]  iiech Relative rank of the first sample
 ** \param[in]  jjech Relative rank of the second sample
 **
 *****************************************************************************/
double Vario::_getBias(int iiech, int jjech)
{
  int nbfl = _model->getNDrift();

  double bias0 = 0.;
  for (int il = 0; il < nbfl; il++)
    for (int jl = 0; jl < nbfl; jl++)
      bias0 += _DRFXA.getValue(iiech, il) * _DRFXGX.getValue(il, jl) * _DRFXA.getValue(jjech, jl);

  double bias1 = 0.;
  for (int il = 0; il < nbfl; il++)
    bias1 += _DRFXA.getValue(iiech, il) * _DRFGX.getValue(jjech, il);

  double bias2 = 0.;
  for (int il = 0; il < nbfl; il++)
    bias2 += _DRFGX.getValue(iiech, il) * _DRFXA.getValue(jjech, il);

  return (bias0 - (bias1 + bias2));
}

/****************************************************************************/
/*!
 **  Linear interpolation
 **
 ** \return  Interpolated value
 **
 ** \param[in]  n      Number of discretization steps
 ** \param[in]  x      Discretized X (sorted increasingly)
 ** \param[in]  y      Discretized Y
 ** \param[in]  x0     Origin
 **
 *****************************************************************************/
double Vario::_linear_interpolate(int n,
                                  const VectorDouble &x,
                                  const VectorDouble &y,
                                  double x0)
{
  if (x0 < x[0]) return (y[0]);
  if (x0 > x[n - 1]) return (y[n - 1]);
  for (int i = 1; i < n; i++)
  {
    if (x0 < x[i - 1]) continue;
    if (x0 > x[i]) continue;
    return (y[i - 1] + (y[i] - y[i - 1]) * (x0 - x[i - 1]) / (x[i] - x[i - 1]));
  }
  return (TEST);
}

/****************************************************************************/
/*!
 **  Update the experimental variogram of the completed variable starting
 **  from the experimental variogram of the truncated variable
 **  This only functions in the monovariate case
 **
 ** \param[in]  nh     Number of Hermite polynomials
 ** \param[in]  ycut   Truncation (lowest) value
 **
 ** \return Error return code
 **
 *****************************************************************************/
int Vario::transformCut(int nh, double ycut)
{
  if (getNVar() != 1)
  {
    messerr("The method 'transformCut' is available in the monovariate case only");
    return 1;
  }
  static double disc = 0.01;
  int ndisc = (int) (2. / disc + 1.);

  /* Core allocation */

  VectorDouble ro(ndisc);
  VectorDouble covyp(ndisc);
  for (int idisc = 0; idisc < ndisc; idisc++)
    ro[idisc] = disc * idisc - 1.;

  /* Calculate the first normalized Hermite polynomials for ycut */

  VectorDouble psic = hermiteCoefLower(ycut, nh);

  /* Variance */

  double variance = 0.;
  for (int ih = 1; ih < nh; ih++)
    variance += psic[ih] * psic[ih];

  for (int idisc = 0; idisc < ndisc; idisc++)
  {
    double sum = 0.;
    for (int ih = 1; ih < nh; ih++)
      sum += psic[ih] * psic[ih] * pow(ro[idisc], ih);
    covyp[idisc] = sum;
  }

  /* Loop on the directions */

  for (int idir = 0, ndir = getNDir(); idir < ndir; idir++)
  {

    /* Loop on the lags */

    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {
      double cyp = variance - getGg(idir, 0, 0, ilag);
      double cyy = _linear_interpolate(ndisc, covyp, ro, cyp);
      setGg(idir, 0, 0, ilag, MAX(0, 1. - cyy));
    }
  }

  /* Update the variance */

  setVar(1., 0, 0);

  return 0;
}

/****************************************************************************/
/*!
 **  Determine the samples used for a variogram in multilayers framework
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  seltab Number of sample definition (0, 1 or 2)
 **
 ** \param[out]  vorder Vario_Order structure
 **
 *****************************************************************************/
int Vario::computeGeometryMLayers(Db *db,
                                  VectorInt &seltab,
                                  Vario_Order *vorder) const
{
   SpaceTarget T1(getSpace(),false);
   SpaceTarget T2(getSpace(),false);
   int iiech, jjech, npair;

   /* Initializations */

   if (db == nullptr) return 1;

   // Local variables to speed up calculations
   bool hasSel = db->hasLocVariable(ELoc::SEL);
   int nech    = db->getNSample();
   int ndir    = getNDir();
   double dist = 0.;

   /* Loop on the directions */

   for (int idir = 0; idir < ndir; idir++)
   {
     const DirParam& dirparam = getDirParam(idir);

     /* Loop on the first point */

     for (int iech = iiech = 0; iech < nech; iech++)
     {
       if (hasSel && !db->isActive(iech)) continue;
       db->getSampleAsSTInPlace(iech, T1);

       if (seltab[iech] == 0) continue;
       for (int ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
       {
         for (int jech = jjech = 0; jech < nech; jech++)
         {
           if (hasSel && !db->isActive(jech)) continue;
           if (seltab[jech] == 0) continue;
           db->getSampleAsSTInPlace(jech, T2);

           for (int jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
           {

             // Reject the point as soon as one BiTargetChecker is not correct
             if (!keepPair(idir, T1, T2, &dist)) continue;

             /* Get the rank of the lag */

             int ilag = dirparam.getLagRank(dist);
             if (IFFFF(ilag)) continue;

             /* Internal storage */

             vario_order_add(vorder, iiech, jjech, &iech, &jech, ilag, idir,
                             ABS(dist));
           }
         }
       }
     }
  }
  vario_order_final(vorder, &npair);
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculates variogram values by sampling a model
 **
 ** \return  Error return code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  mode      CovCalcMode structure
 **
 *****************************************************************************/
int Vario::sampleModel(Model *model, const CovCalcMode*  mode)
{
  int ndim = getNDim();
  int ndir = getNDir();
  int nvar = model->getNVar();

  /* Core allocation */

  VectorDouble d1(ndim,0.);
  MatrixSquareGeneral covtab(nvar);

  setNVar(nvar);
  internalVariableResize();
  internalDirectionResize();

  /* Calculate the C(0) constant term */

  model->evaluateMatInPlace(nullptr, VectorDouble(), covtab, true, 1., mode);
  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++, ecr++)
      setVarIndex(ecr, covtab.getValue(ivar,jvar));

  /* Loop on the directions */

  for (int idir = 0; idir < ndir; idir++)
  {

    /* Loop on the variogram lags */

    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
    {

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          int i = getDirAddress(idir, ivar, jvar, ilag, false, 0);
          setSwByIndex(idir, i, 1.);
          setHhByIndex(idir, i, ilag * getDPas(idir));
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = getHhByIndex(idir, i) * getCodir(idir, idim);
          model->evaluateMatInPlace(nullptr, d1, covtab, true, 1., mode);
          setGgByIndex(idir, i, covtab.getValue(ivar, jvar));
        }
    }
  }
  return 0;
}

void Vario::setVariableName(int ivar, const String &variableName)
{
  if (! _isVariableValid(ivar)) return;
  _variableNames[ivar] = variableName;
}

String Vario::getVariableName(int ivar) const
{
  if (! _isVariableValid(ivar)) return String();
  return _variableNames[ivar];
}

bool Vario::isLagCorrect(int idir, int k) const
{
  double hh = getHhByIndex(idir, k);
  if (isZero(hh) || FFFF(hh)) return false;
  double sw = getSwByIndex(idir, k);
  if (isZero(sw) || FFFF(sw)) return false;
  double gg = getGgByIndex(idir, k);
  return !FFFF(gg);
}

double Vario::getC00(int idir, int ivar, int jvar) const
{
  int iad0           = getDirAddress(idir, ivar, jvar, 0, false, 0);
  int iad            = iad0;
  double c00         = getSwByIndex(idir, iad);
  if (!isZero(c00) || getSwByIndex(idir, iad) > 0) return c00;

  for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
  {
    iad = getDirAddress(idir, ivar, jvar, ilag, false, 1);
    if (!isZero(getGgByIndex(idir, iad)))
      return getGgByIndex(idir, iad);
    iad = getDirAddress(idir, ivar, jvar, ilag, false, -1);
    if (!isZero(getGgByIndex(idir, iad)))
      return getGgByIndex(idir, iad);
  }
  iad = iad0;
  return (getGgByIndex(idir, iad));
}

VectorDouble Vario::computeWeightPerDirection() const
{
  int ndir = getNDir();
  int nvar = getNVar();
  int nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble count(ndir);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    for (int ilag = 0, nlag = getNLag(idir); ilag < nlag; ilag++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * getNLagTotal(idir);
        if (getFlagAsym())
        {
          int iad   = shift + getNLag(idir) + ilag + 1;
          int jad   = shift + getNLag(idir) - ilag - 1;
          double n1 = getSwByIndex(idir, iad);
          double n2 = getSwByIndex(idir, jad);
          if (isLagCorrect(idir, iad)) count[idir] += n1;
          if (isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ilag;
          double nn = getSwByIndex(idir, iad);
          if (isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }
  return count;
}

int Vario::getTotalLagsPerDirection() const
{
  int npatot = 0;
  int ndir   = getNDir();
  for (int idir = 0; idir < ndir; idir++)
    npatot += getNLagTotal(idir);
  return npatot;
}

VectorDouble Vario::computeWeightsFromVario(int wmode) const
{
  int ndir           = getNDir();
  int nvar           = getNVar();
  int npadir         = getTotalLagsPerDirection();
  VectorDouble count = computeWeightPerDirection();
  int nvs2           = nvar * (nvar + 1) / 2;
  VectorDouble wt(npadir * nvs2, 0.);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    int nlag    = getNLag(idir);
    for (int ilag = 0; ilag < nlag; ilag++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * getNLagTotal(idir);
        if (getFlagAsym())
        {
          int iad   = shift + getNLag(idir) + ilag + 1;
          int jad   = shift + getNLag(idir) - ilag - 1;
          double n1 = getSwByIndex(idir, iad);
          double n2 = getSwByIndex(idir, jad);
          if (isLagCorrect(idir, iad)) count[idir] += n1;
          if (isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ilag;
          double nn = getSwByIndex(idir, iad);
          if (isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }

  int ipadir = 0;
  switch (wmode)
  {
    case 1:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int nlag = getNLag(idir);
        for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * getNLagTotal(idir);
            if (getFlagAsym())
            {
              int iad = shift + getNLag(idir) + ilag + 1;
              int jad = shift + getNLag(idir) - ilag - 1;
              if (isLagCorrect(idir, iad) &&
                  isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = count[idir];
            }
            else
            {
              int iad = shift + ilag;
              if (isLagCorrect(idir, iad))
                WT(ijvar, ipadir) = count[idir];
            }
          }
        }
      }
      break;

    case 2:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int nlag = getNLag(idir);
        for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * getNLagTotal(idir);
            if (getFlagAsym())
            {
              int iad = shift + getNLag(idir) + ilag + 1;
              int jad = shift + getNLag(idir) - ilag - 1;
              if (isLagCorrect(idir, iad) ||
                  isLagCorrect(idir, jad))
                continue;
              double n1 = getSwByIndex(idir, iad);
              double n2 = getSwByIndex(idir, jad);
              double d1 = ABS(getHhByIndex(idir, iad));
              double d2 = ABS(getHhByIndex(idir, jad));
              if (d1 > 0 && d2 > 0)
                WT(ijvar, ipadir) =
                  sqrt((n1 + n2) * (n1 + n2) / (n1 * d1 + n2 * d2) / 2.);
            }
            else
            {
              int iad = shift + ilag;
              if (!isLagCorrect(idir, iad)) continue;
              double nn = getSwByIndex(idir, iad);
              double dd = ABS(getHhByIndex(idir, iad));
              if (dd > 0) WT(ijvar, ipadir) = nn / dd;
            }
          }
        }
      }
      break;

    case 3:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int nlag = getNLag(idir);
        for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * getNLagTotal(idir);
            if (getFlagAsym())
            {
              int iad = shift + getNLag(idir) + ilag + 1;
              int jad = shift + getNLag(idir) - ilag - 1;
              if (isLagCorrect(idir, iad) &&
                  isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1. / getNLag(idir);
            }
            else
            {
              int iad = shift + ilag;
              if (isLagCorrect(idir, iad))
                WT(ijvar, ipadir) = 1. / getNLag(idir);
            }
          }
        }
      }
      break;

    default:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int nlag = getNLag(idir);
        for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * getNLagTotal(idir);
            if (getFlagAsym())
            {
              int iad = shift + getNLag(idir) + ilag + 1;
              int jad = shift + getNLag(idir) - ilag - 1;
              if (isLagCorrect(idir, iad) &&
                  isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1.;
            }
            else
            {
              int iad = shift + ilag;
              if (isLagCorrect(idir, iad)) WT(ijvar, ipadir) = 1.;
            }
          }
        }
      }
      break;
  }

  /* Scaling by direction and by variable */

  for (int ijvar = 0; ijvar < nvs2; ijvar++)
  {
    ipadir = 0;
    for (int idir = 0; idir < ndir; idir++)
    {
      double total = 0.;
      int nlag     = getNLag(idir);
      for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          total += WT(ijvar, ipadir);
      }
      if (isZero(total)) continue;
      ipadir -= getNLag(idir);
      for (int ilag = 0, nlag = getNLag(idir); ilag < nlag;
           ilag++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          WT(ijvar, ipadir) /= total;
      }
    }
  }

  /* Scaling by variable variances */

  int ijvar0 = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar0++)
    {
      double ratio =
        (getVar(ivar, jvar) > 0 && getVar(jvar, ivar) > 0)
          ? sqrt(getVar(ivar, jvar) * getVar(jvar, ivar))
          : 1.;
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int nlag = getNLag(idir);
        for (int ilag = 0; ilag < nlag; ilag++, ipadir++)
          if (!FFFF(WT(ijvar0, ipadir))) WT(ijvar0, ipadir) /= ratio;
      }
    }

  // Ultimate check

  double total = VH::cumul(wt);
  if (ABS(total) <= 0.)
  {
    messerr("The sum of the weight is 0. This must be an error");
    wt.clear();
  }
  return wt;
}
