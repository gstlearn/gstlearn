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
#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/Vector.hpp"
#include "Space/SpaceRN.hpp"
#include "MatrixC/MatrixCSSym.hpp"
#include "Basic/AException.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Model/NoStatArray.hpp"

Model::Model(const CovContext& ctxt, bool flagGradient, bool flagLinked)
    : AStringable(),
      ASerializable(),
      _flagGradient(flagGradient),
      _flagLinked(flagLinked),
      _covaList(nullptr),
      _driftList(nullptr),
      _modTrans(),
      _noStat(),
      _ctxt(ctxt),
      generic_cov_function(nullptr)
{
  _create(flagGradient, flagLinked);
}

Model::Model(const Db *db, bool flagGradient, bool flagLinked)
    : AStringable(),
      ASerializable(),
      _flagGradient(flagGradient),
      _flagLinked(flagLinked),
      _covaList(nullptr),
      _driftList(nullptr),
      _modTrans(),
      _noStat(),
      _ctxt(),
      generic_cov_function(nullptr)
{
  _ctxt = CovContext(db);
  _create(flagGradient, flagLinked);
}

Model::Model(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      _flagGradient(false),
      _flagLinked(false),
      _covaList(nullptr),
      _driftList(nullptr),
      _modTrans(),
      _noStat(),
      _ctxt(),
      generic_cov_function(nullptr)
{
  (void) deSerialize(neutralFileName, verbose);
}

Model::Model(const Model &m)
    : _flagGradient(m._flagGradient),
      _flagLinked(m._flagLinked),
      _covaList(dynamic_cast<ACovAnisoList*>(m._covaList->clone())),
      _driftList(dynamic_cast<ADriftList*>(m._driftList->clone())),
      _modTrans(m._modTrans),
      _noStat(m._noStat),
      _ctxt(m._ctxt),
      generic_cov_function(m.generic_cov_function)
{
}

Model& Model::operator= (const Model &m)
{
  if (this != &m)
  {
    _flagGradient = m._flagGradient;
    _flagLinked   = m._flagLinked;
    _covaList     = dynamic_cast<ACovAnisoList*>(m._covaList->clone());
    _driftList    = dynamic_cast<ADriftList*>(m._driftList->clone());
    _modTrans     = m._modTrans;
    _noStat       = m._noStat;
    _ctxt         = m._ctxt;
  }
  return (*this);
}

Model::~Model()
{
  _destroy();
}

String Model::toString(int level) const
{
  std::stringstream sstr;
  int ncov   = _covaList->getCovNumber();
  int ndrift = _driftList->getDriftNumber();

  sstr << toTitle(0,"Model characteristics");
  if (_flagGradient)
    sstr << "(Specific for Handling Gradient)" << std::endl;
  sstr << "Space dimension              = " << getDimensionNumber() << std::endl;
  sstr << "Number of variable(s)        = " << getVariableNumber() << std::endl;
  sstr << "Number of basic structure(s) = " << ncov << std::endl;
  sstr << "Number of drift function(s)  = " << ndrift << std::endl;
  sstr << "Number of drift equation(s)  = " << getDriftEquationNumber() << std::endl;

  /* Covariance part */

  if (ncov > 0)
  {
    sstr << toTitle(1, "Covariance Part");
    sstr << _covaList->toString();
  }

  /* Drift part */

  if (ndrift > 0)
  {
    sstr << toTitle(1, "Drift Part");
    sstr << _driftList->toString();
  }

  // Non-stationary parameters

  if (isNoStat()) sstr << _noStat->toString();

  // Model Transformation Option

  sstr << _modTrans.toString();

  return sstr.str();
}

void Model::delCova(int rank)
{
  _covaList->delCov(rank);
}

void Model::delAllCovas()
{
  _covaList->delAllCov();
}

void Model::addCova(const CovAniso* cov)
{
  _covaList->addCov(cov);
  model_setup(this);
}

void Model::addDrift(const ADriftElem* drift)
{
  _driftList->addDrift(drift);
}

void Model::addDrift(const VectorString& driftSymbols)
{
  ENUM_DRIFTS type;
  int rank;

  for (int i = 0; i < (int) driftSymbols.size(); i++)
  {
    DriftFactory::identifyDrift(driftSymbols[i], &type, &rank,_ctxt);
    ADriftElem* drift = DriftFactory::createDriftFunc(type, _ctxt);
    drift->setRankFex(rank);
    addDrift(drift);
  }
}

void Model::delDrift(int rank)
{
  _driftList->delDrift(rank);
}

void Model::delAllDrifts()
{
  _driftList->delAllDrift();
}

int Model::hasExternalCov() const
{
  for (int icov=0; icov<(int) _covaList->getCovNumber(); icov++)
  {
    if (_covaList->getType(icov) == COV_FUNCTION) return 1;
  }
  return 0;
}

int Model::addNoStat(ANoStat* anostat)
{
  if (getDimensionNumber() > 3)
  {
    messerr("Non stationary model is restricted to Space Dimension <= 3");
    return 1;
  }

  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    int icov = getNoStatElemIcov(ipar);
    ENUM_CONS type = getNoStatElemType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the Model definition

    if (icov < 0 || icov >= getCovaNumber())
    {
      messerr(
          "Invalid Covariance rank (%d) for the Non-Stationary Parameter (%d)",
          icov, ipar);
      return 1;
    }
    if (type == CONS_PARAM)
    {
      messerr(
          "The current methodology does not handle constraint on third parameter");
      return 1;
    }
  }
  _noStat = anostat;
  return 0;
}

int Model::isNoStat() const
{
  return _noStat != nullptr;
}

int Model::getNoStatElemNumber() const
{
  if (! isNoStat()) return 0;
  return _noStat->getNoStatElemNumber();
}

void Model::addNoStatElem(int igrf, int icov, ENUM_CONS type, int iv1, int iv2)
{
  if (! isNoStat()) return;
  _noStat->addNoStatElem(igrf, icov, type, iv1, iv2);
}

void Model::addNoStatElems(const VectorString& codes)
{
  if (! isNoStat()) return;
  _noStat->addNoStatElems(codes);
}

ConsItem Model::getConsItem(int ipar) const
{
  if (! isNoStat())
    my_throw("Nostat is not defined and cannot be returned");
  return _noStat->getItems(ipar);
}

int Model::getNoStatElemIcov(int ipar)
{
  if (! isNoStat())
    my_throw("Nostat is not defined");
  return _noStat->getICov(ipar);
}

ENUM_CONS Model::getNoStatElemType(int ipar)
{
  if (! isNoStat())
    my_throw("Nostat is not defined");
  return _noStat->getType(ipar);
}

double Model::evaluateDrift(const Db* db, int iech, int il, int member) const
{
  if (member != MEMBER_LHS && isDriftFiltered(il))
    return 0.;
  else
  {
    ADriftElem* drift = _driftList->getDrift(il);
    if (drift != nullptr) return drift->eval(db, iech);
  }
  return TEST;
}

/**
 * Sample a Model for given variable(s) and given direction
 * @param hmax   Maximum distance to be sampled
 * @param nh     Number of discretization steps
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param codir  Vector of direction coefficients
 * @param nostd  0 standard; +-1 corr. envelop; ITEST normalized
 *
 * @return
 */
VectorDouble Model::sample(double hmax,
                           int nh,
                           int ivar,
                           int jvar,
                           VectorDouble codir,
                           int nostd)
{
  VectorDouble hh, gg;

  if (ivar < 0 || ivar >= getVariableNumber()) return gg;
  if (jvar < 0 || jvar >= getVariableNumber()) return gg;
  int ndim = getDimensionNumber();
  if (codir.empty())
  {
    codir.resize(ndim);
    (void) ut_angles_to_codir(ndim, 1, VectorDouble(), codir);
  }
  hh.resize(nh);
  gg.resize(nh);

  for (int i = 0; i < nh; i++) hh[i] = hmax * i / nh;

  model_evaluate(this, ivar, jvar, -1, 0, 0, 0, nostd, 0, 0,
                 nh, codir, hh.data(), gg.data());
  return gg;
}

/**
 * Automatic Fitting procedure
 * @param vario   Experimental variogram to be fitted
 * @param types   Vector of ENUM_COVS (treated as 'int' due to compiler problem)
 * @param verbose Verbose option
 * @param mauto   Special parameters for Automatic fitting procedure
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 * @return
 */
int Model::fit(Vario *vario,
               const std::vector<int>& types,
               bool verbose,
               Option_AutoFit mauto,
               const Constraints& constraints,
               Option_VarioFit optvar)
{
  if (vario == (Vario *) NULL) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  _ctxt = CovContext(vario);
  for (int is = 0; is < (int) types.size(); is++)
  {
    ENUM_COVS covtype = ENUM_COVS(types[is]);
    CovAniso cov = CovAniso(covtype,_ctxt);
    addCova(&cov);
  }

  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}

int Model::fit(Vario *vario,
               const std::vector<ENUM_COVS>& types,
               bool verbose,
               Option_AutoFit mauto,
               const Constraints& constraints,
               Option_VarioFit optvar)
{
  if (vario == (Vario *) NULL) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  _ctxt = CovContext(vario);
  for (int is = 0; is < (int) types.size(); is++)
  {
    CovAniso cov = CovAniso(types[is],_ctxt);
    addCova(&cov);
  }
  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}

int Model::deSerialize(const String& filename, bool verbose)
{
  double field, range, value, param, radius;
  int ndim, nvar, ncova, nbfl, type, flag_aniso, flag_rotation;
  VectorDouble aniso_ranges, aniso_rotmat;

  // Open the Neutral File

  if (_fileOpen(filename, "Model", "r", verbose)) return 1;

  // Delete previous Model contents (if any)
  _destroy();

  /* Create the Model structure */

  if (_recordRead("Space Dimension", "%d", &ndim)) return 1;
  if (_recordRead("Number of Variables", "%d", &nvar)) return 1;
  if (_recordRead("Field dimension", "%lf", &field)) return 1;
  if (_recordRead("Radius for Model", "%lf", &radius)) return 1;
  if (_recordRead("Number of Basic Structures", "%d", &ncova)) return 1;
  if (_recordRead("Number of Basic Drift Functions", "%d", &nbfl)) return 1;

  if (ndim <= 0 || (unsigned int)(ndim) != ASpaceObject::getGlobalSpace()->getNDim())
  {
    // TODO : Impose that global space is the same than the one deserialized
    messerr("Wrong space dimension in %s", filename);
    _destroy();
    return 1;
  }
  _ctxt = CovContext(nvar, 2, field);
  _ctxt.setBallRadius(radius);
  _create(false, false);

  /* Reading the covariance part */

  for (int icova = 0; icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    if (_recordRead("Covariance Type", "%d", &type)) return 1;
    if (_recordRead("Isotropic Range", "%lf", &range)) return 1;
    if (_recordRead("Model third Parameter", "%lf", &param)) return 1;
    if (_recordRead("Flag for Anisotropy", "%d", &flag_aniso)) return 1;
    if (flag_aniso)
    {
      aniso_ranges.resize(ndim);
      // In fact, the file contains the anisotropy coefficients
      // After reading, we must turn them into anisotropic ranges
      for (int idim = 0; idim < ndim; idim++)
        if (_recordRead("Anisotropy coefficient", "%lf",
                           &aniso_ranges[idim])) return 1;
      for (int idim = 0; idim < ndim; idim++) aniso_ranges[idim] *= range;

      if (_recordRead("Flag for Anisotropy Rotation", "%d", &flag_rotation))  return 1;
      if (flag_rotation)
      {
        // Warning: the storage in the File is performed by Column
        // whereas the internal storage (Cova) is by column
        aniso_rotmat.resize(ndim * ndim);
        int lec = 0;
        for (int idim = 0; idim < ndim; idim++)
          for (int jdim = 0; jdim < ndim; jdim++)
            if (_recordRead("Anisotropy Rotation Matrix", "%lf",
                               &aniso_rotmat[lec++])) return 1;
      }
    }

    if (isFlagGradient())
    {
      CovGradientNumerical covgrad((ENUM_COVS) type, getContext());
      covgrad.setParam(param);
      if (flag_aniso)
      {
        covgrad.setRanges(aniso_ranges);
        if (flag_rotation) covgrad.setAnisoRotation(aniso_rotmat);
      }
      else
        covgrad.setRange(range);
      addCova(&covgrad);
    }
    else
    {
      CovAniso cova((ENUM_COVS) type, getContext());
      cova.setParam(param);
      if (flag_aniso)
      {
        cova.setRanges(aniso_ranges);
        if (flag_rotation) cova.setAnisoRotation(aniso_rotmat);
      }
      else
        cova.setRange(range);
      addCova(&cova);
    }
  }

  /* Reading the drift part */

  for (int ibfl = 0; ibfl < nbfl; ibfl++)
  {
    if (_recordRead("Drift Function", "%d", &type)) return 1;
    ADriftElem *drift;
    drift = DriftFactory::createDriftFunc((ENUM_DRIFTS) type, getContext());
    drift->setRankFex(0);
    addDrift(drift);
  }

  /* Reading the matrix of means (only if nbfl <= 0) */

  if (nbfl <= 0)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double mean;
      if (_recordRead("Mean of Variable", "%lf", &mean)) return 1;
      setMean(ivar,mean);
    }

  /* Reading the matrices of sills (optional) */

  for (int icova = 0; icova < ncova; icova++)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        if (_recordRead("Matrix of Sills", "%lf", &value)) continue;
        setSill(icova,ivar,jvar,value);
      }
  }

  /* Reading the variance-covariance at the origin (optional) */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      if (_recordRead("Variance-covariance at Origin", "%lf", &value)) continue;
      setCovar0(ivar, jvar, value);
    }

  // Set the default function for calculations
  generic_cov_function = model_calcul_cov_direct;

      // Close the file
  _fileClose(verbose);

  return 0;
}

int Model::serialize(const String& filename, bool verbose) const
{
  if (_fileOpen(filename, "Model", "w", verbose)) return 1;

  /* Write the Model structure */

  _recordWrite("%d",  getDimensionNumber());
  _recordWrite("%d",  getVariableNumber());
  _recordWrite("%lf", getField());
  _recordWrite("%lf", getContext().getBallRadius());
  _recordWrite("#", "General parameters");
  _recordWrite("%d", getCovaNumber());
  _recordWrite("#", "Number of basic covariance terms");
  _recordWrite("%d", getDriftNumber());
  _recordWrite("#", "Number of drift terms");

  /* Writing the covariance part */

  for (int icova = 0; icova < getCovaNumber(); icova++)
  {
    const CovAniso* cova = getCova(icova);
    _recordWrite("%d",  cova->getType());
    _recordWrite("%lf", cova->getRange());
    _recordWrite("%lf", cova->getParam());
    _recordWrite("#", "Covariance characteristics");

    // Writing the Anisotropy information

    _recordWrite("%d", cova->getFlagAniso());
    _recordWrite("#", "Anisotropy Flag");

    if (! cova->getFlagAniso()) continue;
    for (int idim = 0; idim < getDimensionNumber(); idim++)
      _recordWrite("%lf", cova->getAnisoCoeffs(idim));
    _recordWrite("#", "Anisotropy Coefficients");
    _recordWrite("%d", cova->getFlagRotation());
    _recordWrite("#", "Anisotropy Rotation Flag");
    if (! cova->getFlagRotation()) continue;
    // Storing the rotation matrix by Column (compatibility)
    for (int idim = 0; idim < getDimensionNumber(); idim++)
      for (int jdim = 0; jdim < getDimensionNumber(); jdim++)
        _recordWrite("%lf", cova->getAnisoRotMat(jdim,idim));
    _recordWrite("#", "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (int ibfl = 0; ibfl < getDriftNumber(); ibfl++)
  {
    const ADriftElem* drift = getDrift(ibfl);
    _recordWrite("%d", drift->getType());
    _recordWrite("#", "Drift characteristics");
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (getDriftNumber() <= 0)
    for (int ivar = 0; ivar < getVariableNumber(); ivar++)
  {
    _recordWrite("%lf", getContext().getMean(ivar));
    _recordWrite("#", "Mean of Variables");
  }

  /* Writing the matrices of sills (optional) */

  for (int icova = 0; icova < getCovaNumber(); icova++)
  {
    for (int ivar = 0; ivar < getVariableNumber(); ivar++)
      for (int jvar = 0; jvar < getVariableNumber(); jvar++)
        _recordWrite("%lf",getSill(icova,ivar,jvar));
    _recordWrite("#", "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    for (int jvar = 0; jvar < getVariableNumber(); jvar++)
      _recordWrite("%lf", getContext().getCovar0(ivar, jvar));
  _recordWrite("#", "Var-Covar at origin");

  /* Close the file */

  _fileClose(verbose);

  return 0;
}

void Model::_create(bool flagGradient, bool flagLinked)
{
  if (flagGradient)
    _covaList = new CovLMGradient(_ctxt.getSpace());
  else
    _covaList = new CovLMC(_ctxt.getSpace());

  _driftList = new ADriftList(flagLinked);

  // Default function used for calculations
  generic_cov_function = model_calcul_cov_direct;
}

void Model::_destroy()
{
  if (_covaList != nullptr)  _covaList->delAllCov();
  if (_driftList != nullptr) _driftList->delAllDrift();
  delete _covaList;
  delete _driftList;
}
