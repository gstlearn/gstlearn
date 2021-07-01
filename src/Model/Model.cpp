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
#include "Model/NoStatArray.hpp"

Model::Model(const CovContext& ctxt, bool flagGradient, bool flagLinked)
    : _flagGradient(flagGradient),
      _flagLinked(flagLinked),
      _covaList(nullptr),
      _driftList(nullptr),
      _modTrans(),
      _noStat(),
      _ctxt(ctxt),
      generic_cov_function(nullptr)
{
  // Create the multi-Covariance structure adequately
  if (flagGradient)
    _covaList = new CovLMGradient(_ctxt.getSpace());
  else
    _covaList = new CovLMC(_ctxt.getSpace());
  _driftList = new ADriftList(flagLinked);
}

Model::Model(const Db *db, bool flagGradient, bool flagLinked)
    : _flagGradient(flagGradient),
      _flagLinked(flagLinked),
      _covaList(nullptr),
      _driftList(nullptr),
      _modTrans(),
      _noStat(),
      _ctxt(),
      generic_cov_function(nullptr)
{
  _ctxt = CovContext(db);

  // Create the multi-Covariance structure adequately
  if (flagGradient)
    _covaList = new CovLMGradient(_ctxt.getSpace());
  else
    _covaList = new CovLMC(_ctxt.getSpace());
  _driftList = new ADriftList(flagLinked);
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
  _covaList->delAllCov();
  _driftList->delAllDrift();
  delete _covaList;
  delete _driftList;
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

  sstr << _noStat.toString();

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
VectorDouble Model::sampleModel(double hmax,
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
 * @param types   Vector of ENUM_COVS (treated as int due to compiler problem)
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

  // Clean out possible covariances in tRG  he existing model

  delAllCovas();

  // Add the relevant covariances

  CovContext ctxt = CovContext(vario->getVariableNumber());
  for (int is = 0; is < (int) types.size(); is++)
  {
    ENUM_COVS covtype = ENUM_COVS(types[is]);
    CovAniso cov = CovAniso(covtype,ctxt);
    addCova(&cov);
  }

  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}
