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
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Enum/ECov.hpp"
#include "Enum/EModelProperty.hpp"

#include "Anamorphosis/AnamHermite.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/ANoStat.hpp"
#include "Model/CovInternal.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovGradientFunctional.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/ADrift.hpp"

#include "Db/Db.hpp"

#include <math.h>

Model::Model(const CovContext &ctxt)
    : AStringable(),
      ASerializable(),
      _cova(nullptr),
      _driftList(nullptr),
      _ctxt(ctxt)
{
  _create();
}

Model::Model(int nvar, int ndim)
    : AStringable(),
      ASerializable(),
      _cova(nullptr),
      _driftList(nullptr),
      _ctxt()
{
  SpaceRN space = SpaceRN(ndim);
  _ctxt = CovContext(nvar, &space);
  _create();
}

Model::Model(const Model &m)
    : AStringable(m),
      ASerializable(m),
      _cova(nullptr),
      _driftList(nullptr),
      _ctxt(m._ctxt)
{
  ACovAnisoList* mcovalist = dynamic_cast<ACovAnisoList*>(m._cova);
  if (mcovalist != nullptr)
    _cova = dynamic_cast<ACovAnisoList*>(mcovalist->clone());
  if (m._driftList != nullptr)
    _driftList = m._driftList->clone();
}

Model& Model::operator=(const Model &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    ACovAnisoList* mcovalist = dynamic_cast<ACovAnisoList*>(m._cova);
    if (mcovalist != nullptr)
      _cova = dynamic_cast<ACovAnisoList*>(mcovalist->clone());
    if (m._driftList != nullptr)
      _driftList = m._driftList->clone();
    _ctxt = m._ctxt;
  }
  return (*this);
}

Model::~Model()
{
  _clear();
}

int Model::resetFromDb(const Db *db)
{
  int ndim = db->getNDim();
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar <= 0) nvar = 1;
  SpaceRN space = SpaceRN(ndim);
  _ctxt = CovContext(nvar, &space);
  _create();
  return 0;
}

Model* Model::create(const CovContext& ctxt)
{
  return new Model(ctxt);
}

Model* Model::createFromEnvironment(int nvar, int ndim)
{
  return new Model(nvar, ndim);
}

Model* Model::createNugget(int nvar, int ndim, double sill)
{
  Model* model = new Model(nvar, ndim);
  model->addCovFromParam(ECov::NUGGET, 0., sill);
  return model;
}

Model* Model::createFromParam(const ECov& type,
                              double range,
                              double sill,
                              double param,
                              const VectorDouble& ranges,
                              const VectorDouble& sills,
                              const VectorDouble& angles,
                              const ASpace* space,
                              bool flagRange)
{
  int nvar = 1;
  if (! sills.empty())
    nvar = (int)  sqrt(sills.size());

  // TODO: Improve this tedious manipulation
  ASpace* spaceloc = nullptr;
  if (space != nullptr) spaceloc = dynamic_cast<ASpace*>(space->clone());

  if (! ranges.empty())
  {
    delete spaceloc;
    spaceloc = new SpaceRN((int) ranges.size());
  }

  CovContext ctxt = CovContext(nvar,spaceloc);
  Model* model = new Model(ctxt);
  model->addCovFromParam(type, range, sill, param, ranges, sills, angles, flagRange);

  delete spaceloc;
  return model;
}

Model* Model::createFromDb(const Db* db)
{
  Model* model = new Model();
  if (model->resetFromDb(db))
  {
    messerr("Problem when creating Model from Db");
    delete model;
    return nullptr;
  }
  return model;
}

Model* Model::createFromNF(const String &neutralFilename, bool verbose)
{
  Model* model = nullptr;
  std::ifstream is;
  model = new Model();
  bool success = false;
  if (model->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = model->deserialize(is, verbose);
  }

  if (! success)
  {
    delete model;
    model = nullptr;
  }
  return model;
}

String Model::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ncov   = getCovaNumber();
  int ndrift = getDriftNumber();
  if (ncov <= 0 && ndrift <= 0) return sstr.str();

  sstr << toTitle(0, "Model characteristics");
  if (isFlagGradient()) sstr << "(Specific for Handling Gradient)" << std::endl;
  sstr << "Space dimension              = " << getDimensionNumber()
       << std::endl;
  sstr << "Number of variable(s)        = " << getVariableNumber() << std::endl;
  sstr << "Number of basic structure(s) = " << ncov << std::endl;
  sstr << "Number of drift function(s)  = " << ndrift << std::endl;
  sstr << "Number of drift equation(s)  = " << getDriftEquationNumber() << std::endl;

  /* Covariance part */

  if (ncov > 0)
  {
    sstr << toTitle(1, "Covariance Part");
    sstr << _cova->toString();
  }

  /* Drift part */

  if (ndrift > 0)
  {
    sstr << toTitle(1, "Drift Part");
    sstr << _driftList->toString();

    if (isFlagLinked())
      sstr << "Drifts are linked" << std::endl;
  }

  /* Mean Part */

  if (getDriftNumber() <= 0)
  {
    sstr << toVector("Known Mean(s)", getMeans());
    // TODO: could be added but changes all non-regression files
//    sstr << "(Note: Simple Kriging will be used)" << std::endl;
  }

  return sstr.str();
}

void Model::delCova(int icov)
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return;
  covalist->delCov(icov);
}

void Model::delAllCovas()
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return;
  covalist->delAllCov();
}

/**
 * Add a list of Covariances. This operation cleans any previously stored covariance
 * @param covalist List of Covariances to be added
 */
void Model::setCovList(const ACovAnisoList* covalist)
{
  if (covalist == nullptr) return;
  delete _cova;
  _cova = dynamic_cast<ACov*>(covalist->clone());
}

void Model::addCov(const CovAniso *cov)
{
  if (cov == nullptr) return;
  if (! cov->getContext().isEqual(_ctxt))
  {
    messerr("Error: Covariance should share the same Context as 'Model'");
    messerr("Operation is cancelled");
    return;
  }
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return;
  covalist->addCov(cov);
}

void Model::addCovFromParam(const ECov& type,
                            double range,
                            double sill,
                            double param,
                            const VectorDouble& ranges,
                            const VectorDouble& sills,
                            const VectorDouble& angles,
                            bool flagRange)
{
  // Check consistency with parameters of the model

  int ndim = getDimensionNumber();
  if (! ranges.empty())
  {
    if (ndim > 0 && (int) ranges.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)",(int) ranges.size());
      messerr("and the Space dimension stored in the Model (%d)",ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = (int) ranges.size();
  }
  if (! angles.empty())
  {
    if (ndim > 0 && (int) angles.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)",(int) angles.size());
      messerr("and the Space dimension stored in the Model (%d)",ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = (int) angles.size();
  }
  int nvar = getVariableNumber();
  if (! sills.empty())
  {
    if (nvar > 0 && (int) sills.size() != nvar * nvar)
    {
      messerr("Mismatch between the size of 'sills' (%d)",(int) sills.size());
      messerr("and the Number of variables stored in the Model (%d)",nvar);
      messerr("Operation is cancelled");
      return;
    }
    nvar = (int) sqrt((double) sills.size());
  }

  // Define the covariance

  SpaceRN space = SpaceRN(ndim);
  _ctxt = CovContext(nvar, &space);
  CovAniso cov(type, _ctxt);

  double parmax = cov.getParMax();
  if (param > parmax) param = parmax;
  cov.setParam(param);
  if (! ranges.empty())
  {
    if (flagRange)
      cov.setRanges(ranges);
    else
      cov.setScales(ranges);
  }
  else
  {
    if (flagRange)
      cov.setRangeIsotropic(range);
    else
      cov.setScale(range);
  }
  if (! sills.empty())
    cov.setSill(sills);
  else
    cov.setSill(sill);

  if (! angles.empty())
    cov.setAnisoAngles(angles);
  addCov(&cov);
}

/**
 * Add a list of Drifts. This operation cleans any previously stored drift function
 * @param driftlist List of Drifts to be added
 *
 * @remark This method deletes any pre-existing drift functions
 */
void Model::setDriftList(const DriftList* driftlist)
{
  if (driftlist == nullptr) return;
  delete _driftList;
  _driftList = driftlist->clone();

  // Check that the DriftList has the same type of CovContext as the Model
  _driftList->copyCovContext(_ctxt);
}

/**
 * Define the list of drift functions for:
 * - a given degree of the IRF
 * - a given number of external drifts
 * @param order Order of the IRF
 * @param nfex  Number of External Drifts
 *
 * @remark This method deletes any pre-existing drift functions and replaces them by the new definition
 * @remark This replacement is performed accounting for information stored in 'model', such as:
 * - the space dimension
 * - the number of variables
 */
void Model::setDriftIRF(int order, int nfex)
{
  delete _driftList;
  _driftList = DriftFactory::createDriftListFromIRF(order, nfex, _ctxt);
}

void Model::setFlagLinked(bool flagLinked)
{
  if (_driftList == nullptr) return;
  _driftList->setFlagLinked(flagLinked);
}

void Model::addDrift(const ADrift *drift)
{
  if (drift == nullptr) return;
  if (_driftList == nullptr) _driftList = new DriftList(_ctxt);
  ADrift* drift_loc = dynamic_cast<ADrift*>(drift->clone());
  _driftList->addDrift(drift_loc);

  // Check that the DriftList has the same type of CovContext as the Model
  _driftList->copyCovContext(_ctxt);
}

void Model::setDrifts(const VectorString &driftSymbols)
{
  if (_driftList == nullptr)
    _driftList = new DriftList();
  else
    delAllDrifts();

  for (int i = 0; i < (int) driftSymbols.size(); i++)
  {
    ADrift *drift = DriftFactory::createDriftBySymbol(driftSymbols[i]);
    addDrift(drift);
  }
}

void Model::delDrift(int rank)
{
  if (_driftList == nullptr) return;
  _driftList->delDrift(rank);
}

void Model::delAllDrifts()
{
  if (_driftList == nullptr) return;
  _driftList->delAllDrifts();
}

const CovAniso* Model::getCova(int icov) const
{
  if (_cova == nullptr) return nullptr;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return nullptr;
  return covalist->getCova(icov);
}
CovAniso* Model::getCova(int icov)
{
  if (_cova == nullptr) return nullptr;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return nullptr;
  return covalist->getCova(icov);
}
int Model::getCovaNumber(bool skipNugget) const
{
  if (_cova == nullptr) return 0;
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return ITEST;
  return covalist->getCovaNumber(skipNugget);
}
const ECov& Model::getCovaType(int icov) const
{
  if (_cova == nullptr) return ECov::UNKNOWN;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return ECov::UNKNOWN;
  return covalist->getType(icov);
}
const MatrixSquareSymmetric& Model::getSillValues(int icov) const
{
  if (_cova == nullptr) return MatrixSquareSymmetric();
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return MatrixSquareSymmetric();
  return covalist->getSill(icov);
}
double Model::getSill(int icov, int ivar, int jvar) const
{
  if (_cova == nullptr) return TEST;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;
  return covalist->getSill(icov, ivar, jvar);
}
double Model::getRange(int icov) const
{
  if (_cova == nullptr) return TEST;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;
  return covalist->getRange(icov);
}
VectorDouble Model::getRanges(int icov) const
{
  if (_cova == nullptr) return VectorDouble();
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return VectorDouble();
  return covalist->getRanges(icov);
}
double Model::getParam(int icov) const
{
  if (_cova == nullptr) return TEST;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;
  return covalist->getParam(icov);
}
bool Model::isCovaFiltered(int icov) const
{
  if (_cova == nullptr) return false;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return false;
  return covalist->isFiltered(icov);
}
bool Model::isStationary() const
{
  if (_cova == nullptr) return false;
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return false;
  return covalist->isStationary();
}
String Model::getCovName(int icov) const
{
  if (_cova == nullptr) return String();
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return String();
  return covalist->getCovName(icov);
}
int Model::getGradParamNumber(int icov) const
{
  if (_cova == nullptr) return ITEST;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return ITEST;
  return covalist->getGradParamNumber(icov);
}
void Model::setSill(int icov, int ivar, int jvar, double value)
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return;
  covalist->setSill(icov, ivar, jvar, value);
}
void Model::setRangeIsotropic(int icov, double range)
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return;
  covalist->setRangeIsotropic(icov, range);
}
void Model::setMarkovCoeffs(int icov, const VectorDouble& coeffs)
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return;
  covalist->setMarkovCoeffs(icov, coeffs);
}
void Model::updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
{
  if (_cova == nullptr) return;
  _cova->updateCovByPoints(icas1, iech1, icas2, iech2);
}
void Model::updateCovByMesh(int imesh)
{
  if (_cova == nullptr) return;
  _cova->updateCovByMesh(imesh);
}
void Model::setCovaFiltered(int icov, bool filtered)
{
  if (_cova == nullptr) return;
  ACovAnisoList* covalist = _castInCovAnisoList(icov);
  if (covalist == nullptr) return;
  covalist->setFiltered(icov, filtered);
}
int Model::hasExternalCov() const
{
  if (_cova == nullptr) return 0;
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return 0;
  for (int icov = 0; icov < (int) covalist->getCovaNumber(); icov++)
  {
    if (covalist->getType(icov) == ECov::FUNCTION) return 1;
  }
  return 0;
}
double Model::getMaximumDistance() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return TEST;
  return covalist->getMaximumDistance();
}
int Model::getCovaMinIRFOrder() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return ITEST;
  return covalist->getCovaMinIRFOrder();
}
bool Model::hasAnam() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return false;
  return covalist->hasAnam();
}
const AAnam* Model::getAnam() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return nullptr;
  return covalist->getAnam();
}
bool Model::isChangeSupportDefined() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return false;
  return covalist->getAnam()->isChangeSupportDefined();
}
void Model::normalize(double sill)
{
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return;
  covalist->normalize(sill);
}
bool Model::hasNugget() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return false;
  return covalist->hasNugget();
}
int Model::getRankNugget() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return -1;
  return covalist->getRankNugget();
}
VectorInt Model::getActiveCovList() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return VectorInt();
  return covalist->getActiveCovList();
}
VectorInt Model::getAllActiveCovList() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return VectorInt();
  return covalist->getAllActiveCovList();
}
bool Model::isAllActiveCovList() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return false;
  return covalist->isAllActiveCovList();
}
void Model::setActiveFactor(int iclass)
{
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return;
  covalist->setActiveFactor(iclass);
}
int Model::getActiveFactor() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return ITEST;
  return covalist->getActiveFactor();
}
int Model::getAnamNClass() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return ITEST;
  return covalist->getAnamNClass();
}

void Model::evalZAndGradients(const SpacePoint &p1,
                              const SpacePoint &p2,
                              double &covVal,
                              VectorDouble &covGp,
                              VectorDouble &covGG,
                              const CovCalcMode *mode,
                              bool flagGrad) const
{
  CovLMGradient* covgrad = dynamic_cast<CovLMGradient *>(_cova);
  if (covgrad != nullptr)
    covgrad->evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
}
void Model::evalZAndGradients(const VectorDouble &vec,
                              double &covVal,
                              VectorDouble &covGp,
                              VectorDouble &covGG,
                              const CovCalcMode *mode,
                              bool flagGrad) const
{
  CovLMGradient* covgrad = dynamic_cast<CovLMGradient *>(_cova);
  if (covgrad != nullptr)
    covgrad->evalZAndGradients(vec, covVal, covGp, covGG, mode, flagGrad);
}

double Model::evalCov(const VectorDouble &incr,
                      int icov,
                      const ECalcMember &member) const
{
  if (_cova == nullptr) return TEST;
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;

  if (member != ECalcMember::LHS && isCovaFiltered(icov))
    return (0.);
  return getCova(icov)->evalIvarIpas(1., incr);
}

/**
 * Define Non-stationary parameters
 * @param anostat ANoStat pointer will be duplicated
 * @return Error return code
 */
int Model::addNoStat(const ANoStat *anostat)
{
  if (anostat == nullptr) return 0;
  if (_cova == nullptr) return 1;
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return 1;
  return covalist->addNoStat(anostat);
}

/**
 * Switch to a Model dedicated to Gradients
 * (transforms it from CovLMC to CovLMGradient)
 */
void Model::switchToGradient()
{
  // If the Model is already dedicated to Gradient: do nothing
  if (isFlagGradient()) return;

  // If no covariance has been defined yet: do nothing
  if (_cova == nullptr)
  {
    _cova = new CovLMGradient(_ctxt.getSpace());
  }
  else
  {
    const ACovAnisoList* covalist = _castInCovAnisoListConst();
    if (covalist == nullptr) return;
    _cova = new CovLMGradient(*covalist);
  }
}

/**
 * Defining an Anamorphosis information for the Model
 * (in fact, this is added to ACovAnisoList part and transforms it from CovLMC to CovLMCAnamorphosis
 * @param anam Pointer to the anamorphosis
 * @param strcnt Array of covariance description used for IR case
 * @return
 */
int Model::setAnam(const AAnam* anam, const VectorInt& strcnt)
{
  if (anam == nullptr)
  {
    messerr("You must define 'anam' beforehand");
    return 1;
  }
  if (hasAnam())
  {
    // ACovAnisoList is already a covLMCAnamorphosis, simply update the anamorphosis
    CovLMCAnamorphosis* cov = dynamic_cast<CovLMCAnamorphosis*>(_cova);
    if (cov == nullptr)
    {
      messerr("Impossible to reach the internal CovLMCAnamorphosis structure");
      return 1;
    }
    cov->setAnam(anam);
  }
  else
  {
    CovLMC* cov = dynamic_cast<CovLMC*>(_cova);
    if (cov == nullptr)
    {
      messerr("Impossible to add 'anam' to the covariance part of the Model");
      messerr("The original covariance is probably not a 'CovLMC'");
      return 1;
    }

    // Initiate a new CovLMCAnamorphosis class
    CovLMCAnamorphosis* newcov = new CovLMCAnamorphosis(*cov, anam, strcnt);

    // Delete the current ACovAnisoList structure
    delete _cova;

    // Replace it by the newly create one (CovLMCAnamorphosis)
    _cova = newcov;
  }
  return 0;
}

void Model::setTapeRange(double range)
{
  CovLMCTapering* covtape = dynamic_cast<CovLMCTapering*>(_cova);
  if (covtape != nullptr) covtape->setTapeRange(range);
}

int Model::unsetAnam()
{
  if (!hasAnam())
  {
    // ACovAnisoList does not have any Anam: do nothing
    return 0;
  }
  CovLMC* cov = dynamic_cast<CovLMC*>(_cova);
  if (cov == nullptr)
  {
    messerr("Impossible to unset 'anam' from the covariance part of the Model");
    messerr("The original covariance is probably not valid");
    return 1;
  }

  // Initiate a new CovLMC class
  CovLMC* newcov = new CovLMC(*cov);

  // Delete the current ACovAnisoList structure
  delete _cova;

  // Replace it by the newly create one (CovLMC)
  _cova = newcov;
  return 0;
}

void Model::_copyCovContext()
{
  if (_cova == nullptr) return;
  ACovAnisoList *covalist = _castInCovAnisoList();
  if (covalist != nullptr) covalist->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

void Model::setMeans(const VectorDouble& mean)
{
  _ctxt.setMean(mean);
  _copyCovContext();
}
void Model::setMean(double mean, int ivar)
{
  _ctxt.setMean(mean, ivar);
  _copyCovContext();
}
void Model::setCovar0s(const VectorDouble& covar0)
{
  _ctxt.setCovar0(covar0);
  _copyCovContext();
}
void Model::setCovar0(int ivar, int jvar, double covar0)
{
  _ctxt.setCovar0(ivar,jvar,covar0);
  _copyCovContext();
}
void Model::setField(double field)
{
  _ctxt.setField(field);
  _copyCovContext();
}

int Model::getNoStatElemNumber() const
{
  if (!isNoStat()) return 0;
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return 0;
  return covalist->getNoStatElemNumber();
}

int Model::addNoStatElem(int igrf,
                         int icov,
                         const EConsElem &type,
                         int iv1,
                         int iv2)
{
  if (!isNoStat()) return 0;
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return 0;
  return covalist->addNoStatElem(igrf, icov, type, iv1, iv2);
}

int Model::addNoStatElems(const VectorString &codes)
{
  if (!isNoStat()) return 0;
  ACovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return 0;
  return covalist->addNoStatElems(codes);
}

CovParamId Model::getCovParamId(int ipar) const
{
  if (!isNoStat()) return CovParamId();
  const ACovAnisoList* covalist = _castInCovAnisoListConst();
  if (covalist == nullptr) return CovParamId();
  return covalist->getCovParamId(ipar);
}

/****************************************************************************/
/*!
 **  Check if the non-stationary Model has a given non-stationary parameter
 **
 ** \return  1 if the given non-stationary parameter is defined; 0 otherwise
 **
 ** \param[in]   type0    Requested type (EConsElem)
 **
 *****************************************************************************/
bool Model::isNostatParamDefined(const EConsElem &type0) const
{
  if (! isNoStat()) return true;
  const ANoStat *nostat = getNoStat();
  return (nostat->isDefinedByType(type0));
}

const DriftList* Model::getDriftList() const
{
  return _driftList;
}
const ADrift* Model::getDrift(int il) const
{
  if (_driftList == nullptr) return nullptr;
  return _driftList->getDrift(il);
}
int Model::getDriftNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getDriftNumber();
}
int Model::getExternalDriftNumber() const
{
  if (_driftList == nullptr) return 0;
  int nfex = 0;
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (getDrift(il)->isDriftExternal()) nfex++;
  }
  return nfex;
}
int Model::getRankFext(int il) const
{
  if (_driftList == nullptr) return ITEST;
  return _driftList->getRankFex(il);
}
bool Model::isDriftSampleDefined(const Db *db,
                                 int ib,
                                 int nech,
                                 const VectorInt &nbgh,
                                 const ELoc &loctype) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isDriftSampleDefined(db,ib,nech,nbgh,loctype);
}
int Model::getDriftEquationNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getDriftEquationNumber();
}
bool Model::isDriftFiltered(unsigned int il) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isFiltered(il);
}
void Model::setDriftFiltered(int il, bool filtered)
{
  if (_driftList == nullptr) return;
  _driftList->setFiltered(il, filtered);
}
bool Model::isDriftDefined(const VectorInt &powers, int rank_fex) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isDriftDefined(powers, rank_fex);
}
bool Model::isDriftDifferentDefined(const VectorInt &powers, int rank_fex) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isDriftDifferentDefined(powers, rank_fex);
}
void Model::setBetaHat(const VectorDouble &betaHat)
{
  if (_driftList == nullptr) return;
  _driftList->setBetaHat(betaHat);
}
int Model::getDriftMaxIRFOrder(void) const
{
  if (_driftList == nullptr) return -1;
  return _driftList->getDriftMaxIRFOrder();
}
VectorVectorDouble Model::getDrifts(const Db *db, bool useSel)
{
  if (_driftList == nullptr) return VectorVectorDouble();
  return _driftList->getDrifts(db, useSel);
}

/**
 * Evaluate a given drift function for a given sample
 * @param db     Db structure
 * @param iech   Rank of the target sample
 * @param il     Rank of the drift function
 * @param member Member type (used to check filtering)
 * @return
 */
double Model::evalDrift(const Db *db,
                        int iech,
                        int il,
                        const ECalcMember &member) const
{
  if (_driftList == nullptr) return TEST;
  return _driftList->evalDrift(db, iech, il, member);
}

VectorDouble Model::evalDriftBySample(const Db *db,
                                      int iech,
                                      const ECalcMember &member) const
{
  if (_driftList == nullptr) return VectorDouble();
  return _driftList->evalDriftBySample(db, iech, member);
}

void Model::evalDriftBySampleInPlace(const Db *db,
                                     int iech,
                                     const ECalcMember &member,
                                     VectorDouble &drftab) const
{
  if (_driftList == nullptr) return;
  _driftList->evalDriftBySampleInPlace(db, iech, member, drftab);
}

/**
 * Returns the value of the normalized covariance (by the variance/covariance value)
 * for a given pair of variables
 * @param hh    Vector of distances
 * @param ivar  Rank of the first variable
 * @param jvar  Rank of the second variable
 * @param codir Direction coefficients
 * @param mode  CovCalcMode structure
 * @return
 */
VectorDouble Model::sampleUnitary(const VectorDouble &hh,
                                  int ivar,
                                  int jvar,
                                  VectorDouble codir,
                                  const CovCalcMode* mode)
{
  if (ivar < 0 || ivar >= getVariableNumber()) return VectorDouble();
  if (jvar < 0 || jvar >= getVariableNumber()) return VectorDouble();
  if (ivar == jvar) return VectorDouble();
  int ndim = getDimensionNumber();
  if (codir.empty())
  {
    (void) GH::rotationGetDirectionDefault(ndim, codir);
  }
  int nh = (int) hh.size();

  double c00 = eval0(ivar, ivar, mode);
  double c11 = eval0(jvar, jvar, mode);
  c00 = sqrt(c00 * c11);
  VectorDouble gg = sample(hh, codir, ivar, jvar, mode);

  for (int i = 0; i < nh; i++)
    gg[i] /= c00;

  return gg;
}

VectorDouble Model::envelop(const VectorDouble &hh,
                            int ivar,
                            int jvar,
                            int isign,
                            VectorDouble codir,
                            const CovCalcMode* mode)
{
  if (ivar < 0 || ivar >= getVariableNumber()) return VectorDouble();
  if (jvar < 0 || jvar >= getVariableNumber()) return VectorDouble();
  if (ivar == jvar) return VectorDouble();
  if (isign != -1 && isign != 1) return VectorDouble();
  int ndim = getDimensionNumber();
  if (codir.empty())
  {
    (void) GH::rotationGetDirectionDefault(ndim, codir);
  }
  int nh = (int) hh.size();
  VectorDouble gg(nh);
  VectorDouble g1 = sample(hh, codir, ivar, ivar, mode);
  VectorDouble g2 = sample(hh, codir, jvar, jvar, mode);

  for (int i = 0; i < nh; i++)
    gg[i] = isign * sqrt(abs(g1[i] * g2[i]));

  return gg;
}

/**
 * Automatic Fitting procedure
 *
 * @param vario       Experimental variogram to be fitted
 * @param types       Vector of ECov integer values
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 * @param mauto       Special parameters for Automatic fitting procedure
 * @param verbose     Verbose option
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fitFromCovIndices(Vario *vario,
                             const VectorECov &types,
                             const Constraints &constraints,
                             const Option_VarioFit& optvar,
                             const Option_AutoFit& mauto,
                             bool verbose)
{
  if (vario == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  for (int is = 0; is < (int) types.size(); is++)
  {
    CovAniso cov = CovAniso(types[is], _ctxt);
    addCov(&cov);
  }

  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}

/**
 * Automatic Fitting procedure from an experimental Variogram
 *
 * @param vario       Experimental variogram to be fitted
 * @param types       Vector of ECov (see remarks)
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 * @param mauto       Special parameters for Automatic fitting procedure (instance of Option_AutoFit), for exemple wmode (type of weighting function)
 * @param verbose     Verbose option
 *
 * @remarks If no list of specific basic structure is specified, the automatic fitting
 * is performed using a single spherical structure by default.
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fit(Vario* vario,
               const VectorECov& types,
               const Constraints& constraints,
               const Option_VarioFit& optvar,
               const Option_AutoFit& mauto,
               bool verbose)
{
  if (vario == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  for (int is = 0; is < (int) types.size(); is++)
  {
    CovAniso cov = CovAniso(types[is], _ctxt);
    addCov(&cov);
  }
  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}

/**
 * Automatic Fitting procedure from A Variogram Map stored on a DbGrid
 *
 * @param dbmap       DbGrid containing the Variogram Map
 * @param types       Vector of ECov
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 * @param mauto       Special parameters for Automatic fitting procedure (instance of Option_AutoFit), for exemple wmode (type of weighting function)
 * @param verbose     Verbose option
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fitFromVMap(DbGrid* dbmap,
                       const VectorECov& types,
                       const Constraints& constraints,
                       const Option_VarioFit& optvar,
                       const Option_AutoFit& mauto,
                       bool verbose)
{
  if (dbmap == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  for (int is = 0; is < (int) types.size(); is++)
  {
    CovAniso cov = CovAniso(types[is], _ctxt);
    addCov(&cov);
  }
  return vmap_auto_fit(dbmap, this, verbose, mauto, constraints, optvar);
}

bool Model::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim = 0;
  int nvar = 0;
  int ncova = 0;
  int nbfl = 0;
  int type = 0;
  int flag_aniso = 0;
  int flag_rotation = 0;

  double field = 0.;
  double range = 0.;
  double param = 0.;
  double value = 0.;

  VectorDouble aniso_ranges;
  VectorDouble aniso_rotmat;

  /* Create the Model structure */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Number of Variables", nvar);
  ret = ret && _recordRead<double>(is, "Field dimension", field);
  ret = ret && _recordRead<int>(is, "Number of Basic Structures", ncova);
  ret = ret && _recordRead<int>(is, "Number of Basic Drift Functions", nbfl);
  if (! ret) return ret;

  /// TODO : Force SpaceRN creation (deserialization doesn't know yet how to manage other space types)
  _ctxt = CovContext(nvar, ndim);
  _ctxt.setField(field);
  _clear();
  _create();

  /* Reading the covariance part and store it into a CovLMC */

  CovLMC covs(_ctxt.getSpace());
  for (int icova = 0; ret && icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    ret = ret && _recordRead<int>(is, "Covariance Type", type);
    ret = ret && _recordRead<double>(is, "Isotropic Range", range);
    ret = ret && _recordRead<double>(is, "Model third Parameter", param);
    ret = ret && _recordRead(is, "Flag for Anisotropy", flag_aniso);
    if (! ret) return ret;
    if (flag_aniso)
    {
      aniso_ranges.resize(ndim);
      // In fact, the file contains the anisotropy coefficients
      // After reading, we must turn them into anisotropic ranges
      for (int idim = 0; idim < ndim; idim++)
        ret = ret && _recordRead<double>(is, "Anisotropy coefficient", aniso_ranges[idim]);
      if (! ret) return ret;
      for (int idim = 0; idim < ndim; idim++)
        aniso_ranges[idim] *= range;

      ret = ret && _recordRead<int>(is, "Flag for Anisotropy Rotation", flag_rotation);
      if (! ret) return ret;
      if (flag_rotation)
      {
        // Warning: the storage in the File is performed by column
        // whereas the internal storage is by column (TODO : ???)
        aniso_rotmat.resize(ndim * ndim);
        int lec = 0;
        for (int idim = 0; ret && idim < ndim; idim++)
          for (int jdim = 0; ret && jdim < ndim; jdim++)
            ret = ret && _recordRead<double>(is, "Anisotropy Rotation Matrix", aniso_rotmat[lec++]);
      }
    }
    if (! ret) return ret;

    CovAniso cova(ECov::fromValue(type), _ctxt);
    cova.setParam(param);
    if (flag_aniso)
    {
      cova.setRanges(aniso_ranges);
      if (flag_rotation) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRangeIsotropic(range);
    covs.addCov(&cova);
  }
  setCovList(&covs);

  /* Reading the drift part */

  DriftList drifts(_ctxt);
  ADrift* drift;
  for (int ibfl = 0; ret && ibfl < nbfl; ibfl++)
  {
    ret = true; // Reset 'ret' to continue reading after previous error...
    String driftname;
    ret = ret && _recordRead<String>(is, "Drift Identifier", driftname);
    drift = DriftFactory::createDriftByIdentifier(driftname);
    drifts.addDrift(drift);
    delete drift;
  }
  setDriftList(&drifts);

  /* Reading the matrix of means (only if nbfl <= 0) */

  if (nbfl <= 0)
    for (int ivar = 0; ret && ivar < nvar; ivar++)
    {
      double mean = 0.;
      ret = ret && _recordRead<double>(is, "Mean of Variable", mean);
      setMean(mean, ivar);
  }

  /* Reading the matrices of sills (optional) */

  for (int icova = 0; icova < ncova && ret; icova++)
  {
    for (int ivar = 0; ret && ivar < nvar; ivar++)
      for (int jvar = 0; ret && jvar < nvar; jvar++)
      {
        ret = ret && _recordRead<double>(is, "Matrix of Sills", value);
        if (ret) setSill(icova, ivar, jvar, value);
      }
  }

  /* Reading the variance-covariance at the origin (optional) */

  for (int ivar = 0; ret && ivar < nvar; ivar++)
    for (int jvar = 0; ret && jvar < nvar; jvar++)
    {
      ret = ret && _recordRead<double>(is, "Variance-covariance at Origin",
                                       value);
      if (ret) setCovar0(ivar, jvar, value);
    }

  return ret;
}

bool Model::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;

  /* Write the Model structure */

  ret = ret && _recordWrite<int>(os, "", getDimensionNumber());
  ret = ret && _recordWrite<int>(os, "", getVariableNumber());
  ret = ret && _recordWrite<double>(os, "General parameters", getField());
  ret = ret && _recordWrite<int>(os, "Number of basic covariance terms", getCovaNumber());
  ret = ret && _recordWrite<int>(os, "Number of drift terms", getDriftNumber());

  /* Writing the covariance part */

  for (int icova = 0; ret && icova < getCovaNumber(); icova++)
  {
    const CovAniso *cova = getCova(icova);
    ret = ret && _recordWrite<int>(os, "", cova->getType().getValue());
    ret = ret && _recordWrite<double>(os, "", cova->getRange());
    ret = ret && _recordWrite<double>(os, "Covariance characteristics", cova->getParam());

    // Writing the Anisotropy information

    ret = ret && _recordWrite<int>(os, "Anisotropy Flag", cova->getFlagAniso());

    if (!cova->getFlagAniso()) continue;

    for (int idim = 0; ret && idim < getDimensionNumber(); idim++)
      ret = ret && _recordWrite<double>(os, "", cova->getAnisoCoeffs(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<int>(os, "Anisotropy Rotation Flag", cova->getFlagRotation());

    if (!cova->getFlagRotation()) continue;

    // Storing the rotation matrix by Column (compatibility)
    for (int idim = 0; ret && idim < getDimensionNumber(); idim++)
      for (int jdim = 0; ret && jdim < getDimensionNumber(); jdim++)
        ret = ret && _recordWrite<double>(os, "", cova->getAnisoRotMat(jdim, idim));
    ret = ret && _commentWrite(os, "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (int ibfl = 0; ret && ibfl < getDriftNumber(); ibfl++)
  {
    const ADrift *drift = getDrift(ibfl);
    ret = ret && _recordWrite<String>(os,"Drift Identifier", drift->getDriftName());
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (getDriftNumber() <= 0)
    for (int ivar = 0; ret && ivar < getVariableNumber(); ivar++)
    {
      ret = ret && _recordWrite<double>(os, "Mean of Variables", getContext().getMean(ivar));
    }

  /* Writing the matrices of sills (optional) */

  for (int icova = 0; ret && icova < getCovaNumber(); icova++)
  {
    for (int ivar = 0; ret && ivar < getVariableNumber(); ivar++)
      for (int jvar = 0; ret && jvar < getVariableNumber(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getSill(icova, ivar, jvar));
    ret = ret && _commentWrite(os, "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (int ivar = 0; ret && ivar < getVariableNumber(); ivar++)
    for (int jvar = 0; ret && jvar < getVariableNumber(); jvar++)
      ret = ret && _recordWrite<double>(os, "", getContext().getCovar0(ivar, jvar));
  ret = ret && _commentWrite(os, "Var-Covar at origin");

  return ret;
}

void Model::_clear()
{
  delete _cova;
  _cova = nullptr;
  delete _driftList;
  _driftList = nullptr;
}

void Model::_create()
{
  // TODO: The next two lines are there in order to allow direct call to
  // model::addCov() and model::addDrift
  // The defaulted types of CovAnisoList and DriftList are assumed
  _cova = new CovLMC(_ctxt.getSpace());
  _driftList = new DriftList(_ctxt);
}

double Model::getTotalSill(int ivar, int jvar) const
{
  return getCovAnisoList()->getTotalSill(ivar, jvar);
}

MatrixSquareSymmetric Model::getTotalSills() const
{
  return getCovAnisoList()->getTotalSill();
}

/**
 * Returns the Ball radius (from the first covariance of _covaList)
 * @return Value of the Ball Radius (if defined, i.e. for Numerical Gradient calculation)
 */
double Model::getBallRadius() const
{
  if (_cova == nullptr) return TEST;

  // Check is performed on the first covariance
  const ACovAnisoList* covalist = _castInCovAnisoListConst(0);
  if (covalist == nullptr) return ITEST;
  const CovAniso* cova = covalist->getCova(0);
  double ball_radius = cova->getBallRadius();
  if (! FFFF(ball_radius)) return ball_radius;
  return 0.;
}

const AnamHermite* Model::getAnamHermite() const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst(0);
  if (covalist == nullptr) return nullptr;
  const AAnam* anam = covalist->getAnam();
  if (anam == nullptr) return nullptr;
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(anam);
  return anamH;
}

Model* Model::duplicate() const
{
  Model* model = new Model(getContext());

  /* Add the list of Covariances */

  model->setCovList(getCovAnisoList());

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  /* Add non-stationarity information */

  model->addNoStat(getNoStat());

  return model;
}

Model* Model::createReduce(const VectorInt& validVars) const
{
  VectorInt localValidVars = VH::filter(validVars, 0, getVariableNumber());
  int nvar = (int) localValidVars.size();
  if (nvar <= 0)
  {
    messerr("Your new Model has no variable left");
    return nullptr;
  }

  Model* model = new Model(*_ctxt.createReduce(validVars));

  /* Add the list of Covariances */

  model->setCovList(getCovAnisoList()->createReduce(validVars));

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  /* Add non-stationarity information */

  model->addNoStat(getNoStat());

  return model;
}

/**
 * Evaluate the Goodness-of_fit of the Model on the Experimental Variogram
 * It is expressed as the average departure between Model and Variogram
 * scaled to the variance.
 * As this variance may be poorly calculated (< gmax / 5), it may be replaced
 * by the largest value (gmax) divided by 2 (highly non_stationary cases).
 * @param vario Experimental variogram
 * @param verbose Verbose flag

 * @return Value for the Goodness-of_fit (as percentage of the total sill)
 */
double Model::gofToVario(const Vario *vario, bool verbose)
{
  int nvar = getVariableNumber();
  int ndir = vario->getDirectionNumber();

  double total = 0.;

  // Loop on the pair of variables

  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double varij  = vario->getVar(ivar, jvar);
      double varmax = vario->getGmax(ivar, jvar);
      // Modify the normalization as variance seems not consistent
      if (ABS(varij) < varmax / 5)
      {
        if (verbose)
          messerr("Variance seems erroneous. It is replaced by Gmax / 2.");
        varij = varmax / 2.;
      }

      // Loop on the variogram directions

      double totdir = 0.;
      for (int idir = 0; idir < ndir; idir++)
      {

        // Read information from Experimental Variogram

        VectorDouble codir = vario->getCodirs(idir);
        VectorDouble sw = vario->getSwVec(idir, ivar, jvar);
        VectorDouble hh = vario->getHhVec(idir, ivar, jvar);
        VectorDouble gexp = vario->getGgVec(idir, ivar, jvar);

        // Evaluate the Model

        int npas = (int) gexp.size();
        VectorDouble gmod = sample(hh, codir, ivar, jvar, &mode);

        // Evaluate the score

        double totpas = 0;
        double scale = 0.;
        for (int ipas = 0; ipas < npas; ipas++)
        {
          if (sw[ipas] <= 0 || hh[ipas] <= 0.) continue;
          double ecart = sw[ipas] * ABS(gexp[ipas] - gmod[ipas]) / hh[ipas];
          totpas += ecart;
          scale  += sw[ipas] / hh[ipas];
        }
        totpas  = totpas / scale;
        totdir += totpas;
      }
      totdir /= (double) ndir;
      totdir /= varij;
      total  += ABS(totdir);
    }
  total = 100. * total / (double) (nvar * nvar);
  return total;
}

/**
 * Printout of statement concerning the Quality of the GOF
 * @param gof        Value of the Gof
 * @param byValue    true: display GOF value; false: print its quality level
 * @param thresholds Vector giving the Quality thresholds
 */
void Model::gofDisplay(double gof, bool byValue, const VectorDouble& thresholds)
{
  message("Goodness-of-fit (as a percentage of the variance)");
  if (byValue)
  {
    message(" = %5.2lf\n", gof);
    return;
  }
  int nclass = (int)thresholds.size();
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (gof < thresholds[iclass])
    {
      message(" corresponds to level #%d (1 for very good)\n", iclass + 1);
      return;
    }
  }
}

const EModelProperty& Model::getCovMode() const
{
  ACovAnisoList* covs;
  if (_cova == nullptr) return EModelProperty::NONE;

  covs = dynamic_cast<CovLMCTapering*>(_cova);
  if (covs != nullptr) return EModelProperty::TAPE;

  covs = dynamic_cast<CovLMCConvolution*>(_cova);
  if (covs != nullptr) return EModelProperty::CONV;

  covs = dynamic_cast<CovLMCAnamorphosis*>(_cova);
  if (covs != nullptr) return EModelProperty::ANAM;

  covs = dynamic_cast<CovLMGradient*>(_cova);
  if (covs != nullptr) return EModelProperty::GRAD;

  return EModelProperty::NONE;
}

bool Model::isFlagLinked() const
{
  if (_driftList == nullptr) return false;
  return _driftList->isFlagLinked();
}

bool Model::isFlagGradient() const
{
  if (_cova == nullptr) return false;
  return getCovMode() == EModelProperty::GRAD;
}

bool Model::isFlagGradientNumerical() const
{
  if (! isFlagGradient()) return false;

  // Check is performed on the first covariance
  const ACovAnisoList* covalist = _castInCovAnisoListConst(0);
  if (covalist == nullptr) return false;
  const CovGradientNumerical* cova = dynamic_cast<const CovGradientNumerical*>(covalist->getCova(0));
  return (cova != nullptr);
}

bool Model::isFlagGradientFunctional() const
{
  if (! isFlagGradient()) return false;

  // Check is performed on the first covariance
  const ACovAnisoList* covalist = _castInCovAnisoListConst(0);
  if (covalist == nullptr) return false;
  const CovGradientFunctional* cova = dynamic_cast<const CovGradientFunctional*>(covalist->getCova(0));
  return (cova != nullptr);
}

/****************************************************************************/
/*!
 **  Evaluate the drift with a given sample and a given variable
 **  The value is scaled by 'coeffs'
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 ** \param[in]  coeffs  Vector of coefficients
 **
 *****************************************************************************/
double Model::evalDriftVarCoef(const Db *db,
                               int iech,
                               int ivar,
                               const VectorDouble &coeffs) const
{
  if (_driftList == nullptr)
  {
    double mean = getMean(ivar);
    return mean;
  }
  double drift = 0.;
  for (int ib = 0, nfeq = getDriftEquationNumber(); ib < nfeq; ib++)
    drift += evalDriftValue(db, iech, ivar, ib, ECalcMember::LHS) * coeffs[ib];
  return drift;
}

/**
 * A vector of the drift evaluation (for all samples)
 * @param db     Db structure
 * @param coeffs Vector of drift coefficients
 * @param ivar   Variable rank (used for constant drift value)
 * @param useSel When TRUE, only non masked samples are returned
 * @return The vector of values
 *
 * @remark When no drift is defined, a vector is returned filled with the variable mean
 */
VectorDouble Model::evalDriftVarCoefs(const Db *db,
                                      const VectorDouble &coeffs,
                                      int ivar,
                                      bool useSel) const
{
  VectorDouble vec;
  if (_driftList == nullptr)
  {
    if (db == nullptr) return vec;
    int nech = db->getSampleNumber(useSel);
    double mean = getMean(ivar);
    vec = VectorDouble(nech, mean);
  }
  else
  {
    vec = _driftList->evalDriftCoefs(db, coeffs, useSel);
  }
  return vec;
}

double Model::evalDriftValue(const Db *db,
                             int iech,
                             int ivar,
                             int ib,
                             const ECalcMember &member) const
{
  if (_driftList == nullptr) return TEST;
  return _driftList->evalDriftValue(db, iech, ivar, ib, member);
}

VectorECov Model::initCovList(const VectorInt & covranks)
{
  VectorECov list;

  for (int i = 0; i < (int) covranks.size(); i++)
  {
    ECov ec = ECov::fromValue(covranks[i]);
    if (ec == ECov::UNKNOWN)
    {
      ECov::printAll();
      list.clear();
      break;
    }
    list.push_back(ec);
  }
  return list;
}

bool Model::isValid() const
{
  // Covariances: there should be some defined
  if (_cova == nullptr)
  {
    messerr("Model is not valid: no covariance has been defined");
    return false;
  }

  // Drifts: there should be valid
  if (_driftList != nullptr)
  {
    if (! _driftList->isValid()) return false;
  }

  // Check the consistency between the Covariance and the Drift parts
  int irf_drift = getDriftMaxIRFOrder();
  int irf_cova = getCovaMinIRFOrder();
  if (irf_cova > irf_drift)
  {
    messerr("Model if invalid due to IRF degree inconsistency");
    messerr("- Covariance implies a order >= %d", irf_cova);
    messerr("- Drift implies a order %d", irf_drift);
    messerr("(Order -1 stands for strict stationarity)");
    return false;
  }
  return true;
}

const ACovAnisoList* Model::getCovAnisoList() const
{
  return _castInCovAnisoListConst();
}

/**
 * This internal function tries to cast the member '_cova' into a pointer to ACovAnisoList
 * and checks the validity of the argument 'icov' which gives the rank within this list
 * @param icov Rank of the CovAniso (to be checked if >= 0)
 * @return 'nullptr' if not valid cast (the error message is printed internally)
 */
const ACovAnisoList* Model::_castInCovAnisoListConst(int icov) const
{
  // Check the cast procedure
  const ACovAnisoList* covalist = dynamic_cast<const ACovAnisoList*>(_cova);
  if (covalist == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovAnisoList");
    return nullptr;
  }
  if (icov < 0) return covalist;

  // Check the rank
  if (icov >= covalist->getCovaNumber())
  {
    messerr("The rank 'icov' (%d) is not valid. The CovAnisoList contains %d covariances",
            icov, covalist->getCovaNumber());
    return nullptr;
  }
  return covalist;
}

ACovAnisoList* Model::_castInCovAnisoList(int icov)
{
  // Check the cast procedure
  ACovAnisoList* covalist = dynamic_cast<ACovAnisoList*>(_cova);
  if (covalist == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovAnisoList");
    return nullptr;
  }
  if (icov < 0) return covalist;

  // Check the rank
  if (icov >= covalist->getCovaNumber())
  {
    messerr("The rank 'icov' (%d) is not valid. The CovAnisoList contains %d covariances",
            icov, covalist->getCovaNumber());
    return nullptr;
  }
  return covalist;
}

CovAniso Model::extractCova(int icov) const
{
  const ACovAnisoList* covalist = _castInCovAnisoListConst(icov);
  if (covalist == nullptr) return CovAniso(ECov::UNKNOWN, _ctxt);
  return covalist->extractCova(icov);
}

/****************************************************************************/
/*!
 **  Calculate the variogram map from a Model
 **  (presented as Variogram, not Covariance)
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid      Grid structure
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int Model::buildVmapOnDbGrid(DbGrid *dbgrid, const NamingConvention &namconv) const
{
  if (dbgrid == nullptr) return 1;

  /* Initializations */

  int ndim = dbgrid->getNDim();
  int nvar = dbgrid->getLocNumber(ELoc::Z);
  int nv2  = nvar * (nvar + 1) / 2;

  /* Create the variables in the Variogram Map file */

  int iptr = dbgrid->addColumnsByConstant(nv2, 0.);
  if (iptr < 0) return 1;

  /* Loop on the grid nodes */

  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  VectorInt center = dbgrid->getCenterIndices();
  VectorDouble dincr(ndim);
  VectorInt indices(ndim);
  MatrixSquareGeneral mat;
  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (! dbgrid->isActive(iech)) continue;
    dbgrid->rankToIndice(iech, indices);

    for (int idim = 0; idim < ndim; idim++)
      dincr[idim] = (indices[idim] - center[idim]) * dbgrid->getDX(idim);

    // Evaluate the variogram map
    mat = evalNvarIpasIncr(dincr, &mode);

    int ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ecr++)
        dbgrid->setArray(iech, iptr+ecr, mat.getValue(ivar, jvar));
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, iptr, "Model", nv2);
  return 0;
}

/****************************************************************************/
/*!
 **  Stabilize the model (in the monovariate case)
 **
 ** \return  Error returned code
 **
 ** \param[in]  percent  Percentage of nugget effect added
 ** \param[in]  verbose  true for a verbose output
 **
 ** \remark  If the model only contains GAUSSIAN structures, add
 ** \remark  a NUGGET EFFECT structure with a sill equal to a percentage
 ** \remark  of the total sill of the GAUSSIAN component(s)
 **
 ** \remark  This function does not do anything in the multivariate case
 **
 *****************************************************************************/
int Model::stabilize(double percent, bool verbose)
{
  int nvar = getVariableNumber();
  if (nvar > 1) return 0;
  if (percent <= 0.) return 0;
  int ncov = getCovaNumber();

  /* Check if the model only contains GAUSSIAN components */

  double total = 0.;
  for (int icov = 0; icov < ncov; icov++)
  {
    if (getCova(icov)->getType() != ECov::GAUSSIAN) return (0);
    total += getSill(icov, 0, 0);
  }
  total = total * percent / 100.;

  /* Update each Gaussian component */

  for (int icov = 0; icov < ncov; icov++)
    setSill(icov, 0, 0, 1. - total);

  /* Add a NUGGET EFFECT component */

  addCovFromParam(ECov::NUGGET, 0., total);

  /* Printout */

  if (verbose)
  {
    message("The model which only contains Gaussian components\n");
    message("has been stabilized by adding a small Nugget Effect\n");
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Normalize the model
 **
 ** \param[in]  verbose  true for a verbose output
 **
 *****************************************************************************/
int Model::standardize(bool verbose)

{
  int nvar = getVariableNumber();
  int ncov = getCovaNumber();
  VectorDouble total(nvar,0.);

  /* Calculate the total sills for each variable */

  bool flag_norm = false;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    total[ivar] = getTotalSill(ivar, ivar);
    if (isZero(total[ivar])) return 1;
    total[ivar] = sqrt(total[ivar]);
    if (ABS(total[ivar] - 1.) > EPSILON6) flag_norm = true;
  }

  /* Scale the different sills for the different variables */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      for (int icov = 0; icov < ncov; icov++)
      {
        double sill = getSill(icov,ivar, jvar);
        sill /= total[ivar] * total[jvar];
        setSill(icov, ivar, jvar, sill);
      }

  /* Printout */

  if (verbose && flag_norm)
  {
    message("The model has been normalized\n");
    for (int ivar = 0; ivar < nvar; ivar++)
      message("- Variable %d : Scaling factor = %lf\n", ivar + 1,
              total[ivar] * total[ivar]);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances
 **
 ** \return  Array containing the model values
 **
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  codir      Array giving the direction coefficients (optional)
 ** \param[in]  h          Vector of increments
 ** \param[in]  mode       CovCalcMode structure
 ** \param[in]  covint     Non-stationary parameters
 **
 *****************************************************************************/
VectorDouble Model::sample(const VectorDouble &h,
                           const VectorDouble &codir,
                           int ivar,
                           int jvar,
                           const CovCalcMode *mode,
                           const CovInternal *covint)
{
  int nh   = (int) h.size();
  int ndim = getDimensionNumber();
  int nvar = getVariableNumber();

  /* Core allocation */

  VectorDouble d1(ndim);
  MatrixSquareGeneral covtab(nvar);

  /* Get the normalized direction vector */

  VectorDouble codir_loc = codir;
  if (codir_loc.empty())
  {
    (void) GH::rotationGetDirectionDefault(ndim, codir_loc);
  }
  else
  {
    VH::normalizeCodir(ndim, codir_loc);
  }

  /* Loop on the lags */

  VectorDouble g(nh);
  for (int ih = 0; ih < nh; ih++)
  {
    double hh = h[ih];
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = hh * codir_loc[idim];
    evaluateMatInPlace(covint, d1, covtab, true, 1., mode);
    g[ih] = covtab.getValue(ivar, jvar);
  }
  return g;
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances
 **
 ** \return  The model value
 **
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  mode       CovCalcMode structure
 ** \param[in]  codir      Array giving the direction coefficients (optional)
 ** \param[in]  hh         Vector of increments
 **
 *****************************************************************************/
double Model::evaluateOneIncr(double hh,
                              const VectorDouble &codir,
                              int ivar,
                              int jvar,
                              const CovCalcMode *mode)
{
  int ndim = getDimensionNumber();
  int nvar = getVariableNumber();

  /* Core allocation */

  VectorDouble d1(ndim);
  MatrixSquareGeneral covtab(nvar);

  /* Normalize the direction vector codir */

  /* Get the normalized direction vector */

  VectorDouble codir_loc = codir;
  if (codir_loc.empty())
  {
    (void) GH::rotationGetDirectionDefault(ndim, codir_loc);
  }
  else
  {
    VH::normalizeCodir(ndim, codir_loc);
  }

  for (int idim = 0; idim < ndim; idim++)
    d1[idim] = hh * codir_loc[idim];
  evaluateMatInPlace(nullptr, d1, covtab, true, 1., mode);
  return covtab.getValue(ivar, jvar);
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the generic internal function
 **  It can be called for stationary or non-stationary case
 **
 ** \param[in]  covint       Internal structure for non-stationarityAddress for the next term after the drift
 **                          or NULL (for stationary case)
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Multiplicative weight
 ** \param[in]  d1           Distance vector
 ** \param[out] covtab       Covariance array
 **
 *****************************************************************************/
void Model::evaluateMatInPlace(const CovInternal *covint,
                               const VectorDouble &d1,
                               MatrixSquareGeneral &covtab,
                               bool flag_init,
                               double weight,
                               const CovCalcMode *mode)
{
  // Load the non-stationary parameters if needed

  if (isNoStat() && covint != nullptr)
  {
    updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                      covint->getIcas2(), covint->getIech2());
  }

  // Evaluate the Model

  MatrixSquareGeneral mat = evalNvarIpas(1., d1, mode);

  int nvar = getVariableNumber();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = weight * mat.getValue(ivar, jvar);
      if (flag_init)
        covtab.setValue(ivar,jvar,value);
      else
        covtab.updValue(ivar,jvar, EOperator::ADD, value);
      }
}

/*****************************************************************************/
/*!
 **  Returns the covariance for an increment
 **  This is the generic internal function
 **  It can be called for stationary or non-stationary case
 **
 ** \param[in]  covint       Internal structure for non-stationarityAddress for the next term after the drift
 **                          or NULL (for stationary case)
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  weight       Multiplicative weight
 ** \param[in]  d1           Distance vector
 **
 *****************************************************************************/
double Model::evaluateOneGeneric(const CovInternal *covint,
                                 const VectorDouble &d1,
                                 double weight,
                                 const CovCalcMode *mode)
{
  // Load the non-stationary parameters if needed

  if (isNoStat() && covint != nullptr)
  {
    updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                      covint->getIcas2(), covint->getIech2());
  }

  // Return the (weighted) Model value

  return (weight * evalIvarIpas(1, d1, 0, 0, mode));
}

/****************************************************************************/
/*!
 **  Evaluate the model on a Db
 **
 ** \param[in]  db         Db structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  mode       CovCalcMode structure
 **
 *****************************************************************************/
VectorDouble Model::evaluateFromDb(Db *db,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode *mode)
{
  if (getDimensionNumber() != db->getNDim())
  {
    messerr("Dimension of the Db (%d) does not match dimension of the Model (%d)",
            db->getNDim(), getDimensionNumber());
    return VectorDouble();
  }
  int ndim = getDimensionNumber();
  int nvar = getVariableNumber();
  int nech = db->getSampleNumber();

  /* Core allocation */

  VectorDouble d1(ndim,0.);
  MatrixSquareGeneral covtab(nvar);
  VectorDouble gg(nech, TEST);

  /* Loop on the lags */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    db_sample_load(db, ELoc::X, iech, d1.data());
    evaluateMatInPlace(nullptr, d1, covtab, true, 1., mode);
    gg[iech] = covtab.getValue(ivar, jvar);
  }
  return gg;
}

/*****************************************************************************/
/*!
 **  Returns the standard deviation at a given increment for a given model
 **  between two samples of two Dbs
 **
 ** \param[in]  db1         First Db
 ** \param[in]  iech1       Rank in the first Db
 ** \param[in]  db2         Second Db
 ** \param[in]  iech2       Rank in the second Db
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  factor      Multiplicative factor for standard deviation
 ** \param[in]  mode        CovCalcMode structure
 **
 *****************************************************************************/
double Model::calculateStdev(Db *db1,
                             int iech1,
                             Db *db2,
                             int iech2,
                             bool verbose,
                             double factor,
                             const CovCalcMode *mode)
{

  /* Covariance at origin */

  int ndim = db1->getNDim();
  VectorDouble dd(ndim, 0.);
  double c00 = evaluateOneGeneric(nullptr, dd, 1., mode);

  /* Covariance at increment */

  if (db1->getDistanceVec(iech1, iech2, dd, db2)) return TEST;
  double cov = evaluateOneGeneric(nullptr, dd, 1., mode);
  double stdev = factor * sqrt(c00 - cov);

  if (verbose)
  {
    message("Db1(%d) - Db2(%d)", iech1 + 1, iech2 + 1);
    message(" - Incr=");
    for (int idim = 0; idim < ndim; idim++)
      message(" %lf", dd[idim]);
    message(" - c(0)=%lf cov=%lf stdev=%lf\n", c00, cov, stdev);
  }
  return stdev;
}

/*****************************************************************************/
/*!
 **  Update the Model in the case of Non-stationary parameters
 **  This requires the knowledge of the two end-points
 **
 ** \param[in]  covint       Internal structure for non-stationarity
 **                          or NULL (for stationary case)
 **
 *****************************************************************************/
void Model::nostatUpdate(CovInternal *covint)
{
  if (covint == NULL) return;
  updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                    covint->getIcas2(), covint->getIech2());
}

/**
 * Compute the log-likelihood (traditional method)
 *
 * @param db  Db structure where variable are loaded from
 * @param verbose Verbose flag
 *
 * @remarks The calculation considers all the active samples.
 * @remarks It can work in multivariate case with or without drift conditions (linked or not)
 * @remarks The algorithm is stopped (with a message) in the heterotopic case
 * // TODO; improve for heterotopic case
 */
double Model::computeLogLikelihood(Db* db, bool verbose)
{
  int nvar = db->getLocatorNumber(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
    return TEST;
  }
  if (_driftList == nullptr)
  {
    messerr("This function only makes sense when Drift is defined");
    return TEST;
  }
  if (! db->isAllIsotopic())
  {
    messerr("This method is only available for isotopic data set");
    return TEST;
  }
  int nech = db->getSampleNumber(true);
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples   = %d\n", nech);
    message("- Number of variables        = %d\n", nvar);
    message("- Number of drift conditions = %d\n", getDriftEquationNumber());
  }

  // Calculate the covariance matrix C and perform its Cholesky decomposition
  MatrixSquareSymmetric cov = evalCovMatrixSymmetric(db);
  if (cov.computeCholesky())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
    return TEST;
  }

  // Extract the matrix of drifts at samples X
  MatrixRectangular X = evalDriftMatrix(db);

  // Calculate Cm1X = Cm1 * X
  MatrixRectangular Cm1X;
  if (cov.solveCholeskyMat(X, Cm1X))
  {
    messerr("Problem when solving a Linear System after Cholesky decomposition");
    return TEST;
  }

  // Calculate XCm1X = Xt * Cm1X
  MatrixSquareSymmetric* XCm1X = MatrixFactory::prodMatMat<MatrixSquareSymmetric>(&X, &Cm1X, true, false);

  // Establish the vector of multivariate data
  VectorDouble Z = db->getColumnsByLocator(ELoc::Z, true, true);

  // Construct ZCm1X = Zt * Cm1X and perform its Cholesky decomposition
  VectorDouble ZCm1X = Cm1X.prodVecMat(Z);
  if (XCm1X->computeCholesky())
  {
    messerr("Cholesky decomposition of XCm1X matrix failed");
    delete XCm1X;
    return TEST;
  }

  // Calculate beta = (XCm1X)-1 * ZCm1X
  VectorDouble beta;
  if (XCm1X->solveCholesky(ZCm1X, beta))
  {
    messerr("Error when calculating Maximum Likelihood criterion");
    delete XCm1X;
    return TEST;
  }
  setBetaHat(beta);
  delete XCm1X;

  if (verbose)
  {
    VH::display("Optimal Drift coefficients = ", beta);
  }

  // Center the data by the optimal drift: Z = Z - beta * X
  VH::subtractInPlace(Z, X.prodMatVec(beta));

  // Calculate Cm1Zc = Cm1 * Z
  VectorDouble Cm1Zc;
  if (cov.solveCholesky(Z, Cm1Zc))
  {
    messerr("Error when calculating Cm1Zc");
    return TEST;
  }

  // Calculate the log-determinant
  double logdet = cov.computeCholeskyLogDeterminant();

  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(Z, Cm1Zc);

  // Derive the log-likelihood
  double loglike = -0.5 * (logdet + quad + nvar * nech * log(GV_PI));

  // Optional printout
  if (verbose)
  {
    message("Log-Determinant = %lf\n", logdet);
    message("Quadratic term = %lf\n", quad);
    message("Log-likelihood = %lf\n", loglike);
  }
  return loglike;
}
