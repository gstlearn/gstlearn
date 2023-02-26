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
#include <Geometry/GeometryHelper.hpp>
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Enum/ECov.hpp"
#include "Enum/EModelProperty.hpp"

#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovGradientFunctional.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Model/ANoStat.hpp"
#include "Model/NoStatArray.hpp"
#include "Db/Db.hpp"
#include <math.h>

Model::Model(const CovContext &ctxt)
    : AStringable(),
      ASerializable(),
      _covaList(nullptr),
      _driftList(nullptr),
      _noStat(nullptr),
      _ctxt(ctxt)
{
  _create();
}

Model::Model(int nvar, int ndim)
    : AStringable(),
      ASerializable(),
      _covaList(nullptr),
      _driftList(nullptr),
      _noStat(nullptr),
      _ctxt()
{
  SpaceRN space = SpaceRN(ndim);
  _ctxt = CovContext(nvar, &space);
  _create();
}

Model::Model(const Model &m)
    : AStringable(m),
      ASerializable(m),
      _covaList(nullptr),
      _driftList(nullptr),
      _noStat(nullptr),
      _ctxt(m._ctxt)
{
  if (m._covaList != nullptr)
    _covaList = dynamic_cast<ACovAnisoList*>(m._covaList->clone());
  if (m._driftList != nullptr)
    _driftList = m._driftList->clone();
  if (m._noStat != nullptr)
    _noStat = dynamic_cast<ANoStat*>(m._noStat->clone());
}

Model& Model::operator=(const Model &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    if (m._covaList != nullptr)
      _covaList = dynamic_cast<ACovAnisoList*>(m._covaList->clone());
    if (m._driftList != nullptr)
      _driftList = m._driftList->clone();
    if (m._noStat != nullptr)
      _noStat = dynamic_cast<ANoStat*>(m._noStat->clone());
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
  model->addCovFromParam(type, range, sill, param, ranges, sills, angles,
                         flagRange);

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
  sstr << "Number of drift equation(s)  = " << getDriftEquationNumber()
       << std::endl;

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

  if (isNoStat())
  {
    sstr << _noStat->toString();
  }
  return sstr.str();
}

void Model::delCova(int rank)
{
  if (_covaList == nullptr) return;
  _covaList->delCov(rank);
}

void Model::delAllCovas()
{
  if (_covaList == nullptr) return;
  _covaList->delAllCov();
}

/**
 * Add a list of Covariances. This operation cleans any previously stored covariance
 * @param covalist List of Covariances to be added
 */
void Model::setCovList(const ACovAnisoList* covalist)
{
  if (covalist == nullptr) return;
  if (_covaList != nullptr) delete _covaList;
  _covaList = dynamic_cast<ACovAnisoList*>(covalist->clone());
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
  if (_covaList == nullptr) return;
  _covaList->addCov(cov);
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
      cov.setRange(range);
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
  return;
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
  if (_driftList != nullptr) delete _driftList;
  _driftList = driftlist->clone();
}

/**
 * Define the list of drift functions for:
 * - a given degree of the IRF
 * - a given number of external drifts
 * @param order Order of the IRF
 * @param nfex  Number of External Drifts
 *
 * @remark This method deletes any pre-existing drift functions
 */
void Model::setDriftIRF(int order, int nfex)
{
  if (_driftList == nullptr)
    _driftList = new DriftList();
  _driftList->setDriftIRF(order, nfex, _ctxt);
}

void Model::addDrift(const ADriftElem *drift)
{
  if (drift == nullptr) return;
  if (_driftList == nullptr) return;
  _driftList->addDrift(drift);
}

void Model::setDrifts(const VectorString &driftSymbols)
{
  if (_driftList == nullptr)
    _driftList = new DriftList();
  else
    delAllDrifts();

  for (int i = 0; i < (int) driftSymbols.size(); i++)
  {
    int rank_fex = 0;
    EDrift type = DriftFactory::identifyDrift(driftSymbols[i], &rank_fex, _ctxt);
    ADriftElem *drift = DriftFactory::createDriftFunc(type, _ctxt, rank_fex);
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

const CovAniso* Model::getCova(unsigned int icov) const
{
  if (_covaList == nullptr) return nullptr;
  return _covaList->getCova(icov);
}
CovAniso* Model::getCova(unsigned int icov)
{
  if (_covaList == nullptr) return nullptr;
  return _covaList->getCova(icov);
}
int Model::getCovaNumber() const
{
  if (_covaList == nullptr) return 0;
  return _covaList->getCovNumber();
}
const ECov& Model::getCovaType(int icov) const
{
  if (_covaList == nullptr) return ECov::UNKNOWN;
  return _covaList->getType(icov);
}
const MatrixSquareSymmetric& Model::getSill(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getSill(icov);
}
double Model::getSill(int icov, int ivar, int jvar) const
{
  if (_covaList == nullptr) return TEST;
  return _covaList->getSill(icov, ivar, jvar);
}
double Model::getParam(int icov) const
{
  if (_covaList == nullptr) return TEST;
  return _covaList->getParam(icov);
}
bool Model::isCovaFiltered(int icov) const
{
  if (_covaList == nullptr) return false;
  return _covaList->isFiltered(icov);
}
String Model::getCovName(int icov) const
{
  if (_covaList == nullptr) return String();
  return _covaList->getCovName(icov);
}
int Model::getGradParamNumber(int icov) const
{
  if (_covaList == nullptr) return ITEST;
  return _covaList->getGradParamNumber(icov);
}
void Model::setSill(int icov, int ivar, int jvar, double value)
{
  if (_covaList == nullptr) return;
  _covaList->setSill(icov, ivar, jvar, value);
}
void Model::setCovaFiltered(int icov, bool filtered)
{
  if (_covaList == nullptr) return;
  _covaList->setFiltered(icov, filtered);
}
int Model::hasExternalCov() const
{
  if (_covaList == nullptr) return 0;
  for (int icov = 0; icov < (int) _covaList->getCovNumber(); icov++)
  {
    if (_covaList->getType(icov) == ECov::FUNCTION) return 1;
  }
  return 0;
}

/**
 * Define Non-stationary parameters
 * @param anostat ANoStat pointer will be duplicated
 * @return Error return code
 */
int Model::addNoStat(const ANoStat *anostat)
{
  if (anostat == nullptr) return 0;
  if (getDimensionNumber() > 3)
  {
    messerr("Non stationary model is restricted to Space Dimension <= 3");
    return 1;
  }

  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    int icov = getNoStatElemIcov(ipar);
    EConsElem type = getNoStatElemType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the Model definition

    if (icov < 0 || icov >= getCovaNumber())
    {
      messerr("Invalid Covariance rank (%d) for the Non-Stationary Parameter (%d)",
              icov, ipar);
      return 1;
    }
    if (type == EConsElem::PARAM)
    {
      messerr("The current methodology does not handle constraint on third parameter");
      return 1;
    }
  }

  if (_noStat != nullptr) delete _noStat;
  _noStat = dynamic_cast<ANoStat*>(anostat->clone());
  return 0;
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
  if (_covaList == nullptr)
  {
    _covaList = new CovLMGradient(_ctxt.getSpace());
  }
  else
  {
    _covaList = new CovLMGradient(*_covaList);
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
    CovLMCAnamorphosis* cov = dynamic_cast<CovLMCAnamorphosis*>(_covaList);
    if (cov == nullptr)
    {
      messerr("Impossible to reach the internal CovLMCAnamorphosis structure");
      return 1;
    }
    cov->setAnam(anam);
  }
  else
  {
    CovLMC* cov = dynamic_cast<CovLMC*>(_covaList);
    if (cov == nullptr)
    {
      messerr("Impossible to add 'anam' to the covariance part of the Model");
      messerr("The original covariance is probably not a 'CovLMC'");
      return 1;
    }

    // Initiate a new CovLMCAnamorphosis class
    CovLMCAnamorphosis* newcov = new CovLMCAnamorphosis(*cov, anam, strcnt);

    // Delete the current ACovAnisoList structure
    delete _covaList;

    // Replace it by the newly create one (CovLMCAnamorphosis)
    _covaList = newcov;
  }
  return 0;
}

void Model::_copyCovContext()
{
  if (_covaList != nullptr) _covaList->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

void Model::setMeans(const VectorDouble& mean)
{
  _ctxt.setMean(mean);
  _copyCovContext();
}
void Model::setMean(int ivar, double mean)
{
  _ctxt.setMean(ivar, mean);
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

int Model::isNoStat() const
{
  return _noStat != nullptr;
}

int Model::getNoStatElemNumber() const
{
  if (!isNoStat()) return 0;
  return _noStat->getNoStatElemNumber();
}

int Model::addNoStatElem(int igrf,
                         int icov,
                         const EConsElem &type,
                         int iv1,
                         int iv2)
{
  if (!isNoStat()) return 0;
  return _noStat->addNoStatElem(igrf, icov, type, iv1, iv2);
}

int Model::addNoStatElems(const VectorString &codes)
{
  if (!isNoStat()) return 0;
  return _noStat->addNoStatElems(codes);
}

CovParamId Model::getCovParamId(int ipar) const
{
  if (!isNoStat())
    my_throw("Nostat is not defined and cannot be returned");
  return _noStat->getItems(ipar);
}

int Model::getNoStatElemIcov(int ipar)
{
  if (!isNoStat()) return ITEST;
  return _noStat->getICov(ipar);
}
const EConsElem& Model::getNoStatElemType(int ipar)
{
  if (!isNoStat()) return EConsElem::UNKNOWN;
  return _noStat->getType(ipar);
}
const DriftList* Model::getDriftList() const
{
  return _driftList;
}
const ADriftElem* Model::getDrift(int il) const
{
  if (_driftList == nullptr) return nullptr;
  return _driftList->getDrift(il);
}
ADriftElem* Model::getDrift(int il)
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
    if (getDrift(il)->getType() == EDrift::F) nfex++;
  }
  return nfex;
}
const EDrift& Model::getDriftType(int il) const
{
  if (_driftList == nullptr) return EDrift::UNKNOWN;
  return _driftList->getType(il);
}
int Model::getRankFext(int il) const
{
  if (_driftList == nullptr) return ITEST;
  return _driftList->getRankFex(il);
}
const VectorDouble& Model::getCoefDrifts() const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->getCoefDrift();
}
double Model::getCoefDrift(int ivar, int il, int ib) const
{
  if (_driftList == nullptr) return TEST;
  return _driftList->getCoefDrift(ivar, il, ib);
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
bool Model::isDriftDefined(const EDrift& type0) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isDriftDefined(type0);
}
bool Model::isDriftDifferentDefined(const EDrift& type0) const
{
  if (_driftList == nullptr) return false;
  return _driftList->isDriftDifferentDefined(type0);
}
void Model::setCoefDrift(int ivar, int il, int ib, double coeff)
{
  if (_driftList == nullptr) return;
  _driftList->setCoefDrift(ivar, il, ib, coeff);
}
void Model::setCoefDriftByRank(int rank, double coeff)
{
  if (_driftList == nullptr) return;
  _driftList->setCoefDriftByRank(rank, coeff);
}
void Model::setDriftFiltered(int il, bool filtered)
{
  if (_driftList == nullptr) return;
  _driftList->setFiltered(il, filtered);
}
VectorDouble Model::getDrift(const Db *db, int ib, bool useSel)
{
  if (_driftList == nullptr) return VectorDouble();
  return _driftList->getDrift(db, ib, useSel);
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
  if (member != ECalcMember::LHS && isDriftFiltered(il))
    return 0.;
  else
  {
    if (_driftList == nullptr) return TEST;
    ADriftElem *drift = _driftList->getDrift(il);
    if (drift != nullptr) return drift->eval(db, iech);
  }
  return TEST;
}

VectorDouble Model::evalDriftVec(const Db* db,
                                 int iech,
                                 const ECalcMember& member) const
{
  int ndrift = getDriftNumber();
  VectorDouble drftab(ndrift);
  for (int il = 0; il < ndrift; il++)
     drftab[il] = evalDrift(db, iech, il, member);
  return drftab;
}

void Model::evalDriftVecInPlace(const Db* db,
                                int iech,
                                const ECalcMember& member,
                                VectorDouble& drftab) const
{
  int ndrift = getDriftNumber();
  for (int il = 0; il < ndrift; il++)
     drftab[il] = evalDrift(db, iech, il, member);
}

/**
 * A vector of the drift evaluation
 * @param db     Db structure
 * @param coeffs Vector of drift coefficients
 * @param ivar   Variable rank (used for constant drift value)
 * @param useSel When TRUE, only non masked samples are returned
 * @return The vector of values
 * @remark When no drift is defined, a vector is returned filled to mean
 */
VectorDouble Model::evalDrifts(const Db* db,
                               const VectorDouble& coeffs,
                               int ivar,
                               bool useSel) const
{
  VectorDouble vec;
  if (_driftList == nullptr && db != nullptr)
  {
    int nech = db->getSampleNumber(useSel);
    double mean = getMean(ivar);
    vec = VectorDouble(nech,mean);
  }
  else
  {
    vec = _driftList->evalDrifts(db, coeffs, useSel);
  }
  return vec;
}

/**
 * Sample a Model for given variable(s) and given direction
 * @param hh     Vector of distances
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param codir  Vector of direction coefficients
 * @param nostd  0 standard; +-1 corr. envelop; ITEST normalized
 * @param asCov  Produce the result as a Covariance (rather than a Variogram)
 *
 * @return The array of variogram evaluated at discretized positions
 * @return Note that its dimension is 'nh' (if 'addZero' is false and 'nh+1' otherwise)
 */
VectorDouble Model::sample(const VectorDouble& hh,
                           int ivar,
                           int jvar,
                           VectorDouble codir,
                           int nostd,
                           bool asCov)
{
  VectorDouble gg;

  if (ivar < 0 || ivar >= getVariableNumber()) return gg;
  if (jvar < 0 || jvar >= getVariableNumber()) return gg;
  int ndim = getDimensionNumber();
  if (codir.empty())
  {
    codir.resize(ndim);
    (void) GH::rotationGetDirection(ndim, 1, VectorDouble(), codir);
  }
  int nh = (int) hh.size();
  gg.resize(nh);

  model_evaluate(this, ivar, jvar, -1, 0, asCov, 0, nostd, 0, ECalcMember::LHS, nh,
                 codir, hh.data(), gg.data());
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
                             Option_VarioFit optvar,
                             Option_AutoFit mauto,
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
 * @param types       Vector of ECov
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 * @param mauto       Special parameters for Automatic fitting procedure (instance of Option_AutoFit), for exemple wmode (type of weighting function)
 * @param verbose     Verbose option
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fit(Vario *vario,
               const VectorECov &types,
               const Constraints &constraints,
               Option_VarioFit optvar,
               Option_AutoFit mauto,
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
int Model::fitFromVMap(DbGrid *dbmap,
                       const VectorECov &types,
                       const Constraints &constraints,
                       Option_VarioFit optvar,
                       Option_AutoFit mauto,
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
      cova.setRange(range);
    covs.addCov(&cova);
  }
  setCovList(&covs);

  /* Reading the drift part */

  DriftList drifts;
  for (int ibfl = 0; ret && ibfl < nbfl; ibfl++)
  {
    ret = ret && _recordRead<int>(is, "Drift Function", type);
    EDrift dtype = EDrift::fromValue(type);
    int rank_fex = 0;
    if (dtype == EDrift::F)
      ret = ret && _recordRead<int>(is, "External Drift rank", rank_fex);
    ADriftElem *drift = DriftFactory::createDriftFunc(dtype, _ctxt, rank_fex);
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
      setMean(ivar, mean);
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
    const ADriftElem *drift = getDrift(ibfl);
    ret = ret && _recordWrite<int>(os,"Drift characteristics", drift->getType().getValue());
    if (drift->getType() == EDrift::F)
      ret = ret && _recordWrite<int>(os,"External Drift rank", drift->getRankFex());
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
  delete _covaList;
  _covaList = nullptr;
  delete _driftList;
  _driftList = nullptr;
  delete _noStat;
  _noStat = nullptr;
}

void Model::_create()
{
  // TODO: The next two lines are there in order to allow direct call to
  // model::addCov() and model::addDrift
  // The defaulted types of CovAnisoList and DriftList are assumed
  _covaList = new CovLMC(_ctxt.getSpace());
  _driftList = new DriftList(_ctxt.getSpace());
}

double Model::getTotalSill(int ivar, int jvar) const
{
  return getCovAnisoList()->getTotalSill(ivar, jvar);
}

/**
 * Returns the Ball radius (from the first covariance of _covaList)
 * @return Value of the Ball Radius (if defined, i.e. for Numerical Gradient calculation)
 */
double Model::getBallRadius() const
{
  if (_covaList == nullptr) return TEST;

  // Check is performed on the first covariance

  CovAniso* cova = _covaList->getCova(0);
  double ball_radius = cova->getBallRadius();
  if (! FFFF(ball_radius)) return ball_radius;
  return 0.;
}

Model* Model::duplicate() const
{
  Model *model = nullptr;

  model = new Model(getContext());

  /* Add the list of Covariances */

  model->setCovList(getCovAnisoList());

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  /* Add non-stationarity information */

  model->addNoStat(getNoStat());

  return model;
}

/**
 * Calculate the covariance matrix between active samples of Db1
 * and active samples of Db2
 * @param covmat Returned matrix (returned as a vector).
 * @param db1 First Data Base
 * @param db2 Second Data Base (if not provided, the first Db is provided instead)
 * @param ivar0 Rank of the first variable (all variables if not defined)
 * @param jvar0 Rank of the second variable (all variables if not defined)
 * @param flag_norm 1 if the Model must be normalized beforehand
 * @param flag_cov 1 if the Model must be expressed in covariance

 *
 * @remark The returned argument must have been dimensioned beforehand to (nvar * nechA)^2 where:
 * @remark -nvar stands for the number of (active) variables
 * @remark -nechA stands for the number of active samples
 */
void Model::covMatrix(VectorDouble& covmat,
                      Db *db1,
                      Db *db2,
                      int ivar0,
                      int jvar0,
                      int flag_norm,
                      int flag_cov)
{
  model_covmat(this, db1, db2, ivar0, jvar0, flag_norm, flag_cov, covmat.data());
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
        VectorDouble gmod(npas);
        model_evaluate(this, ivar, jvar, -1, 0, 0, 0, 0, 0, ECalcMember::LHS,
                       npas, codir, hh.data(), gmod.data());

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
  else
  {
    int nclass = (int) thresholds.size();
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      if (gof < thresholds[iclass])
      {
        message(" corresponds to level #%d (1 for very good)\n", iclass+1);
        return;
      }
    }
  }
}

const EModelProperty& Model::getCovMode() const
{
  ACovAnisoList* covs;
  if (_covaList == nullptr) return EModelProperty::NONE;

  covs = dynamic_cast<CovLMCTapering*>(_covaList);
  if (covs != nullptr) return EModelProperty::TAPE;

  covs = dynamic_cast<CovLMCConvolution*>(_covaList);
  if (covs != nullptr) return EModelProperty::CONV;

  covs = dynamic_cast<CovLMCAnamorphosis*>(_covaList);
  if (covs != nullptr) return EModelProperty::ANAM;

  covs = dynamic_cast<CovLMGradient*>(_covaList);
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
  if (_covaList == nullptr) return false;
  return getCovMode() == EModelProperty::GRAD;
}

bool Model::isFlagGradientNumerical() const
{
  if (! isFlagGradient()) return false;

  // Check is performed on the first covariance
  CovGradientNumerical* cova = dynamic_cast<CovGradientNumerical*>(_covaList->getCova(0));
  if (cova != nullptr) return true;
  return false;
}

bool Model::isFlagGradientFunctional() const
{
  if (! isFlagGradient()) return false;

  // Check is performed on the first covariance
  CovGradientFunctional* cova = dynamic_cast<CovGradientFunctional*>(_covaList->getCova(0));
  if (cova != nullptr) return true;
  return false;
}

/****************************************************************************/
/*!
 **  Evaluate the drift with a given set of coefficients
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 ** \param[in]  coef    Array of coefficients
 **
 *****************************************************************************/
double Model::_evalDriftCoef(const Db* db,
                             int iech,
                             int ivar,
                             const double* coef) const
{

  VectorDouble drftab = evalDriftVec(db, iech, ECalcMember::LHS);

  /* Check if all the drift terms are defined */

  for (int il = 0; il < getDriftNumber(); il++)
    if (FFFF(drftab[il])) return TEST;

  /* Perform the correction */

  double drift = 0.;
  for (int ib = 0; ib < getDriftEquationNumber(); ib++)
  {
    double value = 0.;
    for (int il = 0; il < getDriftNumber(); il++)
      value += drftab[il] * getCoefDrift(ivar, il, ib);
    drift += value * coef[ib];
  }
  return drift;
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
  if (_covaList == nullptr)
  {
    messerr("Model is not valid: no covariance has been defined");
    return false;
  }
  if (_covaList->getCovNumber() <= 0)
  {
    messerr("Model is not valid: no covariance has been defined");
    return false;
  }

  // Drifts: there should be valid
  if (_driftList != nullptr)
  {
    if (! _driftList->isValid()) return false;
  }
  return true;
}
