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
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/Vector.hpp"
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
#include "Model/NoStatArray.hpp"
#include "Model/EModelProperty.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include <math.h>

Model::Model(const CovContext &ctxt)
    :
    AStringable(),
    ASerializable(),
    _covaList(nullptr),
    _driftList(nullptr),
    _noStat(nullptr),
    _ctxt(ctxt)
{
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
    _driftList = dynamic_cast<DriftList*>(m._driftList->clone());
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
      _driftList = dynamic_cast<DriftList*>(m._driftList->clone());
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
  _ctxt = CovContext(db);
  _create();
  return 0;
}

Model* Model::create(const CovContext& ctxt)
{
  return new Model(ctxt);
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

int Model::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "Model", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

Model* Model::createFromNF(const String &neutralFilename, bool verbose)
{
  Model* model = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "Model", is, verbose))
  {
    model = new Model();
    if (model->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem when reading the Neutral File");
      delete model;
      model = nullptr;
    }
    is.close();
  }
  return model;
}

String Model::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ncov   = getCovaNumber();
  int ndrift = getDriftNumber();

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

  if (isNoStat()) sstr << _noStat->toString();

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
  // TODO: the type of CovAnisoList is defaulted to CovLMC
  if (_covaList == nullptr) _covaList = new CovLMC();
  _covaList->addCov(cov);
}

/**
 * Add a list of Drifts. This operation cleans any previously stored drift function
 * @param driftlist List of Drifts to be added
 */
void Model::setDriftList(const DriftList* driftlist)
{
  if (driftlist == nullptr) return;
  if (_driftList != nullptr) delete _driftList;
  _driftList = dynamic_cast<DriftList*>(driftlist->clone());
}

void Model::addDrift(const ADriftElem *drift)
{
  if (drift == nullptr) return;
  if (_driftList == nullptr)
    _driftList = new DriftList();
  _driftList->addDrift(drift);
}

void Model::addDrift(const VectorString &driftSymbols)
{
  if (_driftList == nullptr)
    my_throw("Model::addDrift: cannot add an element to non-initialized _driftList");

  for (int i = 0; i < (int) driftSymbols.size(); i++)
  {
    int rank = 0;
    EDrift type = DriftFactory::identifyDrift(driftSymbols[i], &rank, _ctxt);
    ADriftElem *drift = DriftFactory::createDriftFunc(type, _ctxt);
    drift->setRankFex(rank);
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
  _driftList->delAllDrift();
}

const CovAniso* Model::getCova(unsigned int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getCova(icov);
}
CovAniso* Model::getCova(unsigned int icov)
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getCova(icov);
}
int Model::getCovaNumber() const
{
  if (_covaList == nullptr) return 0;
  return _covaList->getCovNumber();
}
const ECov& Model::getCovaType(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
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
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getSill(icov, ivar, jvar);
}
double Model::getParam(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getParam(icov);
}
bool Model::isCovaFiltered(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->isFiltered(icov);
}
String Model::getCovName(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getCovName(icov);
}
int Model::getGradParamNumber(int icov) const
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  return _covaList->getGradParamNumber(icov);
}
void Model::setSill(int icov, int ivar, int jvar, double value)
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
  _covaList->setSill(icov, ivar, jvar, value);
}
void Model::setCovaFiltered(int icov, bool filtered)
{
  if (_covaList == nullptr)
    my_throw("Covariance List is empty");
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
  if (!isNoStat())
    my_throw("Nostat is not defined");
  return _noStat->getICov(ipar);
}
const EConsElem& Model::getNoStatElemType(int ipar)
{
  if (!isNoStat())
    my_throw("Nostat is not defined");
  return _noStat->getType(ipar);
}
const DriftList* Model::getDriftList() const
{
  return _driftList;
}
const ADriftElem* Model::getDrift(int il) const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->getDrift(il);
}
ADriftElem* Model::getDrift(int il)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
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
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->getType(il);
}
int Model::getRankFext(int il) const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
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
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->getCoefDrift(ivar, il, ib);
}
int Model::getDriftEquationNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getDriftEquationNumber();
}
bool Model::isDriftFiltered(unsigned int il) const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->isFiltered(il);
}
bool Model::isDriftDefined(const EDrift& type0) const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->isDriftDefined(type0);
}
bool Model::isDriftDifferentDefined(const EDrift& type0) const
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->isDriftDifferentDefined(type0);
}
void Model::setCoefDrift(int ivar, int il, int ib, double coeff)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  _driftList->setCoefDrift(ivar, il, ib, coeff);
}
void Model::setCoefDriftByRank(int rank, double coeff)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  _driftList->setCoefDriftByRank(rank, coeff);
}
void Model::setDriftFiltered(int il, bool filtered)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  _driftList->setFiltered(il, filtered);
}
VectorDouble Model::getDrift(const Db *db, int ib, bool useSel)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
  return _driftList->getDrift(db, ib, useSel);
}
VectorVectorDouble Model::getDrifts(const Db *db, bool useSel)
{
  if (_driftList == nullptr)
    my_throw("Drift List if empty");
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
    if (_driftList == nullptr)
      my_throw("Drift List if empty");
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
 * @param hmax   Maximum distance to be sampled
 * @param nh     Number of discretization steps
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param codir  Vector of direction coefficients
 * @param nostd  0 standard; +-1 corr. envelop; ITEST normalized
 * @param addZero Add the zero distance location
 *
 * @return The array of variogram evaluated at discretized positions
 * @return Note that its dimension is 'nh' (if 'addZero' is false and 'nh+1' otherwise)
 */
VectorDouble Model::sample(double hmax,
                           int nh,
                           int ivar,
                           int jvar,
                           VectorDouble codir,
                           int nostd,
                           bool addZero)
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
  int nhloc = (addZero) ? nh + 1: nh;
  hh.resize(nhloc);
  gg.resize(nhloc);

  int ecr = 0;
  if (addZero)
  {
    hh[ecr] = 0.;
    gg[ecr] = 0.;
  }
  for (int i = 0; i < nh; i++)
    hh[ecr++] = hmax * (i+1) / nh;

  model_evaluate(this, ivar, jvar, -1, 0, 0, 0, nostd, 0, ECalcMember::LHS, nhloc,
                 codir, hh.data(), gg.data());
  return gg;
}

/**
 * Automatic Fitting procedure
 *
 * @param vario       Experimental variogram to be fitted
 * @param types       Vector of ECov integer values (treated as 'int' due to compiler problem)
 * @param verbose     Verbose option
 * @param mauto       Special parameters for Automatic fitting procedure
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fitFromCovIndices(Vario *vario,
                             const VectorInt &types,
                             bool verbose,
                             Option_AutoFit mauto,
                             const Constraints &constraints,
                             Option_VarioFit optvar)
{
  if (vario == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCovas();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  for (int is = 0; is < (int) types.size(); is++)
  {
    ECov covtype = ECov::fromValue(types[is]);
    CovAniso cov = CovAniso(covtype, _ctxt);
    addCov(&cov);
  }

  return model_auto_fit(vario, this, verbose, mauto, constraints, optvar);
}

/**
 * Automatic Fitting procedure
 *
 * @param vario       Experimental variogram to be fitted
 * @param types       Vector of ECov
 * @param verbose     Verbose option
 * @param mauto       Special parameters for Automatic fitting procedure
 * @param constraints Set of Constraints
 * @param optvar      Set of options
 *
 * @return 0 if no error, 1 otherwise
 */
int Model::fit(Vario *vario,
               const std::vector<ECov> &types,
               bool verbose,
               Option_AutoFit mauto,
               const Constraints &constraints,
               Option_VarioFit optvar)
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

int Model::_deserialize(std::istream& is, bool /*verbose*/)
{
  double field, range, value, param;
  int ndim, nvar, ncova, nbfl, type, flag_aniso, flag_rotation;
  VectorDouble aniso_ranges, aniso_rotmat;

  // Delete previous Model contents (if any)
  _clear();

  /* Create the Model structure */

  bool ret = _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Number of Variables", nvar);
  ret = ret && _recordRead<double>(is, "Field dimension", field);
  ret = ret && _recordRead<int>(is, "Number of Basic Structures", ncova);
  ret = ret && _recordRead<int>(is, "Number of Basic Drift Functions", nbfl);
  if (! ret) return 1;

  /// TODO : Force SpaceRN creation (deserialization doesn't know yet how to manage other space types)
  _ctxt = CovContext(nvar, ndim, field);
  _create();

  /* Reading the covariance part */

  CovLMC covs(_ctxt.getSpace());
  for (int icova = 0; icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    ret = ret && _recordRead<int>(is, "Covariance Type", type);
    ret = ret && _recordRead<double>(is, "Isotropic Range", range);
    ret = ret && _recordRead<double>(is, "Model third Parameter", param);
    ret = ret && _recordRead(is, "Flag for Anisotropy", flag_aniso);
    if (! ret) return 1;
    if (flag_aniso)
    {
      aniso_ranges.resize(ndim);
      // In fact, the file contains the anisotropy coefficients
      // After reading, we must turn them into anisotropic ranges
      for (int idim = 0; idim < ndim; idim++)
        ret = ret && _recordRead<double>(is, "Anisotropy coefficient", aniso_ranges[idim]);
      if (! ret) return 1;
      for (int idim = 0; idim < ndim; idim++)
        aniso_ranges[idim] *= range;

      ret = ret && _recordRead<int>(is, "Flag for Anisotropy Rotation", flag_rotation);
      if (! ret) return 1;
      if (flag_rotation)
      {
        // Warning: the storage in the File is performed by column
        // whereas the internal storage is by column (TODO : ???)
        aniso_rotmat.resize(ndim * ndim);
        int lec = 0;
        for (int idim = 0; idim < ndim; idim++)
          for (int jdim = 0; jdim < ndim; jdim++)
            ret = ret && _recordRead<double>(is, "Anisotropy Rotation Matrix", aniso_rotmat[lec++]);
      }
    }
    if (! ret) return 1;

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
  for (int ibfl = 0; ibfl < nbfl; ibfl++)
  {
    ret = ret && _recordRead<int>(is, "Drift Function", type);
    EDrift dtype = EDrift::fromValue(type);
    ADriftElem *drift = DriftFactory::createDriftFunc(dtype, _ctxt);
    drift->setRankFex(0); // TODO : zero? really?
    drifts.addDrift(drift);
    delete drift;
  }
  setDriftList(&drifts);

  /* Reading the matrix of means (only if nbfl <= 0) */

  if (nbfl <= 0) for (int ivar = 0; ivar < nvar; ivar++)
  {
    double mean;
    ret = ret && _recordRead<double>(is, "Mean of Variable", mean);
    setMean(ivar, mean);
  }

  /* Reading the matrices of sills (optional) */

  for (int icova = 0; icova < ncova; icova++)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        if (_recordRead<double>(is, "Matrix of Sills", value))
          setSill(icova, ivar, jvar, value);
      }
  }

  /* Reading the variance-covariance at the origin (optional) */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      if (_recordRead<double>(is, "Variance-covariance at Origin", value))
        setCovar0(ivar, jvar, value);
    }

  return ret ? 0 : 1;
}

int Model::_serialize(std::ostream& os, bool /*verbose*/) const
{
  /* Write the Model structure */

  bool ret = _recordWrite<int>(os, "", getDimensionNumber());
  ret = ret && _recordWrite<int>(os, "", getVariableNumber());
  ret = ret && _recordWrite<double>(os, "General parameters", getField());
  ret = ret && _recordWrite<int>(os, "Number of basic covariance terms", getCovaNumber());
  ret = ret && _recordWrite<int>(os, "Number of drift terms", getDriftNumber());

  /* Writing the covariance part */

  for (int icova = 0; icova < getCovaNumber(); icova++)
  {
    const CovAniso *cova = getCova(icova);
    ret = ret && _recordWrite<int>(os, "", cova->getType().getValue());
    ret = ret && _recordWrite<double>(os, "", cova->getRange());
    ret = ret && _recordWrite<double>(os, "Covariance characteristics", cova->getParam());

    // Writing the Anisotropy information

    ret = ret && _recordWrite<int>(os, "Anisotropy Flag", cova->getFlagAniso());

    if (!cova->getFlagAniso()) continue;
    for (int idim = 0; idim < getDimensionNumber(); idim++)
      ret = ret && _recordWrite<double>(os, "", cova->getAnisoCoeffs(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<int>(os, "Anisotropy Rotation Flag", cova->getFlagRotation());

    if (!cova->getFlagRotation()) continue;
    // Storing the rotation matrix by Column (compatibility)
    for (int idim = 0; idim < getDimensionNumber(); idim++)
      for (int jdim = 0; jdim < getDimensionNumber(); jdim++)
        ret = ret && _recordWrite<double>(os, "", cova->getAnisoRotMat(jdim, idim));
    _commentWrite(os, "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (int ibfl = 0; ibfl < getDriftNumber(); ibfl++)
  {
    const ADriftElem *drift = getDrift(ibfl);
    ret = ret && _recordWrite<int>(os,"Drift characteristics", drift->getType().getValue());
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (getDriftNumber() <= 0)
    for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    {
      ret = ret && _recordWrite<double>(os, "Mean of Variables", getContext().getMean(ivar));
    }

  /* Writing the matrices of sills (optional) */

  for (int icova = 0; icova < getCovaNumber(); icova++)
  {
    for (int ivar = 0; ivar < getVariableNumber(); ivar++)
      for (int jvar = 0; jvar < getVariableNumber(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getSill(icova, ivar, jvar));
    _commentWrite(os, "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (int ivar = 0; ivar < getVariableNumber(); ivar++)
    for (int jvar = 0; jvar < getVariableNumber(); jvar++)
      ret = ret && _recordWrite<double>(os, "", getContext().getCovar0(ivar, jvar));
  _commentWrite(os, "Var-Covar at origin");

  return 0;
}

void Model::_clear()
{
  delete _covaList;
  delete _driftList;
  delete _noStat;
}
void Model::_create()
{
  // TODO: The next two lines are there in order to allow direct call to
  // model::addCova() and model::addDrift
  // The defaulted types of CovAnisoList and DriftList are assumed
  _covaList = new CovLMC();
  _driftList = new DriftList();
}

double Model::getTotalSill(int ivar, int jvar) const
{
  double var = 0.;
  for (int icov=0; icov<getCovaNumber(); icov++)
    var += getSill(icov,ivar,jvar);
  return var;
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
 * scaled to the sill
 * @param vario Experimental variogram
 * @return Value for the Goodness-of_fit (as percentage of the total sill)
 */
double Model::gofToVario(const Vario* vario)
{
  int nvar = getVariableNumber();
  int ndir = vario->getDirectionNumber();

  double total = 0.;

  // Loop on the pair of variables

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double sill = getTotalSill(ivar, jvar);

      // Loop on the variogram directions

      double totdir = 0.;
      for (int idir = 0; idir < ndir; idir++)
      {

        // Read information from Experimental Variogram

        VectorDouble codir = vario->getCodir(idir);
        VectorDouble hh = vario->getHhVec(idir, ivar, jvar);
        VectorDouble gexp = vario->getGgVec(idir, ivar, jvar);

        // Evaluate the Model

        int npas = (int) gexp.size();
        VectorDouble gmod(npas);
        model_evaluate(this, ivar, jvar, -1, 0, 0, 0, 0, 0, ECalcMember::LHS,
                       npas, codir, hh.data(), gmod.data());

        // Evaluate the score

        double totpas = 0;
        for (int ipas = 0; ipas < npas; ipas++)
        {
          double ecart = gexp[ipas] - gmod[ipas];
          totpas += ecart * ecart;
        }
        totpas  = sqrt(totpas) / (double) npas;
        totdir += totpas;
      }
      totdir /= (double) ndir;
      totdir /= sill;
      total  += ABS(totdir);
    }
  total = 100. * total / (double) (nvar * nvar);
  return total;
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
