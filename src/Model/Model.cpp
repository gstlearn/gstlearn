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
#include "Model/Model.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Db/Db.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftList.hpp"
#include "Enum/ECov.hpp"
#include "Model/CovInternal.hpp"
#include "Model/ModelCovList.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_define.h"
#include "geoslib_f.h"
#include <cmath>

namespace gstlrn
{
Model::Model(const CovContext& ctxt)
  : AStringable()
  , ASerializable()
  , ModelCovList(ctxt)
{
  _create();
}

Model::Model(Id nvar, Id ndim)
  : AStringable()
  , ASerializable()
  , ModelCovList()
{
  auto space = SpaceRN::create(ndim);
  _ctxt      = CovContext(nvar, space);
  _create();
}

Model::Model(const Model& m)
  : AStringable(m)
  , ASerializable(m)
  , ModelCovList(m)
{
}

Model& Model::operator=(const Model& m)
{
  if (this != &m)
  {
    ModelCovList::operator=(m);
    AStringable::operator=(m);
    ASerializable::operator=(m);
    setCovAnisoList(dynamic_cast<CovAnisoList*>(m.getCovAnisoList()->clone()));
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

Id Model::resetFromDb(const Db* db)
{
  Id ndim = db->getNDim();
  Id nvar = db->getNLoc(ELoc::Z);
  if (nvar <= 0) nvar = 1;
  auto space = SpaceRN::create(ndim);
  _ctxt      = CovContext(nvar, space);
  _create();
  return 0;
}

bool Model::_isValid() const
{
  // Covariances: there should be some defined
  if (getCovAnisoList() == nullptr)
  {
    messerr("Model is not valid: no covariance has been defined");
    return false;
  }

  // Drifts: there should be valid
  if (_driftList != nullptr)
  {
    if (!_driftList->isValid()) return false;
  }

  // Check the consistency between the Covariance and the Drift parts
  Id irf_drift = getDriftMaxIRFOrder();
  Id irf_cova  = getCovMinIRFOrder();
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

Model* Model::create(const CovContext& ctxt)
{
  return new Model(ctxt);
}

Model* Model::createFromEnvironment(Id nvar, Id ndim)
{
  return new Model(nvar, ndim);
}

Model* Model::createNugget(Id nvar, Id ndim, double sill)
{
  auto* model = new Model(nvar, ndim);
  model->addCovFromParam(ECov::NUGGET, 0., sill);
  return model;
}

Model* Model::createFromParam(const ECov& type,
                              double range,
                              double sill,
                              double param,
                              const VectorDouble& ranges,
                              const MatrixSymmetric& sills,
                              const VectorDouble& angles,
                              const ASpaceSharedPtr& space,
                              bool flagRange)
{
  Id nvar = 1;
  if (!sills.empty()) nvar = sills.getNRows();

  auto spaceloc = ASpace::getDefaultSpaceIfNull(space);

  if (!ranges.empty())
  {
    Id ndim       = static_cast<Id>(spaceloc->getNDim());
    Id ndimRanges = static_cast<Id>(ranges.size());
    if (ndimRanges != 1 && ndimRanges != ndim)
    {
      messerr("Incompatibility between:");
      messerr("Space Dimension = %d", ndim);
      messerr("Dimension of argument 'ranges' = %d", ndimRanges);
      return nullptr;
    }
  }

  CovContext ctxt(nvar, space);
  auto* model = new Model(ctxt);
  model->addCovFromParam(type, range, sill, param, ranges, sills, angles,
                         flagRange);

  return model;
}

Model* Model::createFromParamOldStyle(const ECov& type,
                                      double range,
                                      double sill,
                                      double param,
                                      const VectorDouble& ranges,
                                      const VectorDouble& sills,
                                      const VectorDouble& angles,
                                      const ASpaceSharedPtr& space,
                                      bool flagRange)
{
  Id nvar = 1;
  if (!sills.empty())
    nvar = static_cast<Id>(sqrt(sills.size()));

  auto spaceloc = ASpace::getDefaultSpaceIfNull(space);

  if (!ranges.empty())
  {
    Id ndim       = static_cast<Id>(spaceloc->getNDim());
    Id ndimRanges = static_cast<Id>(ranges.size());
    if (ndimRanges != 1 && ndimRanges != ndim)
    {
      messerr("Incompatibility between:");
      messerr("Space Dimension = %d", ndim);
      messerr("Dimension of argument 'ranges' = %d", ndimRanges);
      return nullptr;
    }
  }

  CovContext ctxt(nvar, spaceloc);
  auto* model = new Model(ctxt);
  model->addCovFromParamOldStyle(type, range, sill, param, ranges, sills,
                                 angles, flagRange);

  return model;
}

Model* Model::createFromDb(const Db* db)
{
  auto* model = new Model();
  if (model->resetFromDb(db) != 0)
  {
    messerr("Problem when creating Model from Db");
    delete model;
    return nullptr;
  }
  return model;
}

Model* Model::createFromNF(const String& NFFilename, bool verbose)
{
  auto* model = new Model();
  if (model->_fileOpenAndDeserialize(NFFilename, verbose)) return model;
  delete model;
  return nullptr;
}

Model* Model::createFromVario(Vario* vario,
                              const VectorECov& types,
                              const Constraints& constraints,
                              const Option_VarioFit& optvar,
                              const Option_AutoFit& mauto,
                              bool verbose)
{
  auto* model = new Model();
  if (model->fit(vario, types, constraints, optvar, mauto, verbose) != 0)
  {
    messerr("Problem when creating Model from fitting an Experimental variogram");
    delete model;
    return nullptr;
  }
  return model;
}

String Model::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  auto ncov   = getNCov();
  auto ndrift = getNDrift();
  if (ncov <= 0 && ndrift <= 0) return sstr.str();

  sstr << toTitle(0, "Model characteristics");
  sstr << "Space dimension              = " << getNDim()
       << std::endl;
  sstr << "Number of variable(s)        = " << getNVar() << std::endl;
  sstr << "Number of basic structure(s) = " << ncov << std::endl;
  sstr << "Number of drift function(s)  = " << ndrift << std::endl;
  sstr << "Number of drift equation(s)  = " << getNDriftEquation() << std::endl;

  /* Covariance part */

  if (ncov > 0)
  {
    sstr << toTitle(1, "Covariance Part");
    sstr << getCovAnisoList()->toString();
  }

  /* Drift part */

  if (ndrift > 0)
    sstr << toTitle(1, "Drift Part");

  sstr << _driftList->toString();

  if (isFlagLinked())
    sstr << "Drifts are linked" << std::endl;

  return sstr.str();
}

/**
 * Add a list of Covariances. This operation cleans any previously stored covariance
 * @param covalist List of Covariances to be added
 */
void Model::setCovAnisoList(const CovAnisoList* covalist)
{
  if (covalist == nullptr)
  {
    messerr("Warning, the covariance is nullptr.");
    return;
  }
  ModelCovList::setCovList(covalist);
}

void Model::addCovFromParamOldStyle(const ECov& type,
                                    double range,
                                    double sill,
                                    double param,
                                    const VectorDouble& ranges,
                                    const VectorDouble& sills,
                                    const VectorDouble& angles,
                                    bool flagRange)
{
  // Check consistency with parameters of the model

  auto ndim = getNDim();
  if (!ranges.empty())
  {
    if (ndim > 0 && ranges.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)", static_cast<Id>(ranges.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = ranges.size();
  }
  if (!angles.empty())
  {
    if (ndim > 0 && angles.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)", static_cast<Id>(angles.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = angles.size();
  }
  auto nvar = getNVar();
  if (!sills.empty())
  {
    if (nvar > 0 && static_cast<Id>(sills.size()) != nvar * nvar)
    {
      messerr("Mismatch between the size of 'sills' (%d)", static_cast<Id>(sills.size()));
      messerr("and the Number of variables stored in the Model (%d)", nvar);
      messerr("Operation is cancelled");
      return;
    }
    nvar = static_cast<Id>(sqrt(static_cast<double>(sills.size())));
  }

  // Define the covariance

  auto space = SpaceRN::create(static_cast<Id>(ndim));
  _ctxt      = CovContext(nvar, space);
  CovAniso cov(type, _ctxt);

  // Define the Third parameter
  double parmax = cov.getParMax();
  if (param > parmax) param = parmax;
  cov.setParam(param);

  // Define the range
  if (!ranges.empty())
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

  // Define the sill
  if (!sills.empty())
    cov.setSill(sills);
  else
  {
    if (nvar <= 1)
      cov.setSill(sill);
    else
    {
      MatrixSymmetric locsills(nvar);
      locsills.setIdentity(sill);
      cov.setSill(locsills);
    }
  }

  if (!angles.empty())
    cov.setAnisoAngles(angles);
  addCov(cov);
}

void Model::addCovFromParam(const ECov& type,
                            double range,
                            double sill,
                            double param,
                            const VectorDouble& ranges,
                            const MatrixSymmetric& sills,
                            const VectorDouble& angles,
                            bool flagRange)
{
  // Check consistency with parameters of the model

  auto ndim = getNDim();
  if (!ranges.empty())
  {
    if (ndim > 0 && ranges.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)",
              static_cast<Id>(ranges.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = ranges.size();
  }
  if (!angles.empty())
  {
    if (ndim > 0 && angles.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)",
              static_cast<Id>(angles.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = angles.size();
  }
  auto nvar = getNVar();
  if (!sills.empty())
  {
    if (nvar > 0 && nvar != sills.getNCols())
    {
      messerr("Mismatch between the number of rows 'sills' (%d)", sills.getNRows());
      messerr("and the Number of variables stored in the Model (%d)", nvar);
      messerr("Operation is cancelled");
      return;
    }
    nvar = static_cast<Id>(sqrt(static_cast<double>(sills.size())));
  }

  // Define the covariance

  auto space = SpaceRN::create(static_cast<Id>(ndim));
  _ctxt      = CovContext(nvar, space);
  CovAniso cov(type, _ctxt);

  // Define the Third parameter
  double parmax = cov.getParMax();
  if (param > parmax) param = parmax;
  cov.setParam(param);

  // Define the range
  if (!ranges.empty())
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

  // Define the sill
  if (!sills.empty())
    cov.setSill(sills);
  else
  {
    if (nvar <= 1)
      cov.setSill(sill);
    else
    {
      MatrixSymmetric locsills(nvar);
      locsills.setIdentity(sill);
      cov.setSill(locsills);
    }
  }

  _ctxt.setNVar(cov.getNVar());
  _copyCovContext();
  if (!angles.empty()) cov.setAnisoAngles(angles);
  addCovAniso(cov);
}

double Model::evalCovFromIncr(const VectorDouble& incr,
                              Id icov,
                              const ECalcMember& member) const
{
  if (_cova == nullptr) return TEST;
  const CovAnisoList* covalist = castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;

  if (member != ECalcMember::LHS && covalist->isFiltered(icov))
    return (0.);
  return getCovAniso(icov)->evalIvarIpas(1., incr);
}

/**
 * Defining an Anamorphosis information for the Model
 * (in fact, this is added to CovAnisoList part and transforms it from CovAnisoList to CovLMCAnamorphosis
 * @param anam Pointer to the anamorphosis
 * @param strcnt Array of covariance description used for IR case
 * @return
 */
Id Model::setAnam(const gstlrn::AAnam* anam, const VectorInt& strcnt)
{
  if (anam == nullptr)
  {
    messerr("You must define 'anam' beforehand");
    return 1;
  }
  if (hasAnam())
  {
    // CovAnisoList is already a covLMCAnamorphosis, simply update the anamorphosis
    auto* cov = dynamic_cast<CovLMCAnamorphosis*>(_cova.get());
    if (cov == nullptr)
    {
      messerr("Impossible to reach the internal CovLMCAnamorphosis structure");
      return 1;
    }
    cov->setAnam(anam);
  }
  else
  {
    CovAnisoList* cov = dynamic_cast<CovAnisoList*>(getCovAnisoListModify());
    if (cov == nullptr)
    {
      messerr("Impossible to add 'anam' to the covariance part of the Model");
      messerr("The original covariance is probably not a 'CovAnisoList'");
      messerr("The original covariance is probably not a 'CovAnisoList'");
      return 1;
    }

    // Initiate a new CovLMCAnamorphosis class
    auto* newcov = new CovLMCAnamorphosis(*cov, anam, strcnt);

    // Replace the current list by the newly create one (CovLMCAnamorphosis)
    ModelCovList::setCovList(newcov);
  }
  return 0;
}

Id Model::unsetAnam()
{
  if (!hasAnam())
  {
    // CovAnisoList does not have any Anam: do nothing
    return 0;
  }
  auto* cov = dynamic_cast<CovAnisoList*>(_cova.get());
  if (cov == nullptr)
  {
    messerr("Impossible to unset 'anam' from the covariance part of the Model");
    messerr("The original covariance is probably not valid");
    return 1;
  }

  // Initiate a new CovAnisoList class
  auto* newcov = new CovAnisoList(*cov);

  // Replace the current list by the newly create one (CovLMCAnamorphosis)

  ModelCovList::setCovList(newcov);

  return 0;
}

void Model::_copyCovContext()
{
  if (_cova == nullptr) return;
  CovAnisoList* covalist = _castInCovAnisoList();
  if (covalist != nullptr) covalist->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
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
Id Model::fitFromCovIndices(Vario* vario,
                            const VectorECov& types,
                            const Constraints& constraints,
                            const Option_VarioFit& optvar,
                            const Option_AutoFit& mauto,
                            bool verbose)
{
  if (vario == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCov();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  for (Id is = 0; is < static_cast<Id>(types.size()); is++)
  {
    CovAniso cov(types[is], _ctxt);
    addCov(cov);
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
Id Model::fit(Vario* vario,
              const VectorECov& types,
              const Constraints& constraints,
              const Option_VarioFit& optvar,
              const Option_AutoFit& mauto,
              bool verbose)
{
  if (vario == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCov();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  _driftList->copyCovContext(_ctxt);

  for (Id is = 0; is < static_cast<Id>(types.size()); is++)
  {
    CovAniso cov(types[is], _ctxt);
    addCov(cov);
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
Id Model::fitFromVMap(DbGrid* dbmap,
                      const VectorECov& types,
                      const Constraints& constraints,
                      const Option_VarioFit& optvar,
                      const Option_AutoFit& mauto,
                      bool verbose)
{
  if (dbmap == nullptr) return 1;

  // Clean out possible covariances in the existing model

  delAllCov();

  // Add the relevant covariances

  for (Id is = 0; is < static_cast<Id>(types.size()); is++)
  {
    CovAniso cov(types[is], _ctxt);
    addCov(cov);
  }
  return vmap_auto_fit(dbmap, this, verbose, mauto, constraints, optvar);
}

bool Model::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  Id ndim          = 0;
  Id nvar          = 0;
  Id ncova         = 0;
  Id nbfl          = 0;
  Id type          = 0;
  Id flag_aniso    = 0;
  Id flag_rotation = 0;

  double field = 0.;
  double range = 0.;
  double param = 0.;
  double value = 0.;

  VectorDouble aniso_ranges;
  VectorDouble aniso_rotmat;

  /* Create the Model structure */

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Space Dimension", ndim);
  ret      = ret && _recordRead<Id>(is, "Number of Variables", nvar);
  ret      = ret && _recordRead<double>(is, "Field dimension", field);
  ret      = ret && _recordRead<Id>(is, "Number of Basic Structures", ncova);
  ret      = ret && _recordRead<Id>(is, "Number of Basic Drift Functions", nbfl);
  if (!ret) return ret;

  /// TODO : Force SpaceRN creation (deserialization doesn't know yet how to manage other space types)
  _ctxt = CovContext(nvar, ndim);
  _ctxt.setField(field);
  _clear();
  _create();

  /* Reading the covariance part and store it into a CovAnisoList */

  CovAnisoList covs(_ctxt);
  for (Id icova = 0; ret && icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    ret = ret && _recordRead<Id>(is, "Covariance Type", type);
    ret = ret && _recordRead<double>(is, "Isotropic Range", range);
    ret = ret && _recordRead<double>(is, "Model third Parameter", param);
    ret = ret && _recordRead(is, "Flag for Anisotropy", flag_aniso);
    if (!ret) return ret;
    if (flag_aniso != 0)
    {
      aniso_ranges.resize(ndim);
      // In fact, the file contains the anisotropy coefficients
      // After reading, we must turn them into anisotropic ranges
      for (Id idim = 0; idim < ndim; idim++)
        ret = ret && _recordRead<double>(is, "Anisotropy coefficient", aniso_ranges[idim]);
      if (!ret) return ret;
      for (Id idim = 0; idim < ndim; idim++)
        aniso_ranges[idim] *= range;

      ret = ret && _recordRead<Id>(is, "Flag for Anisotropy Rotation", flag_rotation);
      if (!ret) return ret;
      if (flag_rotation != 0)
      {
        // Warning: the storage in the File is performed by column
        // whereas the internal storage is by column (TODO : ???)
        aniso_rotmat.resize(ndim * ndim);
        Id lec = 0;
        for (Id idim = 0; ret && idim < ndim; idim++)
          for (Id jdim = 0; ret && jdim < ndim; jdim++)
            ret = ret && _recordRead<double>(is, "Anisotropy Rotation Matrix", aniso_rotmat[lec++]);
      }
    }
    if (!ret) return ret;

    CovAniso cova(ECov::fromValue(type), _ctxt);
    cova.setParam(param);
    if (flag_aniso != 0)
    {
      cova.setRanges(aniso_ranges);
      if (flag_rotation != 0) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRangeIsotropic(range);
    covs.addCov(cova);
  }
  setCovAnisoList(&covs);

  /* Reading the drift part */

  DriftList drifts(_ctxt);
  ADrift* drift;
  for (Id ibfl = 0; ret && ibfl < nbfl; ibfl++)
  {
    ret = true; // Reset 'ret' to continue reading after previous error...
    String driftname;
    ret   = ret && _recordRead<String>(is, "Drift Identifier", driftname);
    drift = DriftFactory::createDriftByIdentifier(driftname);
    drifts.addDrift(drift);
    delete drift;
  }
  setDriftList(&drifts);

  /* Reading the matrix of means (only if nbfl <= 0) */

  if (nbfl <= 0)
    for (Id ivar = 0; ret && ivar < nvar; ivar++)
    {
      double mean = 0.;
      ret         = ret && _recordRead<double>(is, "Mean of Variable", mean);
      setMean(mean, ivar);
    }

  /* Reading the matrices of sills (optional) */

  for (Id icova = 0; icova < ncova && ret; icova++)
  {
    for (Id ivar = 0; ret && ivar < nvar; ivar++)
      for (Id jvar = 0; ret && jvar < nvar; jvar++)
      {
        ret = ret && _recordRead<double>(is, "Matrix of Sills", value);
        if (ret) setSill(icova, ivar, jvar, value);
      }
  }

  /* Reading the variance-covariance at the origin (optional) */

  for (Id ivar = 0; ret && ivar < nvar; ivar++)
    for (Id jvar = 0; ret && jvar < nvar; jvar++)
    {
      ret = ret && _recordRead<double>(is, "Variance-covariance at Origin",
                                       value);
      if (ret) setCovar0(ivar, jvar, value);
    }

  return ret;
}

bool Model::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;

  /* Write the Model structure */

  ret = ret && _recordWrite<Id>(os, "", static_cast<Id>(getNDim()));
  ret = ret && _recordWrite<Id>(os, "", getNVar());
  ret = ret && _recordWrite<double>(os, "General parameters", getField());
  ret = ret && _recordWrite<Id>(os, "Number of basic covariance terms", getNCov());
  ret = ret && _recordWrite<Id>(os, "Number of drift terms", getNDrift());

  /* Writing the covariance part */

  for (Id icova = 0; ret && icova < getNCov(); icova++)
  {
    const CovAniso* cova = getCovAniso(icova);
    ret                  = ret && _recordWrite<Id>(os, "", cova->getType().getValue());
    ret                  = ret && _recordWrite<double>(os, "", cova->getRangeIso());
    ret                  = ret && _recordWrite<double>(os, "Covariance characteristics", cova->getParam());

    // Writing the Anisotropy information

    ret = ret && _recordWrite<Id>(os, "Anisotropy Flag", static_cast<Id>(cova->getFlagAniso()));

    if (!cova->getFlagAniso()) continue;

    for (Id idim = 0; ret && idim < static_cast<Id>(getNDim()); idim++)
      ret = ret && _recordWrite<double>(os, "", cova->getAnisoCoeff(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<Id>(os, "Anisotropy Rotation Flag", static_cast<Id>(cova->getFlagRotation()));

    if (!cova->getFlagRotation()) continue;

    // Storing the rotation matrix by Column (compatibility)
    Id ndim = static_cast<Id>(getNDim());
    for (Id idim = 0; ret && idim < ndim; idim++)
      for (Id jdim = 0; ret && jdim < ndim; jdim++)
        ret = ret && _recordWrite<double>(os, "", cova->getAnisoRotMatElement(jdim, idim));
    ret = ret && _commentWrite(os, "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (Id ibfl = 0; ret && ibfl < getNDrift(); ibfl++)
  {
    const ADrift* drift = getDrift(ibfl);
    ret                 = ret && _recordWrite<String>(os, "Drift Identifier", drift->getDriftName());
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (getNDrift() <= 0)
    for (Id ivar = 0; ret && ivar < getNVar(); ivar++)
    {
      ret = ret && _recordWrite<double>(os, "Mean of Variables", getMean(ivar));
    }

  /* Writing the matrices of sills (optional) */

  for (Id icova = 0; ret && icova < getNCov(); icova++)
  {
    for (Id ivar = 0; ret && ivar < getNVar(); ivar++)
      for (Id jvar = 0; ret && jvar < getNVar(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getSill(icova, ivar, jvar));
    ret = ret && _commentWrite(os, "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (Id ivar = 0; ret && ivar < getNVar(); ivar++)
    for (Id jvar = 0; ret && jvar < getNVar(); jvar++)
      ret = ret && _recordWrite<double>(os, "", getContext()->getCovar0(ivar, jvar));
  ret = ret && _commentWrite(os, "Var-Covar at origin");

  return ret;
}

void Model::_clear()
{
  _cova = nullptr;
  delete _driftList;
  _driftList = nullptr;
}

void Model::_create()
{
  // TODO: The next two lines are there in order to allow direct call to
  // model::addCov() and model::addDrift
  // The defaulted types of CovAnisoList and DriftList are assumed

  CovAnisoList tmp(_ctxt);
  delete _driftList;
  _driftList = new DriftList(_ctxt);
  setCovAnisoList(&tmp);
}

void Model::addCovAniso(const CovAniso& cov)
{
  ModelCovList::addCov(cov);
}

Model* Model::duplicate() const
{
  auto* model = new Model(*getContext());

  /* Add the list of Covariances */

  model->setCovAnisoList(getCovAnisoList());

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  return model;
}

Model* Model::createReduce(const VectorInt& validVars) const
{
  VectorInt localValidVars = VH::filter(validVars, 0, getNVar());
  Id nvar                  = static_cast<Id>(localValidVars.size());
  if (nvar <= 0)
  {
    messerr("Your new Model has no variable left");
    return nullptr;
  }

  auto* model = new Model(*_ctxt.createReduce(validVars)); // TODO LEAK

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  /* Add the list of Covariances */

  model->setCovAnisoList(getCovAnisoList()->createReduce(validVars));

  return model;
}

VectorECov Model::initCovList(const VectorInt& covranks)
{
  VectorECov list;

  for (Id i = 0; i < static_cast<Id>(covranks.size()); i++)
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
    if (!_driftList->isValid()) return false;
  }

  // Check the consistency between the Covariance and the Drift parts
  Id irf_drift = getDriftMaxIRFOrder();
  Id irf_cova  = getCovMinIRFOrder();
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

const CovAnisoList* Model::getCovAnisoList() const
{
  return castInCovAnisoListConst();
}

CovAnisoList* Model::getCovAnisoListModify()
{
  return _castInCovAnisoList();
}

/**
 * This internal function tries to cast the member '_cova' into a pointer to CovAnisoList
 * and checks the validity of the argument 'icov' which gives the rank within this list
 * @param icov Rank of the CovAniso (to be checked if >= 0)
 * @return 'nullptr' if not valid cast (the error message is printed internally)
 */
const CovAnisoList* Model::castInCovAnisoListConst(Id icov) const
{
  // Check the cast procedure
  const auto* covalist = dynamic_cast<const CovAnisoList*>(_cova.get());
  if (covalist == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovAnisoList");
    return nullptr;
  }
  if (icov < 0) return covalist;

  // Check the rank
  if (icov >= covalist->getNCov())
  {
    messerr("The rank 'icov' (%d) is not valid. The CovAnisoList contains %d covariances",
            icov, covalist->getNCov());
    return nullptr;
  }
  return covalist;
}

CovLMCTapering* Model::_castInCovLMCTapering()
{
  auto* covtape = dynamic_cast<CovLMCTapering*>(_cova.get());
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCTapering");
    return nullptr;
  }
  return covtape;
}

const CovLMCTapering* Model::castInCovLMCTaperingConst() const
{
  const auto* covtape = dynamic_cast<const CovLMCTapering*>(_cova.get());
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCTapering");
    return nullptr;
  }
  return covtape;
}
const CovLMCAnamorphosis* Model::castInCovLMCAnamorphosisConst() const
{
  const auto* covtape = dynamic_cast<const CovLMCAnamorphosis*>(_cova.get());
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCAnamorphosis");
    return nullptr;
  }
  return covtape;
}

CovLMCAnamorphosis* Model::_castInCovLMCAnamorphosis()
{
  auto* covtape = dynamic_cast<CovLMCAnamorphosis*>(_cova.get());
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCAnamorphosis");
    return nullptr;
  }
  return covtape;
}

CovAnisoList* Model::_castInCovAnisoList(Id icov)
{
  // Check the cast procedure
  auto* covalist = dynamic_cast<CovAnisoList*>(_getCovModify());
  if (covalist == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovAnisoList");
    return nullptr;
  }
  if (icov < 0) return covalist;

  // Check the rank
  if (icov >= covalist->getNCov())
  {
    messerr("The rank 'icov' (%d) is not valid. The CovAnisoList contains %d covariances",
            icov, covalist->getNCov());
    return nullptr;
  }
  return covalist;
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
Id Model::stabilize(double percent, bool verbose)
{
  auto nvar = getNVar();
  if (nvar > 1) return 0;
  if (percent <= 0.) return 0;
  auto ncov = getNCov();

  /* Check if the model only contains GAUSSIAN components */

  double total = 0.;
  for (Id icov = 0; icov < ncov; icov++)
  {
    if (getCovAniso(icov)->getType() != ECov::GAUSSIAN) return (0);
    total += getSill(icov, 0, 0);
  }
  total = total * percent / 100.;

  /* Update each Gaussian component */

  for (Id icov = 0; icov < ncov; icov++)
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
Id Model::standardize(bool verbose)

{
  auto nvar = getNVar();
  auto ncov = getNCov();
  VectorDouble total(nvar, 0.);

  /* Calculate the total sills for each variable */

  bool flag_norm = false;
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    total[ivar] = getTotalSill(ivar, ivar);
    if (isZero(total[ivar])) return 1;
    total[ivar] = sqrt(total[ivar]);
    if (ABS(total[ivar] - 1.) > EPSILON6) flag_norm = true;
  }

  /* Scale the different sills for the different variables */

  for (Id ivar = 0; ivar < nvar; ivar++)
    for (Id jvar = 0; jvar < nvar; jvar++)
      for (Id icov = 0; icov < ncov; icov++)
      {
        double sill = getSill(icov, ivar, jvar);
        sill /= total[ivar] * total[jvar];
        setSill(icov, ivar, jvar, sill);
      }

  /* Printout */

  if (verbose && flag_norm)
  {
    message("The model has been normalized\n");
    for (Id ivar = 0; ivar < nvar; ivar++)
      message("- Variable %d : Scaling factor = %lf\n", ivar + 1,
              total[ivar] * total[ivar]);
  }
  return 0;
}

void Model::gofDisplay(double gof,
                       bool byValue,
                       const VectorDouble& thresholds)
{
  ACov::gofDisplay(gof, byValue, thresholds);
}

Id Model::getNVar() const
{
  Id nvar = _cova->getNVar();
  if (nvar <= 0)
    nvar = _ctxt.getNVar();
  return nvar;
}

Model* Model::createFillRandom(Id ndim,
                               Id nvar,
                               const std::vector<ECov>& types,
                               double hmax,
                               Id order,
                               Id nfex,
                               Id seed)
{
  // Create the Covariance Part
  Model* model  = Model::create(CovContext(nvar, ndim));
  Id ncov       = static_cast<Id>(types.size());
  Id seed_local = seed;
  for (Id icov = 0; icov < ncov; icov++)
  {
    MatrixSymmetric* sill =
      MatrixSymmetric::createRandomDefinitePositive(nvar, seed_local);
    double range = (hmax * (1. + icov)) / (2. * ncov);
    model->addCovFromParam(types[icov], range, 0., 1., VectorDouble(), *sill);
    delete sill;
    seed_local = 0;
  }

  // Create the Drift part
  if (order < 0)
  {
    VectorDouble means = VH::simulateGaussian(nvar);
    model->setMeans(means);
  }
  else
  {
    model->setDriftIRF(order, nfex);
  }
  return model;
}

#ifdef HDF5
bool Model::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto modelG = SerializeHDF5::getGroup(grp, "Model");
  if (!modelG) return false;

  bool ret     = true;
  Id ndim      = 0;
  Id nvar      = 0;
  Id ncov      = 0;
  Id ndrift    = 0;
  double field = 0.;

  ret = ret && SerializeHDF5::readValue(*modelG, "NDim", ndim);
  ret = ret && SerializeHDF5::readValue(*modelG, "NVar", nvar);
  ret = ret && SerializeHDF5::readValue(*modelG, "Field", field);
  ret = ret && SerializeHDF5::readValue(*modelG, "NCov", ncov);
  ret = ret && SerializeHDF5::readValue(*modelG, "NDrift", ndrift);
  if (!ret) return ret;

  // Process the Context
  _ctxt = CovContext(nvar, ndim);
  _ctxt.setField(field);
  _clear();
  _create();

  // Process the Covariances
  auto covsG = SerializeHDF5::getGroup(*modelG, "Covs");
  if (!covsG) return false;
  CovAnisoList covs(_ctxt);
  for (Id icov = 0; ret && icov < ncov; icov++)
  {
    String locName = "Covariance" + std::to_string(icov);
    auto covG      = SerializeHDF5::getGroup(*covsG, locName);
    if (!covG) return false;

    // General characteristics
    Id vartype       = 0;
    double range     = 0.;
    double param     = 0.;
    Id flag_aniso    = 0;
    Id flag_rotation = 0;
    VectorDouble aniso_ranges;
    VectorDouble aniso_rotmat;
    VectorDouble sills;

    ret = ret && SerializeHDF5::readValue(*covG, "Type", vartype);
    ret = ret && SerializeHDF5::readValue(*covG, "Range", range);
    ret = ret && SerializeHDF5::readValue(*covG, "Param", param);

    // Anisotropy
    ret = ret && SerializeHDF5::readValue(*covG, "FlagAniso", flag_aniso);
    if (flag_aniso)
    {
      ret = ret && SerializeHDF5::readVec(*covG, "Aniso", aniso_ranges);

      // Rotation
      ret = ret && SerializeHDF5::readValue(*covG, "FlagRotation", flag_rotation);
      if (flag_rotation)
        ret = ret && SerializeHDF5::readVec(*covG, "Rotation", aniso_rotmat);
    }

    // Sills
    ret = ret && SerializeHDF5::readVec(*covG, "Sills", sills);
    if (!ret) return ret;

    CovAniso cova(ECov::fromValue(vartype), _ctxt);
    cova.setParam(param);
    cova.setSill(sills);
    if (flag_aniso)
    {
      for (Id idim = 0; idim < ndim; idim++)
        aniso_ranges[idim] *= range;
      cova.setRanges(aniso_ranges);
      if (flag_rotation) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRangeIsotropic(range);
    covs.addCov(cova);
  }
  setCovAnisoList(&covs);

  // Process the drift part
  auto driftsG = SerializeHDF5::getGroup(*modelG, "Drifts");
  if (!driftsG) return false;
  DriftList drifts(_ctxt);
  ADrift* drift;
  String driftname;
  for (Id ibfl = 0; ret && ibfl < ndrift; ibfl++)
  {
    String locName = "Drift" + std::to_string(ibfl);
    auto driftG    = SerializeHDF5::getGroup(*driftsG, locName);
    if (!driftG) return false;

    ret = ret && SerializeHDF5::readValue(*driftG, "Name", driftname);

    drift = DriftFactory::createDriftByIdentifier(driftname);
    drifts.addDrift(drift);
    delete drift;
  }
  setDriftList(&drifts);

  // Process the Means
  if (ndrift <= 0)
  {
    VectorDouble means;
    ret = ret && SerializeHDF5::readVec(*modelG, "Means", means);
    setMeans(means);
  }

  // Process the covariance matrix
  VectorDouble covar0s;
  ret = ret && SerializeHDF5::readVec(*modelG, "Covar", covar0s);
  setCovar0s(covar0s);

  return ret;
}

bool Model::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto modelG = grp.createGroup("Model");

  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(modelG, "NDim", static_cast<Id>(getNDim()));
  ret      = ret && SerializeHDF5::writeValue(modelG, "NVar", getNVar());
  ret      = ret && SerializeHDF5::writeValue(modelG, "Field", getField());
  ret      = ret && SerializeHDF5::writeValue(modelG, "NCov", getNCov());
  ret      = ret && SerializeHDF5::writeValue(modelG, "NDrift", getNDrift());

  // Writing the covariance part
  auto covsG = modelG.createGroup("Covs");
  for (Id icov = 0, ncov = getNCov(); ret && icov < ncov; icov++)
  {
    const CovAniso* cova = getCovAniso(icov);
    String locName       = "Covariance" + std::to_string(icov);
    auto covG            = covsG.createGroup(locName);

    // General characteristics
    ret = ret && SerializeHDF5::writeValue(covG, "Type", cova->getType().getValue());
    ret = ret && SerializeHDF5::writeValue(covG, "Range", cova->getRangeIso());
    ret = ret && SerializeHDF5::writeValue(covG, "Param", cova->getParam());

    // Anisotropy
    ret = ret && SerializeHDF5::writeValue(covG, "FlagAniso", static_cast<Id>(cova->getFlagAniso()));
    if (cova->getFlagAniso())
    {
      ret = ret && SerializeHDF5::writeVec(covG, "Aniso", cova->getAnisoCoeffs());

      ret = ret && SerializeHDF5::writeValue(covG, "FlagRotation", static_cast<Id>(cova->getFlagRotation()));
      if (cova->getFlagRotation())
        ret = ret && SerializeHDF5::writeVec(covG, "Rotation", cova->getAnisoRotMat().getValues());
    }

    // Sills
    ret = ret && SerializeHDF5::writeVec(covG, "Sills", cova->getSill().getValues());
  }

  // Writing the drift part
  auto driftsG = modelG.createGroup("Drifts");
  for (Id ibfl = 0, nbfl = getNDrift(); ret && ibfl < nbfl; ibfl++)
  {
    const ADrift* drift = getDrift(ibfl);
    String locName      = "Drift" + std::to_string(ibfl);
    auto driftG         = driftsG.createGroup(locName);

    ret = ret && SerializeHDF5::writeValue(driftG, "Name", drift->getDriftName());
  }

  // Writing the matrix of means (if nbfl <= 0)
  if (getNDrift() <= 0)
    ret = ret && SerializeHDF5::writeVec(modelG, "Means", getMeans());

  /// Writing the variance-covariance at the origin (optional)
  ret = ret && SerializeHDF5::writeVec(modelG, "Covar", getCovar0());

  return ret;
}
#endif
} // namespace gstlrn
