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
#include "Basic/AStringable.hpp"
#include "Model/ModelCovList.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"
#include "geoslib_f.h"

#include "Enum/ECov.hpp"
#include "Enum/EModelProperty.hpp"

#include "Anamorphosis/AnamHermite.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/CovInternal.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovGradientFunctional.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/ADrift.hpp"
#include "LinearOp/CholeskyDense.hpp"

#include "Db/Db.hpp"

#include <math.h>

Model::Model(const CovContext &ctxt)
    : AStringable(),
      ASerializable(),
      ModelCovList(ctxt)
{
  _create();
}

Model::Model(int nvar, int ndim)
    : AStringable(),
      ASerializable(),
      ModelCovList()
{
  auto space = SpaceRN::create(ndim);
  _ctxt = CovContext(nvar, space);
  _create();
}

Model::Model(const Model &m)
    : AStringable(m),
      ASerializable(m),
      ModelCovList(m._ctxt)
{
  CovAnisoList* mcovalist = dynamic_cast<CovAnisoList*>(m._covList);
  if (mcovalist != nullptr)
    ModelCovList::setCovList(dynamic_cast<CovAnisoList*>(mcovalist->clone()));
  if (m._driftList != nullptr)
    _driftList = m._driftList->clone();
}

Model& Model::operator=(const Model &m)
{
  if (this != &m)
  { 
    AStringable::operator=(m);
    ASerializable::operator=(m);
   setCovAnisoList(dynamic_cast<CovAnisoList*>(m._covList));
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
  int nvar = db->getNLoc(ELoc::Z);
  if (nvar <= 0) nvar = 1;
  auto space = SpaceRN::create(ndim);
  _ctxt = CovContext(nvar, space);
  _create();
  return 0;
}


bool Model::_isValid() const
{
  // Covariances: there should be some defined
  if (_covList == nullptr)
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
  int irf_drift = getDriftMaxIRFOrder();
  int irf_cova  = getCovaMinIRFOrder();
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
                              const MatrixSquareSymmetric& sills,
                              const VectorDouble& angles,
                              const ASpaceSharedPtr& space,
                              bool flagRange)
{
  int nvar = 1;
  if (!sills.empty()) nvar = sills.getNRows();

  auto spaceloc = ASpace::getDefaultSpaceIfNull(space);

  if (!ranges.empty())
  {
    int ndim       = spaceloc->getNDim();
    int ndimRanges = (int)ranges.size();
    if (ndimRanges != 1 && ndimRanges != ndim)
    {
      messerr("Incompatibility between:");
      messerr("Space Dimension = %d", ndim);
      messerr("Dimension of argument 'ranges' = %d", ndimRanges);
      return nullptr;
    }
  }

  CovContext ctxt = CovContext(nvar, space);
  Model* model    = new Model(ctxt);
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
  int nvar = 1;
  if (! sills.empty())
    nvar = (int)  sqrt(sills.size());

  auto spaceloc = ASpace::getDefaultSpaceIfNull(space);
 
  if (! ranges.empty())
  {
    int ndim = spaceloc->getNDim();
    int ndimRanges = (int)ranges.size();
    if (ndimRanges != 1 && ndimRanges != ndim)
    {
      messerr("Incompatibility between:");
      messerr("Space Dimension = %d", ndim);
      messerr("Dimension of argument 'ranges' = %d", ndimRanges);
      return nullptr;
    }
  }

  CovContext ctxt = CovContext(nvar,spaceloc);
  Model* model = new Model(ctxt);
  model->addCovFromParamOldStyle(type, range, sill, param, ranges, sills,
                                 angles, flagRange);

  return model;
}

Model* Model::createFromDb(const Db* db)
{
  Model* model = new Model();
  if (model->resetFromDb(db) != 0)
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

Model* Model::createFromVario(Vario* vario,
                              const VectorECov& types,
                              const Constraints& constraints,
                              const Option_VarioFit& optvar,
                              const Option_AutoFit& mauto,
                              bool verbose)
{
  Model* model = new Model();
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
  int ncov   = getNCov();
  int ndrift = getNDrift();
  if (ncov <= 0 && ndrift <= 0) return sstr.str();

  sstr << toTitle(0, "Model characteristics");
  if (isFlagGradient()) sstr << "(Specific for Handling Gradient)" << std::endl;
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
    sstr << _covList->toString();
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
  
  delete _covList;
  ModelCovList::setCovList(covalist->clone());
}

void Model::addCov(const CovAniso *cov)
{
  if (cov == nullptr) 
  {
    messerr("Error: Covariance is nullptr");
    return;
  }
    
  if (! cov->getContext().isEqual(_ctxt))
  {
    messerr("Error: Covariance should share the same Context as 'Model'");
    messerr("Operation is cancelled");
    return;
  }
  if (_covList == nullptr)
  {
    messerr("Error: Covariance List is nullptr");
    return;
  }
  CovAnisoList* covalist = _castInCovAnisoList();
  if (covalist == nullptr) return;
  covalist->addCovAniso(cov);
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

  int ndim = getNDim();
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
  int nvar = getNVar();
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

  auto space = SpaceRN::create(ndim);
  _ctxt = CovContext(nvar, space);
  CovAniso cov(type, _ctxt);

  // Define the Third parameter
  double parmax = cov.getParMax();
  if (param > parmax) param = parmax;
  cov.setParam(param);

  // Define the range
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

  // Define the sill
  if (!sills.empty())
    cov.setSill(sills);
  else
  {
    if (nvar <= 1)
      cov.setSill(sill);
    else
    {
      MatrixSquareSymmetric locsills(nvar);
      locsills.setIdentity(sill);
      cov.setSill(locsills);
    }
  }

  if (! angles.empty())
    cov.setAnisoAngles(angles);
  addCov(&cov);
}

void Model::addCovFromParam(const ECov& type,
                            double range,
                            double sill,
                            double param,
                            const VectorDouble& ranges,
                            const MatrixSquareSymmetric& sills,
                            const VectorDouble& angles,
                            bool flagRange)
{
  // Check consistency with parameters of the model

  int ndim = getNDim();
  if (!ranges.empty())
  {
    if (ndim > 0 && (int)ranges.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)",
              (int)ranges.size());
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = (int)ranges.size();
  }
  if (!angles.empty())
  {
    if (ndim > 0 && (int)angles.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)",
              (int)angles.size());
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return;
    }
    ndim = (int)angles.size();
  }
  int nvar = getNVar();
  if (!sills.empty())
  {
    if (nvar > 0 && nvar != sills.getNCols())
    {
      messerr("Mismatch between the number of rows 'sills' (%d)", sills.getNRows());
      messerr("and the Number of variables stored in the Model (%d)", nvar);
      messerr("Operation is cancelled");
      return;
    }
    nvar = (int)sqrt((double)sills.size());
  }

  // Define the covariance

  auto space = SpaceRN::create(ndim);
  _ctxt         = CovContext(nvar, space);
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
      MatrixSquareSymmetric locsills(nvar);
      locsills.setIdentity(sill);
      cov.setSill(locsills);
    }
  }

 
  _ctxt.setNVar(cov.getNVar());
  _copyCovContext();
  if (!angles.empty()) cov.setAnisoAngles(angles);
  addCov(&cov);
}



double Model::evalCov(const VectorDouble &incr,
                      int icov,
                      const ECalcMember &member) const
{
  if (_cova == nullptr) return TEST;
  const CovAnisoList* covalist = castInCovAnisoListConst(icov);
  if (covalist == nullptr) return TEST;

  if (member != ECalcMember::LHS && covalist->isFiltered(icov))
    return (0.);
  return getCova(icov)->evalIvarIpas(1., incr);
}


/**
 * Switch to a Model dedicated to Gradients
 * (transforms it from CovAnisoList to CovLMGradient)
 */
void Model::switchToGradient()
{
  // If the Model is already dedicated to Gradient: do nothing
  if (isFlagGradient()) return;

  // If no covariance has been defined yet: do nothing
  if (_cova == nullptr)
  {
    ModelCovList::setCovList(new CovLMGradient(_ctxt));
  }
  else
  {
    const CovAnisoList* covalist = castInCovAnisoListConst();
    if (covalist == nullptr) return;
    ModelCovList::setCovList(new CovLMGradient(*covalist));
  }
}

/**
 * Defining an Anamorphosis information for the Model
 * (in fact, this is added to CovAnisoList part and transforms it from CovAnisoList to CovLMCAnamorphosis
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
    // CovAnisoList is already a covLMCAnamorphosis, simply update the anamorphosis
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
    CovAnisoList* cov = dynamic_cast<CovAnisoList*>(_covList);
    if (cov == nullptr)
    {
      messerr("Impossible to add 'anam' to the covariance part of the Model");
      messerr("The original covariance is probably not a 'CovAnisoList'");
      messerr("The original covariance is probably not a 'CovAnisoList'");
      return 1;
    }

    // Initiate a new CovLMCAnamorphosis class
    CovLMCAnamorphosis* newcov = new CovLMCAnamorphosis(*cov, anam, strcnt);

    // Replace the current list by the newly create one (CovLMCAnamorphosis)

    ModelCovList::setCovList(newcov);

  }
  return 0;
}

int Model::unsetAnam()
{
  if (!hasAnam())
  {
    // CovAnisoList does not have any Anam: do nothing
    return 0;
  }
    CovAnisoList* cov = dynamic_cast<CovAnisoList*>(_cova);
    if (cov == nullptr)
    {
      messerr("Impossible to unset 'anam' from the covariance part of the Model");
      messerr("The original covariance is probably not valid");
      return 1;
    }

  // Initiate a new CovAnisoList class
  CovAnisoList* newcov = new CovAnisoList(*cov);

  // Replace the current list by the newly create one (CovLMCAnamorphosis)

  ModelCovList::setCovList(newcov);

  return 0;
}

void Model::_copyCovContext()
{
  if (_cova == nullptr) return;
  CovAnisoList *covalist = _castInCovAnisoList();
  if (covalist != nullptr) covalist->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

void Model::setField(double field)
{
  _ctxt.setField(field);
  _copyCovContext();
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

  delAllCov();

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

  delAllCov();

  // Add the relevant covariances

  _ctxt = CovContext(vario); /// TODO : What to do with that ?
  _driftList->copyCovContext(_ctxt);

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

  delAllCov();

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

  /* Reading the covariance part and store it into a CovAnisoList */

  CovAnisoList covs(_ctxt);
  for (int icova = 0; ret && icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    ret = ret && _recordRead<int>(is, "Covariance Type", type);
    ret = ret && _recordRead<double>(is, "Isotropic Range", range);
    ret = ret && _recordRead<double>(is, "Model third Parameter", param);
    ret = ret && _recordRead(is, "Flag for Anisotropy", flag_aniso);
    if (! ret) return ret;
    if (flag_aniso != 0)
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
      if (flag_rotation != 0)
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
    if (flag_aniso != 0)
    {
      cova.setRanges(aniso_ranges);
      if (flag_rotation != 0) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRangeIsotropic(range);
    covs.addCovAniso(&cova);
  }
  setCovAnisoList(&covs);

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

  ret = ret && _recordWrite<int>(os, "", getNDim());
  ret = ret && _recordWrite<int>(os, "", getNVar());
  ret = ret && _recordWrite<double>(os, "General parameters", getField());
  ret = ret && _recordWrite<int>(os, "Number of basic covariance terms", getNCov());
  ret = ret && _recordWrite<int>(os, "Number of drift terms", getNDrift());

  /* Writing the covariance part */

  for (int icova = 0; ret && icova < getNCov(); icova++)
  {
    const CovAniso *cova = getCova(icova);
    ret = ret && _recordWrite<int>(os, "", cova->getType().getValue());
    ret = ret && _recordWrite<double>(os, "", cova->getRange());
    ret = ret && _recordWrite<double>(os, "Covariance characteristics", cova->getParam());

    // Writing the Anisotropy information

    ret = ret && _recordWrite<int>(os, "Anisotropy Flag", (int) cova->getFlagAniso());

    if (!cova->getFlagAniso()) continue;

    for (unsigned int idim = 0; ret && idim < getNDim(); idim++)
      ret = ret && _recordWrite<double>(os, "", cova->getAnisoCoeffs(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<int>(os, "Anisotropy Rotation Flag", (int) cova->getFlagRotation());

    if (!cova->getFlagRotation()) continue;

    // Storing the rotation matrix by Column (compatibility)
    for (unsigned int idim = 0; ret && idim < getNDim(); idim++)
      for (unsigned int jdim = 0; ret && jdim < getNDim(); jdim++)
        ret = ret && _recordWrite<double>(os, "", cova->getAnisoRotMat(jdim, idim));
    ret = ret && _commentWrite(os, "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (int ibfl = 0; ret && ibfl < getNDrift(); ibfl++)
  {
    const ADrift *drift = getDrift(ibfl);
    ret = ret && _recordWrite<String>(os,"Drift Identifier", drift->getDriftName());
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (getNDrift() <= 0)
    for (int ivar = 0; ret && ivar < getNVar(); ivar++)
    {
      ret = ret && _recordWrite<double>(os, "Mean of Variables", getMean(ivar));
    }

  /* Writing the matrices of sills (optional) */

  for (int icova = 0; ret && icova < getNCov(); icova++)
  {
    for (int ivar = 0; ret && ivar < getNVar(); ivar++)
      for (int jvar = 0; ret && jvar < getNVar(); jvar++)
        ret = ret && _recordWrite<double>(os, "", getSill(icova, ivar, jvar));
    ret = ret && _commentWrite(os, "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (int ivar = 0; ret && ivar < getNVar(); ivar++)
    for (int jvar = 0; ret && jvar < getNVar(); jvar++)
      ret = ret && _recordWrite<double>(os, "", getContext()->getCovar0(ivar, jvar));
  ret = ret && _commentWrite(os, "Var-Covar at origin");

  return ret;
}

void Model::_clear()
{
  _cova = nullptr;
  _covList = nullptr;
  delete _driftList;
  _driftList = nullptr;
}

void Model::_create()
{
  // TODO: The next two lines are there in order to allow direct call to
  // model::addCov() and model::addDrift
  // The defaulted types of CovAnisoList and DriftList are assumed

  setCovAnisoList(new CovAnisoList(_ctxt));
  _driftList = new DriftList(_ctxt);
}


Model* Model::duplicate() const
{
  Model* model = new Model(*getContext());

  /* Add the list of Covariances */

  model->setCovAnisoList(getCovAnisoList());

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  return model;
}

Model* Model::createReduce(const VectorInt& validVars) const
{
  VectorInt localValidVars = VH::filter(validVars, 0, getNVar());
  int nvar = (int) localValidVars.size();
  if (nvar <= 0)
  {
    messerr("Your new Model has no variable left");
    return nullptr;
  }

  Model* model = new Model(*_ctxt.createReduce(validVars));

  /* Add the list of Covariances */

  model->setCovAnisoList(getCovAnisoList()->createReduce(validVars));

  /* Add the list of Drifts */

  model->setDriftList(getDriftList());

  return model;
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
  const CovAnisoList* covalist = castInCovAnisoListConst(0);
  if (covalist == nullptr) return false;
  const CovGradientNumerical* cova = dynamic_cast<const CovGradientNumerical*>(covalist->getCova(0));
  return (cova != nullptr);
}

bool Model::isFlagGradientFunctional() const
{
  if (! isFlagGradient()) return false;

  // Check is performed on the first covariance
  const CovAnisoList* covalist = castInCovAnisoListConst(0);
  if (covalist == nullptr) return false;
  const CovGradientFunctional* cova = dynamic_cast<const CovGradientFunctional*>(covalist->getCova(0));
  return (cova != nullptr);
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
const CovAnisoList* Model::castInCovAnisoListConst(int icov) const
{
  // Check the cast procedure
  const CovAnisoList* covalist = dynamic_cast<const CovAnisoList*>(_cova);
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
  CovLMCTapering* covtape = dynamic_cast<CovLMCTapering*>(_cova);
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCTapering");
    return nullptr;
  }
  return covtape;
}

CovLMGradient* Model::_castInCovLMGradient()
{
  CovLMGradient* covg = dynamic_cast<CovLMGradient*>(_cova);
  if (covg == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMGradient");
    return nullptr;
  }
  return covg;
}

const CovLMGradient* Model::castInCovLMGradientConst() const
{
  const CovLMGradient* covg = dynamic_cast<const CovLMGradient*>(_cova);
  if (covg == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMGradient");
    return nullptr;
  }
  return covg;
}

const CovLMCTapering* Model::castInCovLMCTaperingConst() const
{
  const CovLMCTapering* covtape = dynamic_cast<const CovLMCTapering*>(_cova);
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCTapering");
    return nullptr;
  }
  return covtape;
}
const CovLMCAnamorphosis* Model::castInCovLMCAnamorphosisConst() const
{
  const CovLMCAnamorphosis* covtape = dynamic_cast<const CovLMCAnamorphosis*>(_cova);
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCAnamorphosis");
    return nullptr;
  }
  return covtape;
}

CovLMCAnamorphosis* Model::_castInCovLMCAnamorphosis() 
{
  CovLMCAnamorphosis* covtape = dynamic_cast<CovLMCAnamorphosis*>(_cova);
  if (covtape == nullptr)
  {
    messerr("The member '_cova' in this model cannot be converted into a pointer to CovLMCAnamorphosis");
    return nullptr;
  }
  return covtape;
}

CovAnisoList* Model::_castInCovAnisoList(int icov)
{
  // Check the cast procedure
  CovAnisoList* covalist = dynamic_cast<CovAnisoList*>(_covList);
  if (covalist == nullptr)
  {
    messerr("The member '_covList' in this model cannot be converted into a pointer to CovAnisoList");
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
int Model::stabilize(double percent, bool verbose)
{
  int nvar = getNVar();
  if (nvar > 1) return 0;
  if (percent <= 0.) return 0;
  int ncov = getNCov();

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
  int nvar = getNVar();
  int ncov = getNCov();
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

void Model::gofDisplay(double gof,
                       bool byValue,
                      const VectorDouble& thresholds)
{
  ACov::gofDisplay(gof, byValue, thresholds);
}

int Model::getNVar() const
{
    // TODO/ the strange next line have been commented out.
    // There should be either validated or suppressed
    //if (isFlagGradient())
    //      return 3; // This strange number of variables is linked to the Gradient calculation
    //    else
    // However, note used for Gradient (Functional type) in Potential
    int nvar = _cova->getNVar();
    if (nvar <= 0)
      nvar = _ctxt.getNVar();
    return nvar;
  }