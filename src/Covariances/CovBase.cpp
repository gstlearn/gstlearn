/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c3) MINES Paris / ARMINES                                       */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovBase.hpp"
#include "Basic/Iterators.hpp"
#include "Basic/ListParams.hpp"
#include "Basic/ParamInfo.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/TabNoStatSills.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include <cmath>
#include <cstddef>
#include <memory>

namespace gstlrn
{

  ParamInfo CovBase::createParamInfoForCholSill()
  {
    return ParamInfo(
      String("Cholesky sill"),
      TEST,
      {-INF, INF},
      String("Term of the Cholesky decomposition of the sill matrix"));
  }

  CovBase::CovBase(ACov* cor, const MatrixSymmetric& sill)
    : ACov(cor == nullptr ? CovContext() : cor->getContext())
    , _cholSillsInfo(MatrixT<ParamInfo>(sill.getNRows(),
                                        sill.getNCols(),
                                        createParamInfoForCholSill()))
    , _cholSills(MatrixSymmetric(sill.getNRows()))
    , _sillCur(sill)
    , _cor(std::dynamic_pointer_cast<ACov>(cor == nullptr ? nullptr
                                                          : cor->cloneShared()))
    , _itRange(sill.getNRows())
  {

    createNoStatTab();

    _ctxt.setNVar(sill.getNSize());
    _aic.clear();

    if (cor != nullptr)
    {
      _ctxt = cor->getContextCopy();
    }
    _ctxt.setNVar(sill.getNSize());
    _workMat.resize(_ctxt.getNVar(), _ctxt.getNVar());
    _workMat.setIdentity();
  }

  CovBase::CovBase(const CovBase& r)
    : ACov(r)
    , _itRange(r._cholSills.getNRows())
  {
    _cholSillsInfo = r._cholSillsInfo;
    _cholSills = r._cholSills;
    _sillCur = r._sillCur;
    _workMat = r._workMat;
    _aic = r._aic;
    _cor = std::dynamic_pointer_cast<ACov>(r._cor->cloneShared());
  }

  CovBase& CovBase::operator=(const CovBase& r)
  {
    if (this != &r)
    {
      _cholSillsInfo = r._cholSillsInfo;
      _cholSills = r._cholSills;
      _sillCur = r._sillCur;
      _workMat = r._workMat;
      _aic = r._aic;
      _cor = std::dynamic_pointer_cast<ACov>(r._cor->cloneShared());
      _itRange = LowerTriangularRange(r._cholSills.getNRows());
    }
    return *this;
  }

  CovBase::~CovBase() {}

  void CovBase::setCor(ACov* cor)
  {
    _cor = std::dynamic_pointer_cast<ACov>(cor->cloneShared());
    auto nvar = getNVar();
    if (cor != nullptr)
    {
      _ctxt = cor->getContextCopy();
      _ctxt.setNVar(nvar);
    }
  }

  void CovBase::_setContext(const CovContext& ctxt)
  {
    _cor->setContext(ctxt);
    _updateFromContext();
  }

  void CovBase::setSill(double sill) const
  {
    auto nvar = getNVar();
    if (nvar > 0 && nvar != 1)
    {
      messerr("Number of provided sill doesn't match number of variables");
      return;
    }
    _sillCur.resetFromValue(1, 1, sill);
  }

  void CovBase::setSill(const MatrixSymmetric& sill) const
  {
    auto nvar = getNVar();
    if (nvar > 0 && nvar != sill.getNCols())
    {
      messerr("Number of provided sills doesn't match number of variables");
      return;
    }
    _sillCur = sill;
  }

  void CovBase::setSill(const VectorDouble& sill) const
  {
    Id size = static_cast<Id>(sill.size());
    auto nvar = getNVar();
    if (size != nvar * nvar)
    {
      messerr("Number of provided sills doesn't match number of variables");
      return;
    }
    _sillCur.setValues(sill);
  }

  void CovBase::setSill(Id ivar, Id jvar, double sill) const
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
    /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
    if (!_sillCur.isValid(ivar, jvar)) return;
    _sillCur.setValue(ivar, jvar, sill);
  }

  double CovBase::getAic(Id ivar, Id jvar) const
  {
    if (!_isVariableValid(ivar)) return 0.;
    if (!_isVariableValid(jvar)) return 0.;
    return _aic.getValue(ivar, jvar);
  }

  void CovBase::setCholSill(Id ivar, Id jvar, double val) const
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
    if (ivar < jvar)
    {
      messerr(
        "The Cholesky decomposition of the sill matrix is lower triangular");
      return;
    }
    _cholSills.setValue(ivar, jvar, val);
  }

  void CovBase::setAic(Id ivar, Id jvar, double val)
  {
    if (!_isVariableValid(ivar)) return;
    if (!_isVariableValid(jvar)) return;
    _aic.setValue(ivar, jvar, val);
  }

  /**
   * @brief Calculate the square root of the sill matrix.
   * This decomposition requires the Sill matrix to be definite positive.
   */
  void CovBase::computeAic()
  {
    // Check that the sill matrix has been defined before squaring it
    if (!_sillCur.empty()) _aic = _sillCur.squareRoot();
  }

  void CovBase::initializeAic()
  {
    if (_sillCur.empty()) return;
    Id nvar = _sillCur.getNCols();
    _aic = MatrixSymmetric(nvar);
    _aic.setIdentity();
  }

  bool CovBase::_isVariableValid(Id ivar) const
  {
    return checkArg("Rank of the Variable", ivar, getNVar());
  }

  void CovBase::_initFromContext()
  {
    _cor->initFromContext();
    _sillCur.reset(_ctxt.getNVar(), _ctxt.getNVar());
    setOptimEnabled(true);
    _ctxt.setSpace(_cor->getSpace());
  }

  void CovBase::initSill(double value)
  {
    _sillCur.fill(value);
  }

  bool CovBase::isConsistent(const ASpace* space) const
  {
    return _cor->isConsistent(space);
  }

  Id CovBase::addEvalCovVecRHSInPlace(vect vec,
                                      const VectorInt& index1,
                                      Id iech2,
                                      const KrigOpt& krigopt,
                                      SpacePoint& pin,
                                      SpacePoint& pout,
                                      VectorDouble& tabwork,
                                      double lambda,
                                      const ECalcMember& calcMember) const
  {
    DECLARE_UNUSED(lambda)
    return _cor->addEvalCovVecRHSInPlace(vec,
                                         index1,
                                         iech2,
                                         krigopt,
                                         pin,
                                         pout,
                                         tabwork,
                                         getSill(0, 0),
                                         calcMember);
  }

  double CovBase::_eval(const SpacePoint& p1,
                        const SpacePoint& p2,
                        Id ivar,
                        Id jvar,
                        const CovCalcMode* mode) const
  {
    return getSill(ivar, jvar) * _cor->evalCov(p1, p2, ivar, jvar, mode);
  }

  double CovBase::getSill(Id ivar, Id jvar) const
  {
    return _sillCur.getValue(ivar, jvar);
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
  void CovBase::nostatUpdate(CovInternal* covint) const
  {
    if (covint == nullptr) return;
    updateCovByPoints(covint->getIcas1(),
                      covint->getIech1(),
                      covint->getIcas2(),
                      covint->getIech2());
  }

  void CovBase::_copyCovContext(const CovContext& ctxt)
  {
    _ctxt.copyCovContext(ctxt);
    _cor->copyCovContext(ctxt);
  }

  /**
   * Transform a set of Space Points using the anisotropy tensor
   * The set of resulting Space Points are stored as private member of this.
   * Note that ALL samples are processed, independently from the presence of a selection
   * or checking for heterotopy.
   * @param mode 1 for p1As; 2for p2As
   * @param ps vector of SpacePoints
   */
  void CovBase::_optimizationPreProcess(Id mode,
                                        const std::vector<SpacePoint>& ps) const
  {
    _cor->optimizationPreProcess(mode, ps);
  }

  SpacePoint& CovBase::_optimizationLoadInPlace(Id iech, Id mode, Id rank) const
  {
    return _cor->optimizationLoadInPlace(iech, mode, rank);
  }

  /**
   * Checks that the Optimization has already been initiated, by:
   * - checking that the storage (for Sample Points projected in the Covariance
   * rotation system) is already allocated
   * - checking that the dimension of this storage is correct (only if 'db' is provided):
   * in particular, this check is not necessary when freeing this storage.
   */
  bool CovBase::isOptimizationInitialized(const Db* db) const
  {
    if (_p1As.empty()) return false;
    if (db == nullptr) return true;
    Id n = static_cast<Id>(_p1As.size());
    return n == db->getNSample();
  }

  // Set of functions to make parameters no stationary (or to make them back stationary).
  // There is two types of non stationarities : NoStatDb in which the parameters are read in a
  // DbGrid or NoStatFunctional for which you have to provide a function of the coordinates.
  // Each parameter can have its own type of No stationarity and its own DbGrid in case
  // of NoStatDb.
  // For specifying the NoStat DbGrid, you can first attach it by using attachNoStatDb.
  // If not, you have to specify the DbGrid when you make the first parameter non stationary.

  void CovBase::_attachNoStatDb(const Db* db)
  {
    DECLARE_UNUSED(db)
    if (_cor->getDbNoStat() == nullptr)
    {
      std::shared_ptr<const Db> dbptr = _tabNoStat->getDbNoStatRef();
      _cor->setNoStatDbIfNecessary(dbptr);
    }
  }

  Id CovBase::makeElemNoStat(const EConsElem& econs,
                             Id iv1,
                             Id iv2,
                             const AFunctional* func,
                             const Db* db,
                             const String& namecol)
  {
    if (ACov::makeElemNoStat(econs, iv1, iv2, func, db, namecol)) return 1;

    return _cor->makeElemNoStat(econs, iv1, iv2, func, db, namecol);
  }

  ///////////////////// Sill ////////////////////////

  void CovBase::makeSillNoStatDb(const String& namecol,
                                 Id ivar,
                                 Id jvar,
                                 const Db* db)
  {
    if (!_checkSill(ivar, jvar)) return;
    makeElemNoStat(EConsElem::SILL, ivar, jvar, nullptr, db, namecol);
    //_cor->checkAndManageNoStatDb(db, namecol);
  }

  void
    CovBase::makeSillNoStatFunctional(const AFunctional* func, Id ivar, Id jvar)
  {
    if (!_checkSill(ivar, jvar)) return;
    makeElemNoStat(EConsElem::SILL, ivar, jvar, func);
  }

  void CovBase::makeSillsStationary(bool silent)
  {
    if (getTabNoStatSills().empty() && !silent)
    {
      messerr("All the sills are already stationary!");
      return;
    }
    _tabNoStat->clear();
  }

  void CovBase::makeSillStationary(Id ivar, Id jvar)
  {
    if (!_checkSill(ivar, jvar)) return;
    if (_tabNoStat->removeElem(EConsElem::SILL, ivar, jvar) == 0)
    {
      messerr("This parameter was already stationary!");
    }
  }

  /////////////////////////// Check functions ////////////////////:

  bool CovBase::_checkSill(Id ivar, Id jvar) const
  {
    auto nvar = getNVar();
    if ((ivar > nvar) || (jvar > nvar))
    {
      messerr("Your model has only %d variables.", nvar);
      return false;
    }
    return true;
  }

  bool CovBase::_checkDims(Id idim, Id jdim) const
  {
    auto ndim = getNDim();
    if ((idim > ndim) || (jdim > ndim))
    {
      messerr("Your model is only in dimension %d.", ndim);
      return false;
    }
    return true;
  }

  /////////////  Functions to attach no stat information on various supports ////////
  void CovBase::informMeshByMesh(const AMesh* amesh) const
  {
    _tabNoStat->informMeshByMesh(amesh);
    _cor->informMeshByMesh(amesh);
  }

  void CovBase::informMeshByApex(const AMesh* amesh) const
  {
    _tabNoStat->informMeshByApex(amesh);
    _cor->informMeshByApex(amesh);
  }

  void CovBase::informDbIn(const Db* dbin) const
  {
    _tabNoStat->informDbIn(dbin);
    _cor->informDbIn(dbin);
  }

  void CovBase::informDbOut(const Db* dbout) const
  {
    _tabNoStat->informDbOut(dbout);
    _cor->informDbOut(dbout);
  }

  double CovBase::getValue(const EConsElem& econs, Id iv1, Id iv2) const
  {
    double val = _cor->getValue(econs, iv1, iv2);
    if (val == TEST)
    {
      if (econs == EConsElem::SILL) return getSill(iv1, iv2);
    }
    return val;
  }

  VectorDouble CovBase::informCoords(const VectorVectorDouble& coords,
                                     const EConsElem& econs,
                                     Id iv1,
                                     Id iv2) const
  {
    if (econs == EConsElem::SILL)
    {
      VectorDouble result(coords[0].size(), getValue(econs, iv1, iv2));
      _tabNoStat->informCoords(coords, econs, iv1, iv2, result);
      return result;
    }

    return _cor->informCoords(coords, econs, iv1, iv2);
  }

  void CovBase::informMeshByMeshForSills(const AMesh* amesh) const
  {
    _tabNoStat->informMeshByMesh(amesh, EConsElem::SILL);
  }

  void CovBase::informMeshByApexForSills(const AMesh* amesh) const
  {
    _tabNoStat->informMeshByApex(amesh, EConsElem::SILL);
  }

  void CovBase::informDbInForSills(const Db* dbin) const
  {
    _tabNoStat->informDbIn(dbin, EConsElem::SILL);
  }

  void CovBase::informDbOutForSills(const Db* dbout) const
  {
    _tabNoStat->informDbOut(dbout, EConsElem::SILL);
  }

  /**
   * Update the Model according to the Non-stationary parameters
   * @param icas1 Type of first Db: 1 for Input; 2 for Output
   * @param iech1 Rank of the target within Db1 (or -1)
   * @param icas2 Type of first Db: 1 for Input; 2 for Output
   * @param iech2 Rank of the target within Dbout (or -2)
   */
  void CovBase::updateCovByPoints(Id icas1, Id iech1, Id icas2, Id iech2) const
  {
    // If no non-stationary parameter is defined, simply skip
    if (!isNoStat()) return;
    double val1, val2;

    const auto paramsnostat = _tabNoStat->getTable();
    // Loop on the elements that can be updated one-by-one

    for (const auto& e: paramsnostat)
    {
      EConsElem type = e.first.getType();
      e.second->getValuesOnDb(icas1, iech1, &val1, icas2, iech2, &val2);

      if (type == EConsElem::SILL)
      {
        Id iv1 = e.first.getIV1();
        Id iv2 = e.first.getIV2();
        setSill(iv1, iv2, sqrt(val1 * val2));
      }
    }
    _cor->updateCovByPoints(icas1, iech1, icas2, iech2);
  }

  void CovBase::updateCovByMesh(Id imesh, bool aniso) const
  {
    // If no non-stationary parameter is defined, simply skip
    if (!isNoStat()) return;

    // Loop on the elements that can be updated one-by-one
    if (!aniso)
    {
      const auto paramsnostat = _tabNoStat->getTable();
      for (const auto& e: paramsnostat)
      {
        EConsElem type = e.first.getType();
        if (type == EConsElem::SILL)
        {
          double sill = e.second->getValueOnMeshByApex(imesh);
          Id iv1 = e.first.getIV1();
          Id iv2 = e.first.getIV2();
          setSill(iv1, iv2, sill);
        }
      }
    }
    _cor->updateCovByMesh(imesh, aniso);
  }

  std::unique_ptr<TabNoStat> CovBase::_createNoStatTab()
  {
    return std::make_unique<TabNoStatSills>();
  }

  void CovBase::_makeStationary()
  {
    _cor->makeStationary();
  }

  void CovBase::_manage(const Db* db1, const Db* db2) const
  {
    if (db1 != nullptr) informDbIn(db1);
    if (db2 != nullptr) informDbOut(db2);
    _cor->manage(db1, db2);
  }

  void CovBase::_optimizationPostProcess() const
  {
    _cor->optimizationPostProcess();
  }

  void CovBase::_updateFromContext()
  {
    _cor->updateFromContext();
  }

  void CovBase::_load(const SpacePoint& p, bool case1) const
  {
    _cor->load(p, case1);
  }

  void CovBase::_optimizationSetTarget(SpacePoint& pt) const
  {
    _cor->optimizationSetTarget(pt);
  }

  bool CovBase::_isNoStat() const
  {
    return _cor->isNoStat() || isNoStatForVariance();
  }

  Id CovBase::getNSills() const
  {
    TabNoStatSills& tabnostat = getTabNoStatSills();
    return tabnostat.getNSills();
  }

  bool CovBase::isNoStatForVariance() const
  {
    TabNoStatSills& tabnostat = getTabNoStatSills();
    return tabnostat.isDefinedForVariance();
  }

  void
    CovBase::_multiplyCorDerivativesBySills(Id oldSize,
                                            std::vector<covmaptype>* gradFuncs)
  {
    // Multiply the derivatives of the correlation functions by the sills
    // This is done to ensure that the gradient functions reflect the sill scaling

    const auto newSize = gradFuncs->size(); // snapshot après append

    for (size_t i = oldSize; i < newSize; ++i)
    {
      const covmaptype f = (*gradFuncs)[i];
      (*gradFuncs)[i] = [f, this](const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  Id ivar,
                                  Id jvar,
                                  const CovCalcMode* mode)
      {
        double result = f(p1, p2, ivar, jvar, mode) * this->getSill(ivar, jvar);
        return result;
      };
    }
  }

  static double softplus(double x)
  {
    // Softplus function to ensure positive values
    return std::log1p(std::exp(x));
    // stable version of softplus:  log(1 + exp(x))
  }

  static double softplusinv(double x)
  {
    // Inverse of the softplus function
    if (x <= 0) return -std::numeric_limits<double>::infinity();
    return x + std::log1p(-std::exp(-x));
    // stable version of softplus inverse: log(exp(x) - 1)
  }

  static double softplusDerivative(double x)
  {
    // Derivative of the softplus function
    return 1. / (1.0 + exp(-x));
  }

  void CovBase::appendParams(ListParams& listParams,
                             std::vector<covmaptype>* gradFuncs)
  {
    const auto oldSize = gradFuncs->size();
    _cor->appendParams(listParams, gradFuncs);
    _multiplyCorDerivativesBySills(static_cast<Id>(oldSize), gradFuncs);

    for (const auto& [ivar, jvar]: _itRange)
      listParams.addParam(_cholSillsInfo(ivar, jvar));

    for (const auto& pair: _itRange)
    {
      const Id ivard = pair.first;
      const Id jvard = pair.second;
      if (_cholSillsInfo(ivard, jvard).isFixed())
        continue; // Skip fixed parameters
      gradFuncs->emplace_back(
        [ivard, jvard, this](const SpacePoint& p1,
                             const SpacePoint& p2,
                             Id ivar,
                             Id jvar,
                             const CovCalcMode* mode) -> double
        {
          MatrixSymmetric dSillDChol(this->getNVar());
          dSillDChol.fill(0.);
          for (Id i = jvard; i < this->getNVar(); i++)
          {
            double val = this->_cholSillsInfo.getValue(i, jvard).getValue();
            if (i == static_cast<Id>(jvard))
            {
              double grad_softplus = softplusDerivative(val);
              val = grad_softplus
                  * softplus(val); // Apply softplus to ensure positive values
            }
            if (i == static_cast<Id>(ivard))
            {
              val *= 2; // Derivative of the diagonal element is 2
            }
            dSillDChol.setValue(i, ivard, val);
            dSillDChol.setValue(ivard, i, val);
          }
          double cor = this->_eval(p1, p2, ivar, jvar, mode);
          return dSillDChol.getValue(ivar, jvar) * cor;
        });
    }
  }

  void CovBase::initParams(const MatrixSymmetric& vars, double href)
  {
    CholeskyDense chol(vars);
    for (const auto& [ivar, jvar]: _itRange)
    {
      double value = chol.getLowerTriangle(ivar, jvar);
      if (ivar == jvar)
      {
        _cholSillsInfo(ivar, jvar)
          .setMaxValue(softplusinv(
            2.
            * ABS(value))); // Protection against too large values for diagonal
        value = softplusinv(ABS(value));
      }
      _cholSillsInfo(ivar, jvar).setValue(value);
    }

    _cor->initParams(vars, href);
  }

  void CovBase::updateCov()
  {
    _cor->updateCov();
    Id nvaroptim = 0;
    for (const auto& [ivar, jvar]: _itRange)
    {
      if (_cholSillsInfo(ivar, jvar).isFixed()) continue;
      nvaroptim++;
      double val = _cholSillsInfo(ivar, jvar).getValue();
      if (ivar == jvar)
        val = softplus(val); // Apply softplus to ensure positive values

      _cholSills.setValue(ivar, jvar, val);
    }

    if (nvaroptim > 0)
      _sillCur.prodMatMatInPlace(&_cholSills, &_cholSills, false, true);
  }

#ifdef HDF5
  bool CovBase::deserializeH5(H5::Group& grp)
  {
    bool ret = true;

    // Retrieve the Sills
    VectorDouble sills;
    ret = ret && SerializeHDF5::readVec(grp, "Sill Matrix", sills);
    setSill(sills);

    // Retreive the CorAniso
    ret = ret && _cor->deserializeH5(grp);

    // Non stationary case.  Look for "NonStat" paragraph. It not found, simply return ... silently
    auto nostatG = SerializeHDF5::getGroup(grp, "NoStatSills", false);
    if (nostatG) ret = ret && _tabNoStat->deserializeH5(*nostatG);

    return ret;
  }

  bool CovBase::serializeH5(H5::Group& grp) const
  {
    bool ret = true;

    // Serialize the sills
    ret =
      ret && SerializeHDF5::writeVec(grp, "Sill Matrix", getSill().getValues());

    // Serialize the CorAniso
    ret = ret && _cor->serializeH5(grp);

    // Non stationary case
    if (isNoStat() && _tabNoStat->size() > 0)
    {
      auto nonstatG = grp.createGroup("NoStatSills");
      ret = ret && _tabNoStat->serializeH5(nonstatG);
    }

    return ret;
  }
#endif

} // namespace gstlrn
