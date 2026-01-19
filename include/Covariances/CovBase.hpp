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
#pragma once
#include "Basic/AFunctional.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Iterators.hpp"
#include "Basic/ParamInfo.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Covariances/TabNoStatSills.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/CovInternal.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
class AFunctional;
class CovInternal;

namespace gstlrn
{
class GSTLEARN_EXPORT CovBase: public ACov
{
public:
  CovBase(ACov* cor = nullptr, const MatrixSymmetric& sills = MatrixSymmetric());
  CovBase(const CovBase& r);
  CovBase& operator=(const CovBase& r);
  virtual ~CovBase();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovBase)

  /// Interface to ASerializable
  String getNFName() const override { return "CovBase"; }
#ifdef HDF5
  bool deserializeH5(H5::Group& grp) override;
  bool serializeH5(H5::Group& grp) const override;
#endif

  static ParamInfo createParamInfoForCholSill();

  bool isConsistent(const ASpace* space) const override;
  Id getNVar() const override { return _ctxt.getNVar(); }
  bool isOptimizationInitialized(const Db* db = nullptr) const;

  void setCholSill(Id ivar, Id jvar, double val) const;
  virtual void setSill(double sill) const; /// Only valid when there is only one variable (in the context)
  virtual void setSill(const MatrixSymmetric& sill) const;
  virtual void setSill(const VectorDouble& sill) const;
  virtual void setSill(Id ivar, Id jvar, double sill) const;
  void initSill(double value = 0.);
  void setAic(const MatrixSquare& aic) { _aic = aic; }
  void setAic(Id ivar, Id jvar, double val);
  void computeAic();
  void initializeAic();

  const MatrixSymmetric& getSill() const { return _sillCur; }
  const MatrixSquare& getAics() const { return _aic; }
  double getAic(Id ivar, Id jvar) const;
  virtual void setCor(ACov* cor);
  const ACov* getCor() const { return _cor.get(); }
  ACov* getCorModify() { return _cor.get(); }

  double getSill(Id ivar, Id jvar) const;

  void makeSillNoStatDb(const String& namecol, Id ivar = 0, Id jvar = 0, const Db* db = nullptr);
  void makeSillStationary(Id ivar = 0, Id jvar = 0);
  void makeSillsStationary(bool silent = false);
  void makeSillNoStatFunctional(const AFunctional* func, Id ivar = 0, Id jvar = 0);

  TabNoStatSills& getTabNoStatSills() const { return static_cast<TabNoStatSills&>(*_tabNoStat); }

  Id getNSills() const;
  bool isNoStatForVariance() const;

  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  VectorDouble informCoords(const VectorVectorDouble& coords,
                            const EConsElem& econs,
                            Id iv1 = 0,
                            Id iv2 = 0) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;
  void informMeshByMeshForSills(const AMesh* amesh) const;
  void informMeshByApexForSills(const AMesh* amesh) const;
  void informDbInForSills(const Db* dbin) const;
  void informDbOutForSills(const Db* dbout) const;

  /// Tell if the use of Optimization is enabled or not

  void updateCovByPoints(Id icas1, Id iech1, Id icas2, Id iech2) const override;
  void updateCovByMesh(Id imesh, bool aniso = true) const override;

  double getValue(const EConsElem& econs, Id iv1, Id iv2) const override;
  void nostatUpdate(CovInternal* covint) const;
#ifndef SWIG
  Id addEvalCovVecRHSInPlace(vect vect,
                             const VectorInt& index1,
                             Id iech2,
                             const KrigOpt& krigopt,
                             SpacePoint& pin,
                             SpacePoint& pout,
                             VectorDouble& tabwork,
                             double lambda                 = 1.,
                             const ECalcMember& calcMember = ECalcMember::RHS) const override;
#endif
  void setOptimEnabled(bool flag) const override
  {
    _optimEnabled = flag;
    _cor->setOptimEnabled(flag);
  }
  Id makeElemNoStat(const EConsElem& econs, Id iv1, Id iv2, const AFunctional* func = nullptr, const Db* db = nullptr, const String& namecol = String()) override;

  void appendParams(ListParams& listParams,
                    std::vector<covmaptype>* gradFuncs = nullptr) override;
  void updateCov() override;
  void initParams(const MatrixSymmetric& vars,
                  double href = 1.) override;
  ParamInfo& getParamInfoCholSills(Id ivar, Id jvar) { return _cholSillsInfo(ivar, jvar); }

protected:
  void _attachNoStatDb(const Db* db) override;

  void _manage(const Db* db1, const Db* db2) const override;

  bool _checkSill(Id ivar = 0, Id jvar = 0) const;
  bool _checkDims(Id idim, Id jdim) const;

  bool _isVariableValid(Id ivar) const;

  /// Update internal parameters consistency with the context
  void _updateFromContext() override;
  void _initFromContext() override;
  void _copyCovContext(const CovContext& ctxt) override;

private:
  void _makeStationary() override;
  std::unique_ptr<TabNoStat> _createNoStatTab() override;

  bool _isNoStat() const override;
  void _setContext(const CovContext& ctxt) override;

  void _optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const override;
  SpacePoint& _optimizationLoadInPlace(Id iech, Id mode, Id rank) const override;
  void _optimizationPostProcess() const override;

  void _evalOptim(SpacePoint* p1A,
                  SpacePoint* p2A,
                  MatrixSquare& mat,
                  const CovCalcMode* mode) const;
  void _load(const SpacePoint& p, bool case1) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _multiplyCorDerivativesBySills(Id oldSize, std::vector<covmaptype>* gradFuncs);

protected:
  MatrixT<ParamInfo> _cholSillsInfo;
  mutable MatrixSquare _cholSills;
  mutable MatrixSymmetric _sillCur;
  mutable MatrixSquare _workMat;
  mutable MatrixSymmetric _aic;

private:
  std::shared_ptr<ACov> _cor;
  LowerTriangularRange _itRange;
};
} // namespace gstlrn
