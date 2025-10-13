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

#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ECov.hpp"
#include "Model/AModelFitSills.hpp"
#include "Space/ASpace.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <vector>

namespace gstlrn
{
class ASpace;
class SpacePoint;
class MatrixSquare;

class AStringFormat;
class AModelFitSills;
class AAnam;
class CovBase;
class CovContext;

/**
 * \brief
 * This class describes the **Covariance** as a list of elementary covariances (see CovBase.hpp for more details)
 * where the calculation rule is simple: the returned value is the **sum** of each elementary (active) covariance function.
 *
 * This class also carry two other important informations:
 * - a vector giving the status of each elementary covariance item: it may be *active* or *filtered*
 * - a complex structure allowing each parameter (range, sill, anisotropy angle, ...) of each of the elementary covariances
 * to be non-stationary (to have a value which depends on the location). For more details, see ANoStat.hpp.
 */
class GSTLEARN_EXPORT CovList: public ACov
{
public:
  CovList(const CovContext& ctxt = CovContext());
  CovList(const CovList& r);
  CovList& operator=(const CovList& r);
  virtual ~CovList();

  /// Interface for ACov
  Id getNVar() const override;
  bool isIndexable() const override { return true; }
  double eval0(Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

  void updateCovByPoints(Id icas1, Id iech1, Id icas2, Id iech2) const override;

  /// Interface for AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// CovList Interface
  virtual void addCov(const CovBase& cov);

  void addCovList(const CovList& covs);
  // Remove an elementary covariance structure
  void delCov(Id icov);
  // Remove all elementary covariance structures
  void delAllCov();
#ifndef SWIG
  Id addEvalCovVecRHSInPlace(vect vect,
                             const VectorInt& index1,
                             Id iech2,
                             const KrigOpt& krigopt,
                             SpacePoint& pin,
                             SpacePoint& pout,
                             VectorDouble& tabwork,
                             double lambda                 = 1,
                             const ECalcMember& calcMember = ECalcMember::RHS) const override;
#endif
  void setCovFiltered(Id icov, bool filtered);
  Id getNCov() const;
  Id getNCovNuggetExcluded() const;
  bool isFiltered(Id icov) const;
  virtual double getTotalSill(Id ivar = 0, Id jvar = 0) const;
  MatrixSymmetric getTotalSills() const;
  bool isAllActiveCovList() const;

  void setOptimEnabled(bool flag) const override;
  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovBase* getCov(Id icov) const;
  CovBase* getCovModify(Id icov);
  virtual String getCovName(Id icov) const;
  virtual const ECov& getCovType(Id icov) const;
  virtual void setCov(Id icov, const CovBase* covs);
  void setSill(Id icov, Id ivar, Id jvar, double value);
  void setSills(Id icov, const MatrixSymmetric& sills);
  const MatrixSymmetric& getSills(Id icov) const;
  double getSill(Id icov, Id ivar, Id jvar) const;

  // Methods necessary for Optimization
  void _optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationPostProcess() const override;
  SpacePoint& _optimizationLoadInPlace(Id iech, Id mode, Id rank) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

  void setActiveCovListFromOne(Id keepOnlyCovIdx) const;
  void setActiveCovListFromInterval(Id inddeb, Id indto) const;
  void setActiveCovList(const VectorInt& activeCovList, bool allActiveCov) const;

  void copyCovContext(const CovContext& ctxt) override;
  void normalize(double sill = 1., Id ivar = 0, Id jvar = 0);

  Id makeElemNoStat(const EConsElem& econs,
                    Id iv1,
                    Id iv2,
                    const AFunctional* func = nullptr,
                    const Db* db            = nullptr,
                    const String& namecol   = String()) override;
  void makeSillNoStatDb(Id icov, const String& namecol, Id ivar = 0, Id jvar = 0);
  void makeSillStationary(Id icov, Id ivar = 0, Id jvar = 0);
  void makeSillsStationary(Id icov, bool silent = false);
  void makeSillNoStatFunctional(Id icov, const AFunctional* func, Id ivar = 0, Id jvar = 0);

  void appendParams(ListParams& listParams,
                    std::vector<covmaptype>* gradFuncs = nullptr) override;
  void updateCov() override;
  void initParams(const MatrixSymmetric& vars, double href = 1.) override;
  void deleteFitSills() const;

  void setFitSills(AModelFitSills* amopts) const;
  AModelFitSills* getFitSills() const;
  Id getNitergCum() const { return _itergCum; }
  void setAllCovActive();

  bool isValidForSpectral() const override;
  MatrixDense simulateSpectralOmega(int ns) const override;
protected:
  bool _isCovarianceIndexValid(Id icov) const;
  void _load(const SpacePoint& p, bool case1) const override;

protected:
  const VectorInt& _getListActiveCovariances(const CovCalcMode* mode) const;
  void _updateLists();

  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

private:
  void _attachNoStatDb(const Db* db) override;

  void _makeStationary() override;

  bool _isNoStat() const override;
  void _setContext(const CovContext& ctxt) override;
  virtual void _delCov(Id icov) { DECLARE_UNUSED(icov) };
  // Remove all elementary covariance structures
  virtual void _delAllCov() {};
  void _manage(const Db* db1, const Db* db2) const override;

#ifndef SWIG

protected:
  std::vector<std::shared_ptr<CovBase>> _covs; /// Vector of elementary covariances
  VectorBool _filtered;                        /// Vector of filtered flags (Dimension = nbcova)
  mutable bool _allActiveCov;                  /*! True if all covariances are active */
  mutable VectorInt _allActiveCovList;         /*! List of indices of all covariances */
  mutable VectorInt _activeCovList;            /*! List of indices of the active covariances */
#endif

private:
  mutable AModelFitSills* _modelFitSills; /* Model fitting procedure for Sills */
  mutable Id _itergCum;
};
} // namespace gstlrn
