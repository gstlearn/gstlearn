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
#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Enum/ECov.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

#include <vector>

class ASpace;
class SpacePoint;
class MatrixSquareGeneral;
class CovBase;
class CovContext;
class AStringFormat;
class AAnam;

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

  /// Interface for ASpaceObject
  virtual bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  virtual int getNVar() const override;
  virtual bool isIndexable() const override { return true; }
  virtual double eval0(int ivar                = 0,
                       int jvar                = 0,
                       const CovCalcMode* mode = nullptr) const override;

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2) const override;

  /// Interface for AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// CovList Interface
  virtual void addCov(const CovBase* cov);

  void addCovList(const CovList* covs);
  // Remove an elementary covariance structure
  void delCov(int icov);
  // Remove all elementary covariance structures
  void delAllCov();
#ifndef SWIG
  int addEvalCovVecRHSInPlace(vect vect,
                              const VectorInt& index1,
                              int iech2,
                              const KrigOpt& krigopt,
                              SpacePoint& pin,
                              SpacePoint& pout,
                              VectorDouble& tabwork,
                              double lambda = 1) const override;
#endif
  void setCovFiltered(int icov, bool filtered);
  int getNCov() const;
  bool isFiltered(int icov) const;
  virtual double getTotalSill(int ivar = 0, int jvar = 0) const;
  MatrixSquareSymmetric getTotalSills() const;
  bool isAllActiveCovList() const;

  void setOptimEnabled(bool flag) const override;
  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovBase* getCov(int icov) const;
  CovBase* getCovModify(int icov);
  virtual String getCovName(int icov) const;
  virtual const ECov& getCovType(int icov) const;
  virtual void setCov(int icov, const CovBase* covs);
  void setSill(int icov, int ivar, int jvar, double value);
  void setSills(int icov, const MatrixSquareSymmetric& sills);
  const MatrixSquareSymmetric& getSills(int icov) const;
  double getSill(int icov, int ivar, int jvar) const;

  // Methods necessary for Optimization
  void _optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationPostProcess() const override;
  SpacePoint& _optimizationLoadInPlace(int iech, int mode, int rank) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

  void setActiveCovListFromOne(int keepOnlyCovIdx) const;
  void setActiveCovListFromInterval(int inddeb, int indto) const;
  void setActiveCovList(const VectorInt& activeCovList, bool allActiveCov) const;

  void copyCovContext(const CovContext& ctxt) override;
  void normalize(double sill = 1., int ivar = 0, int jvar = 0);

  int makeElemNoStat(const EConsElem& econs,
                     int iv1,
                     int iv2,
                     const AFunctional* func = nullptr,
                     const Db* db            = nullptr,
                     const String& namecol   = String()) override;
  void makeSillNoStatDb(int icov, const String& namecol, int ivar = 0, int jvar = 0);
  void makeSillStationary(int icov, int ivar = 0, int jvar = 0);
  void makeSillsStationary(int icov,bool silent = false);
  void makeSillNoStatFunctional(int icov, const AFunctional* func, int ivar = 0, int jvar = 0);

protected:
  bool _isCovarianceIndexValid(int icov) const;
  void _load(const SpacePoint& p, bool case1) const override;

protected:
  const VectorInt& _getListActiveCovariances(const CovCalcMode* mode) const;
  void _updateLists();

  virtual double _eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar                = 0,
                       int jvar                = 0,
                       const CovCalcMode* mode = nullptr) const override;

private:
  void _attachNoStatDb(const Db* db) override;

  void _makeStationary() override;

  bool _isNoStat() const override;
  void _setContext(const CovContext& ctxt) override;
  virtual void _delCov(int icov) { DECLARE_UNUSED(icov) };
  // Remove all elementary covariance structures
  virtual void _delAllCov() {};
  void _manage(const Db* db1, const Db* db2) const override;

#ifndef SWIG

protected:
  std::vector<CovBase*> _covs;         /// Vector of elementary covariances
  VectorBool _filtered;                /// Vector of filtered flags (size is nb. cova)
  mutable bool _allActiveCov;          /*! True if all covariances are active */
  mutable VectorInt _allActiveCovList; /*! List of indices of all covariances */
  mutable VectorInt _activeCovList;    /*! List of indices of the active covariances */
#endif
};
