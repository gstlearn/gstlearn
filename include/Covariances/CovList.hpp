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
 * This class describes the **Covariance** as a list of elementary covariances (see CovAniso.hpp for more details)
 * where the calculation rule is simple: the returned value is the **sum** of each elementary (active) covariance function.
 *
 * This class also carry two other important informations:
 * - a vector giving the status of each elementary covariance item: it may be *active* or *filtered*
 * - a complex structure allowing each parameter (range, sill, anisotropy angle, ...) of each of the elementary covariances
 * to be non-stationary (to have a value which depends on the location). For more details, see ANoStat.hpp.
 */
class GSTLEARN_EXPORT CovList : public ACov
{
public:
  CovList(const ASpaceSharedPtr& space = ASpaceSharedPtr());
  CovList(const CovList &r) = delete;
  CovList& operator= (const CovList &r) = delete;
  virtual ~CovList();

  /// Interface for ASpaceObject
  virtual bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  virtual int    getNVariables() const override;
  virtual bool   isIndexable() const override { return true; }
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual void addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
  virtual void _addEvalCovMatBiPointInPlace(
                              MatrixSquareGeneral &mat,
                              const SpacePoint &p1,
                              const SpacePoint &p2,
                              const CovCalcMode *mode = nullptr) const override;
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
  // Filter a covariance
  void setFiltered(int icov, bool filtered);

  int             getCovaNumber() const;
  bool            isFiltered(int icov) const;
  virtual double  getTotalSill(int ivar = 0, int jvar = 0) const;
  MatrixSquareSymmetric getTotalSills() const;
  VectorInt       getActiveCovList() const;
  VectorInt       getAllActiveCovList() const;
  bool            isAllActiveCovList() const;
  bool            isNoStat() const override;
  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovBase*    getCova(int icov) const;
  virtual String getCovName(int icov) const;
  virtual const ECov& getType(int icov) const;
  virtual void      setCova(int icov, const CovBase* covs);
  void               setSill(int icov, int ivar, int jvar, double value);
  const MatrixSquareSymmetric& getSill(int icov) const;
  double             getSill(int icov, int ivar, int jvar) const;

  // Methods necessary for Optimization
  void _optimizationPreProcess(const std::vector<SpacePoint> &p) const override;
  void _optimizationPostProcess() const override ;
  void _optimizationSetTarget(const SpacePoint &pt) const override;
  void optimizationSetTargetByIndex(int iech) const override;

protected:
  bool _isCovarianceIndexValid(int icov) const;
  void _loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const override;
  double _loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const;

protected:
 static bool _considerAllCovariances(const CovCalcMode* mode);

private:
  virtual void _delCov(int icov)
  {
    DECLARE_UNUSED(icov)
  };
  // Remove all elementary covariance structures
  virtual void _delAllCov(){};
  void _manage(const Db* db1,const Db* db2) const override;
  

#ifndef SWIG
protected:
  std::vector<const CovBase*> _covs;      /// Vector of elementary covariances
  VectorBool             _filtered; /// Vector of filtered flags (size is nb. cova)
#endif
};
