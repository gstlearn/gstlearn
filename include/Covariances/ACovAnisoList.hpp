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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ECov.hpp"

#include "Basic/ICloneable.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Model/ANoStat.hpp"

#include <vector>

class ASpace;
class SpacePoint;
class MatrixSquareSymmetric;
class CovAniso;
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
class GSTLEARN_EXPORT ACovAnisoList : public ACov, public ICloneable
// TODO : rename CovAnisoList (this is not an abstract class)
{
public:
  ACovAnisoList(const ASpace* space = nullptr);
  ACovAnisoList(const ACovAnisoList &r);
  ACovAnisoList& operator= (const ACovAnisoList &r);
  virtual ~ACovAnisoList();

  /// ICloneable interface
  IMPLEMENT_CLONING(ACovAnisoList)

  /// Interface for ASpaceObject
  virtual bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  virtual int    getNVariables() const override;
  virtual bool   isIndexable() const override { return true; }
  virtual bool   isNoStat() const override { return _noStat != nullptr; }
  virtual const ANoStat* getNoStat() const override { return _noStat; }
  virtual ANoStat* getNoStatModify() const override { return _noStat; } // TODO: to be suppressed
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual void eval0MatInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
  virtual void evalMatInPlace(const SpacePoint &p1,
                              const SpacePoint &p2,
                              MatrixSquareGeneral &mat,
                              const CovCalcMode *mode = nullptr) const override;
  virtual void evalMatOptimInPlace(int icas1,
                                   int iech1,
                                   int icas2,
                                   int iech2,
                                   MatrixSquareGeneral &mat,
                                   const CovCalcMode *mode = nullptr) const override;
  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2) override;
  virtual void updateCovByMesh(int imesh) override;

  /// Interface for AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ACovAnisoList Interface
  virtual void addCov(const CovAniso* cov);
  virtual bool hasAnam() const { return false; }
  virtual const AAnam* getAnam() const { return nullptr; }
  virtual void setActiveFactor(int /*iclass*/) { return; }
  virtual int getActiveFactor() const { return 0; }
  virtual int getAnamNClass() const { return 0; }

  void addCovList(const ACovAnisoList* covs);
  // Remove an elementary covariance structure
  void delCov(unsigned int i);
  // Remove all elementary covariance structures
  void delAllCov();
  // Filter a covariance
  void setFiltered(unsigned int i, bool filtered);

  int             getCovaNumber(bool skipNugget = false) const;
  bool            isFiltered(unsigned int i) const;
  bool            hasRange() const;
  bool            isStationary() const;
  double          getMaximumDistance() const;
  double          getTotalSill(int ivar, int jvar) const;
  MatrixSquareGeneral getTotalSill() const;
  void            normalize(double sill = 1., int ivar=0, int jvar=0);
  VectorInt       getActiveCovList() const;
  VectorInt       getAllActiveCovList() const;
  bool            isAllActiveCovList() const;

  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovAniso*    getCova(int icov) const;
  CovAniso*          getCova(int icov); // TODO : beurk :(
  void               setCova(int icov, CovAniso* covs);
  const ECov&        getType(int icov) const;
  String             getCovName(int icov) const;
  double             getParam(unsigned int icov) const;
  double             getRange(int icov) const;
  VectorDouble       getRanges(int icov) const;
  const MatrixSquareSymmetric& getSill(unsigned int icov) const;
  double             getSill(unsigned int icov, int ivar, int jvar) const;
  int                getGradParamNumber(unsigned int icov) const;
  void               setSill(unsigned int icov, int ivar, int jvar, double value);
  void               setType(unsigned int icov, const ECov& type);
  void               setParam(unsigned int icov, double value);
  CovAniso           extractCova(int icov) const;
  int                getCovaMinIRFOrder() const;

  // Methods necessary for Optimization
  void optimizationPreProcess(const std::vector<SpacePoint> &vec) const;
  void optimizationPostProcess() const;
  void optimizationSetTarget(const SpacePoint &pt) const;
  void evalOptimInPlace(VectorDouble &res,
                        int ivar = 0,
                        int jvar = 0,
                        const CovCalcMode *mode = nullptr) const;
  VectorVectorDouble evalCovMatrixOptim(const Db *db1,
                                        const Db *db2,
                                        int ivar,
                                        int jvar,
                                        const CovCalcMode *mode) const;
  ////////////////////////////////////////////////

  void copyCovContext(const CovContext& ctxt);
  bool hasNugget() const;

  const ACovAnisoList* createReduce(const VectorInt &validVars) const;

  int addNoStat(const ANoStat *anostat);
  void delNoStat();
  int getNoStatElemNumber() const;
  const EConsElem& getNoStatElemType(int ipar) const;
  int addNoStatElem(int igrf,
                    int icov,
                    const EConsElem &type,
                    int iv1,
                    int iv2);
  int addNoStatElems(const VectorString &codes);
  CovParamId getCovParamId(int ipar) const;
  int getNoStatElemIcov(int ipar) const;

protected:
  bool   _isCovarianceIndexValid(unsigned int i) const;

private:
  bool _considerAllCovariances(const CovCalcMode* mode) const;

#ifndef SWIG
protected:
  std::vector<CovAniso*> _covs;     /// Vector of elementary covariances
  VectorBool             _filtered; /// Vector of filtered flags (size is nb. cova)
  ANoStat*               _noStat;   /// Description of Non-stationary Model
#endif
};
