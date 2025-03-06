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

#include "Anamorphosis/AnamHermite.hpp"
#include "Enum/EModelProperty.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ECov.hpp"

#include "Basic/ICloneable.hpp"
#include "Covariances/CovList.hpp"
#include "Covariances/CovCalcMode.hpp"

class ASpace;
class SpacePoint;
class MatrixSquareGeneral;
class CovAniso;
class CovContext;
class AStringFormat;
class AAnam;
class AnamHermite;
class EModelProperty;

/**
 * \brief
 * This class describes the **Covariance** as a list of elementary covariances 
 * (see CovAniso.hpp for more details)
 * where the calculation rule is simple: the returned value is the **sum** of each elementary (active) covariance function.
 */
class GSTLEARN_EXPORT CovAnisoList : public CovList
// TODO : rename CovAnisoList (this is not an abstract class)
{
public:
  CovAnisoList(const CovContext& ctxt = CovContext());
  CovAnisoList(const CovAnisoList &r);
  CovAnisoList& operator= (const CovAnisoList &r);
  virtual ~CovAnisoList();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovAnisoList)

  /// Interface for ASpaceObject
  virtual bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  virtual int    getNVar() const override;
  virtual bool   isIndexable() const override { return true; }
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;

  /// Interface for AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// CovAnisoList Interface
  void addCov(const CovBase* cov) override;
  const AnamHermite* getAnamHermite() const;

  const EModelProperty& getCovMode() const;
  virtual bool hasAnam() const { return false; }
  virtual const AAnam* getAnam() const { return nullptr; }
  virtual void setActiveFactor(int /*iclass*/) { }
  virtual int getActiveFactor() const { return 0; }
  virtual int getAnamNClass() const { return 0; }

  void addCovList(const CovAnisoList* covs);

  int             getNCov(bool skipNugget = false) const;
  bool            hasRange() const;
  bool            isStationary() const;
  double          getMaximumDistance() const;
  double          getTotalSill(int ivar = 0, int jvar = 0) const override;
  bool            isNoStat() const override;
  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovAniso*    getCovAniso(int icov) const;
  CovAniso*          getCovAniso(int icov); // TODO : beurk :(
  void               setCov(int icov, const CovBase* cov) override;
  const ECov&        getCovType(int icov) const override;
  String             getCovName(int icov) const override;
  void               setRangeIsotropic(int icov, double range);
  void               setType(int icov, const ECov& type);
  void               setParam(int icov, double value);
  void               setMarkovCoeffs(int icov, const VectorDouble& coeffs);
  double             getParam(int icov) const;
  double             getRange(int icov) const;
  VectorDouble       getRanges(int icov) const;
  VectorDouble       getAngles(int icov) const;
  int                getNGradParam(int icov) const;
  CovAniso           extractCova(int icov) const;
  int                getCovMinIRFOrder() const;
  double             getBallRadius() const;
  int                hasExternalCov() const;
  bool               isChangeSupportDefined() const;
  // Methods necessary for Optimization

  bool hasNugget() const;
  int  getRankNugget() const;
  const CovAnisoList* createReduce(const VectorInt& validVars) const;
  void setOptimEnabled(bool status);

private:
  // Returns a pointer on an existing Cov and cast it to CovAniso
  const CovAniso* _getCovAniso(int icov) const;
  CovAniso* _getCovAnisoModify(int icov);

protected:
  bool _isCovarianceIndexValid(int icov) const;
};
