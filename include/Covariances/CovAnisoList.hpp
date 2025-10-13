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

#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Enum/EModelProperty.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

#include "Basic/ICloneable.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovList.hpp"

namespace gstlrn
{
class AAnam; // Forward declaration
class CovAniso;
class CovContext;
class AnamHermite;
class ASpace;
class SpacePoint;
class MatrixSquare;

class AStringFormat;
class EModelProperty;

/**
 * \brief
 * This class describes the **Covariance** as a list of elementary covariances
 * (see CovAniso.hpp for more details)
 * where the calculation rule is simple: the returned value is the **sum** of each elementary (active) covariance function.
 */
class GSTLEARN_EXPORT CovAnisoList: public CovList
// TODO : rename CovAnisoList (this is not an abstract class)
{
public:
  CovAnisoList(const CovContext& ctxt = CovContext());
  CovAnisoList(const CovAnisoList& r);
  CovAnisoList& operator=(const CovAnisoList& r);
  virtual ~CovAnisoList();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovAnisoList)

  /// Interface for ASpaceObject
  bool isConsistent(const ASpace* space) const override;

  /// Interface for ACov
  Id getNVar() const override;
  bool isIndexable() const override { return true; }
  double eval0(Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

  /// Interface for AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// CovAnisoList Interface
  void addCov(const CovBase& cov) override;
  const AnamHermite* getAnamHermite() const;

  const EModelProperty& getCovMode() const;
  virtual bool hasAnam() const { return false; }
  virtual const AAnam* getAnam() const { return nullptr; }
  virtual void setActiveFactor(Id /*iclass*/) {}
  virtual Id getActiveFactor() const { return 0; }
  virtual Id getAnamNClass() const { return 0; }

  void addCovList(const CovAnisoList& covs);

  Id getNCov(bool skipNugget = false) const;
  bool hasRange() const;
  bool isStationary() const;
  double getMaximumDistance() const;
  double getTotalSill(Id ivar = 0, Id jvar = 0) const override;
  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const CovAniso* getCovAniso(Id icov) const;
  CovAniso* getCovAniso(Id icov); // TODO : beurk :(
  void setCov(Id icov, const CovBase* cov) override;
  const ECov& getCovType(Id icov) const override;
  String getCovName(Id icov) const override;
  void setRangeIsotropic(Id icov, double range);
  void setType(Id icov, const ECov& type);
  void setParam(Id icov, double value);
  void setMarkovCoeffs(Id icov, const VectorDouble& coeffs);
  double getParam(Id icov) const;
  double getRange(Id icov) const;
  VectorDouble getRanges(Id icov) const;
  VectorDouble getAngles(Id icov) const;
  Id getNGradParam(Id icov) const;
  CovAniso extractCova(Id icov) const;
  Id getCovMinIRFOrder() const;
  double getBallRadius() const;
  Id hasExternalCov() const;
  bool isChangeSupportDefined() const;
  void appendParams(ListParams& listParams,
                    std::vector<covmaptype>* gradFuncs = nullptr) override;
  // Methods necessary for Optimization
  bool hasNugget() const;
  Id getRankNugget() const;
  const CovAnisoList* createReduce(const VectorInt& validVars) const;

  // Non-stationary parameters
  void makeRangeNoStatDb(Id icov, const String& namecol, Id idim = 0);
  void makeScaleNoStatDb(Id icov, const String& namecol, Id idim = 0);
  void makeAngleNoStatDb(Id icov, const String& namecol, Id idim = 0);

  void makeTensorNoStatDb(Id icov, const String& namecol, Id idim = 0, Id jdim = 0);
  void makeParamNoStatDb(Id icov, const String& namecol);
  void makeRangeNoStatFunctional(Id icov, const AFunctional* func, Id idim = 0);
  void makeScaleNoStatFunctional(Id icov, const AFunctional* func, Id idim = 0);
  void makeAngleNoStatFunctional(Id icov, const AFunctional* func, Id idim = 0);
  void makeTensorNoStatFunctional(Id icov, const AFunctional* func, Id idim = 0, Id jdim = 0);
  void makeParamNoStatFunctional(Id icov, const AFunctional* func);
  void makeRangeStationary(Id icov, Id idim = 0);
  void makeScaleStationary(Id icov, Id idim = 0);
  void makeAngleStationary(Id icov, Id idim = 0);

  void makeTensorStationary(Id icov, Id idim, Id jdim);
  void makeParamStationary(Id icov);

  bool getSameRotation() const { return _sameRotation; }
  void setSameRotation(bool samerot) { _sameRotation = samerot; }

private:
  // Returns a pointer on an existing Cov and cast it to CovAniso
  const CovAniso* _getCovAniso(Id icov) const;
  CovAniso* _getCovAnisoModify(Id icov);

protected:
  bool _isCovarianceIndexValid(Id icov) const;

private:
  bool _sameRotation;
};
} // namespace gstlrn