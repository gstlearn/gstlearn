/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EAnam.hpp"
#include "Enum/EConvDir.hpp"
#include "Enum/EConvType.hpp"

#include "Covariances/CovLMC.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"

class ASpace;
class SpacePoint;
class CovAniso;
class Model;

class GSTLEARN_EXPORT CovLMCAnamorphosis : public CovLMC
{
public:
  CovLMCAnamorphosis(const AAnam* anam,
                     const VectorInt& strcnt = VectorInt(),
                     const ASpace* space = nullptr);
  CovLMCAnamorphosis(const CovLMC& lmc,
                     const AAnam* anam,
                     const VectorInt& strcnt);
  CovLMCAnamorphosis(const CovLMCAnamorphosis &r);
  CovLMCAnamorphosis& operator= (const CovLMCAnamorphosis &r);
  virtual ~CovLMCAnamorphosis();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMCAnamorphosis)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  /// Interface for ACovAnisoList
  void addCov(const CovAniso* cov) override;
  bool hasAnam() const override { return true; }
  const AAnam* getAnam() override { return _anam; }
  int setActiveFactor(int iclass) override;
  int getActiveFactor() const override { return _activeFactor; }
  int getAnamNClass() const override { return _anam->getNClass(); }

  int init(const VectorInt& strcnt = VectorInt());
  const EAnam getAnamType() const;
  void setAnam(const AAnam*& anam) { _anam = anam; }

private:
  double _evalHermite(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode) const;
  double _evalDiscreteDD(int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode& mode) const;
  double _evalDiscreteIR(int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode& mode) const;
  double _evalHermite0(int ivar, int jvar, const CovCalcMode& mode) const;
  double _evalDiscreteDD0(int ivar, int jvar, const CovCalcMode& mode) const;
  double _evalDiscreteIR0(int ivar, int jvar, const CovCalcMode& mode) const;
  void   _transformCovCalcModeIR(CovCalcMode& mode, int iclass) const;

private:
  int    _activeFactor;       /* Target factor (-1: Raw; 1: Gaussian; n: rank of factor) */
  VectorInt _anamStrCount;    /* List of covariances in the Model (for RI only) */
  const AAnam* _anam;
};
