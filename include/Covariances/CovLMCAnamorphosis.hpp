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

#include "Basic/Vector.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/EConvType.hpp"
#include "Covariances/EConvDir.hpp"
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
                     int anam_iclass,
                     int anam_var,
                     const VectorInt& strcnt = VectorInt(),
                     const ASpace* space = nullptr);

  CovLMCAnamorphosis(const CovLMCAnamorphosis &r);
  CovLMCAnamorphosis& operator= (const CovLMCAnamorphosis &r);
  virtual ~CovLMCAnamorphosis();

  /// Interface from IClonable
  virtual IClonable* clone() const override { return new CovLMCAnamorphosis(*this); };

  /// Interface from ATringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(int anam_iclass, int anam_var, const VectorInt& strcnt);

  int getAnamIClass() const { return _anamIClass; }
  int getAnamPointBlock() const { return _anamPointBlock; }
  const EAnam getAnamType() const;

  int setAnamIClass(int anamIClass);
  void setAnamPointBlock(int anamPointBlock) { _anamPointBlock = anamPointBlock; }

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
  double _covSumResidualIR(const VectorDouble& covs, int icut0) const;

private:
  int    _anamIClass;         /* Target factor (-1: discretized grade) */
  int    _anamPointBlock;     /* Type of point / block covariance */
  VectorInt _anamStrCount;    /* List of covariances in the Model (for RI only) */
  const AAnam* _anam;
};
