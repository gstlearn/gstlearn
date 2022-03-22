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
  CovLMCAnamorphosis(const EAnam& anam_type,
                     int anam_nclass,
                     int anam_iclass,
                     int anam_var,
                     double anam_coefr,
                     double anam_coefs,
                     VectorDouble& anam_strcnt,
                     VectorDouble& anam_stats,
                     const ASpace* space = nullptr);
  CovLMCAnamorphosis(const CovLMCAnamorphosis &r);
  CovLMCAnamorphosis& operator= (const CovLMCAnamorphosis &r);
  virtual ~CovLMCAnamorphosis();

  virtual IClonable* clone() const override { return new CovLMCAnamorphosis(*this); };
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(const EAnam& anam_type,
           int anam_nclass,
           int anam_iclass,
           int anam_var,
           double anam_coefr,
           double anam_coefs,
           VectorDouble& anam_strcnt,
           VectorDouble& anam_stats);

  int getAnamIClass() const { return _anamIClass; }
  const VectorDouble& getAnamMeans() const { return _anamMeans; }
  double getAnamMeans(int iclass) const { return _anamMeans[iclass]; }
  int getAnamNClass() const { return _anamNClass; }
  int getAnamPointBlock() const { return _anamPointBlock; }
  const VectorDouble& getAnamStrCount() const { return _anamStrCount; }
  const EAnam getAnamType() const { return _anamType; }

  int setAnamIClass(int anamIClass);
  void setAnamPointBlock(int anamPointBlock) { _anamPointBlock = anamPointBlock; }
  AAnam* getAnam() const { return _anam; }

private:
  double _evalHermite(AnamHermite *anam,
                      int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode) const;
  double _evalDiscreteDD(AnamDiscreteDD *anam,
                         int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode& mode) const;
  double _evalDiscreteIR(AnamDiscreteIR *anam,
                         int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode& mode) const;
  double _st_cov_residual(AnamDiscreteIR *anam,
                          CovCalcMode& mode,
                          int icut0,
                          int ivar,
                          int jvar,
                          const SpacePoint& p1,
                          const SpacePoint& p2) const;
  double _st_covsum_residual(AnamDiscreteIR* anam,
                             CovCalcMode& mode,
                             int icut0,
                             int ivar,
                             int jvar,
                             const SpacePoint& p1,
                             const SpacePoint& p2) const;

private:
  EAnam  _anamType;
  int    _anamIClass;         /* Target factor (-1: discretized grade) */
  int    _anamNClass;         /* Number of indicator classes */
  int    _anamPointBlock;     /* Type of point / block covariance */
  VectorDouble _anamStrCount; /* Array of structure count per model (IR)  */
  VectorDouble _anamMeans;    /* Array of statistics per class */
  AAnam*       _anam;
};
