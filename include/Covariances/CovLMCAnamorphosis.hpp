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

#include "Covariances/ACovAnisoList.hpp"
#include "gstlearn_export.hpp"

#include "Enum/EAnam.hpp"
#include "Anamorphosis/AAnam.hpp"

class ASpace;
class SpacePoint;
class CovAniso;
class Model;

class GSTLEARN_EXPORT CovLMCAnamorphosis : public ACovAnisoList
{
public:
  CovLMCAnamorphosis(const AAnam* anam,
                     const VectorInt& strcnt = VectorInt(),
                     const ASpace* space = nullptr);
  CovLMCAnamorphosis(const ACovAnisoList& lmc,
                     const AAnam* anam,
                     const VectorInt& strcnt);
  CovLMCAnamorphosis(const CovLMCAnamorphosis &r);
  CovLMCAnamorphosis& operator= (const CovLMCAnamorphosis &r);
  virtual ~CovLMCAnamorphosis();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMCAnamorphosis)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ACov Interface
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;                   
  void addCov(const CovAniso* cov) override;
  bool hasAnam() const override { return true; }
  const AAnam* getAnam() const override { return _anam; }
  void setActiveFactor(int iclass) override;
  int getActiveFactor() const override { return _activeFactor; }
  int getAnamNClass() const override { return _anam->getNClass(); }

  int init(const VectorInt& strcnt = VectorInt());
  const EAnam getAnamType() const;
  void setAnam(const AAnam*& anam) { _anam = anam; }

protected:
    void _loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const override;

    void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                        const SpacePoint& pwork1, 
                        const SpacePoint& pwork2, 
                        const CovCalcMode *mode) const override;
    void _optimizationSetTarget(const SpacePoint &pt) const override
    {
      ACov::_optimizationSetTarget(pt);
    }

private:
  
  double _evalHermite(int ivar,
                      int jvar,
                      double h,
                      const CovCalcMode* mode) const;
  double _evalDiscreteDD(int ivar,
                         int jvar,
                         double h,
                         const CovCalcMode* mode) const;
  double _evalDiscreteIR(int ivar,
                         int jvar,
                         double h,
                         const CovCalcMode* mode) const;

  double _evalHermite(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode* mode) const;
  double _evalDiscreteDD(int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode* mode) const;
  double _evalDiscreteIR(int ivar,
                         int jvar,
                         const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode* mode) const;
  double _evalHermite0(int ivar, int jvar, const CovCalcMode* mode) const;
  double _evalDiscreteDD0(int ivar, int jvar, const CovCalcMode* mode) const;
  double _evalDiscreteIR0(int ivar, int jvar, const CovCalcMode* mode) const;
  void   _transformCovCalcModeIR(CovCalcMode* mode, int iclass) const;

private:
  int       _activeFactor;    /* Target factor (-1: Raw; 1: Gaussian; n: rank of factor) */
  VectorInt _anamStrCount;    /* List of covariances in the Model (for RI only) */
  const AAnam* _anam;
};
