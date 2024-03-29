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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ALinearOpMulti.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"

#include <vector>

class Chebychev;
/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMultiConditional : public ALinearOpMulti {

public:
  PrecisionOpMultiConditional();
  PrecisionOpMultiConditional(const PrecisionOpMultiConditional &m)= delete;
  PrecisionOpMultiConditional& operator= (const PrecisionOpMultiConditional &m)= delete;
  virtual ~PrecisionOpMultiConditional();

  /// Interface for PrecisionOpMultiConditional
  virtual void makeReady();
  virtual int push_back(PrecisionOp *pmatElem, IProjMatrix *projDataElem = nullptr);
  virtual double computeLogDetOp(int nbsimu = 1, int seed = 123) const;

  /// Interface for ALinearOpMulti
  int  sizes() const override { return static_cast<int> (_multiPrecisionOp.size()); }
  int  size(int i) const override { return _multiPrecisionOp[i]->getSize(); }

  VectorDouble getAllVarianceData() const {return _varianceData;}
  double getVarianceData(int iech)const {return  _varianceData[iech];}
  void setVarianceData(double nugg){ _varianceData = VectorDouble(_ndat,nugg);}
  void setVarianceDataVector(const VectorDouble& nugg){_varianceData = nugg;}
  VectorVectorDouble computeRhs(const VectorDouble& datVal) const;
  void computeRhsInPlace(const VectorDouble& datVal,VectorVectorDouble& rhs) const;
  void simulateOnMeshings(VectorVectorDouble &result) const;
  void simulateOnMeshing(VectorDouble& result,int icov = 0) const;
  void simulateOnDataPointFromMeshings(const VectorVectorDouble& simus,VectorDouble& result) const;
  void evalInvCov(const VectorDouble& inv, VectorDouble& result) const;
  std::pair<double,double> computeRangeEigenVal() const;
  std::pair<double,double> rangeEigenValQ() const;
  double getMaxEigenValProj() const;
  double sumLogVar() const;

  double computeLogDetQ(int nbsimu = 1, int seed = 123) const;
  double computeTotalLogDet(int nbsimu = 1, int seed = 123) const;
  double computeQuadratic(const VectorDouble& x) const;
  void preparePoly(Chebychev& logPoly) const;
  void AtA(const VectorVectorDouble& inv,VectorVectorDouble& outv) const;
  VectorDouble computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const;
  const ProjMatrix* getProjMatrix(int i = 0) const { return (ProjMatrix*) _multiProjData[i];}
  const PrecisionOp* getMultiPrecisionOp(int i = 0) const { return _multiPrecisionOp[i]; }

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

protected:
  void _evalDirect(const VectorVectorDouble& inv, VectorVectorDouble& outv) const override;
  void _allocate(int i) const;

private:
  std::vector<PrecisionOp*>  _multiPrecisionOp; // Pointers are simply stored; do not delete
  std::vector<IProjMatrix*>  _multiProjData; // Pointers are simply stored; do not delete
  VectorDouble               _varianceData; // Dimension: _ndat
  int                        _ndat;
  int                        _ncova;
  mutable VectorDouble       _work1;
  mutable VectorDouble       _work1bis;
  mutable VectorDouble       _work1ter;
  mutable VectorDouble       _workdata;
  mutable VectorVectorDouble _work2;
  mutable VectorVectorDouble _work3;
};
