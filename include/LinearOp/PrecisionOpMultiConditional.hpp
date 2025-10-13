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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ALinearOpMulti.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include <vector>


namespace gstlrn
{
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
  virtual void makeReady(){};
  virtual Id push_back(PrecisionOp *pmatElem, IProj *projDataElem = nullptr);
  virtual double computeLogDetOp(Id nbsimu = 1) const;

  /// Interface for ALinearOpMulti
  Id  sizes() const override { return static_cast<Id> (_multiPrecisionOp.size()); }
  Id  size(Id i) const override { return _multiPrecisionOp[i]->getSize(); }

  VectorDouble getAllVarianceData() const {return _varianceData;}
  double getVarianceData(Id iech)const {return  _varianceData[iech];}
  void setVarianceData(double nugg){ _varianceData = VectorDouble(_ndat,nugg);}
  void setVarianceDataVector(const VectorDouble& nugg){_varianceData = nugg;}
  
  std::pair<double,double> computeRangeEigenVal() const;
  std::pair<double,double> rangeEigenValQ() const;
  double getMaxEigenValProj() const;
  double computeLogDetNoise() const;

  double computeLogDetQ(Id nMC = 1) const;
  double computeTotalLogDet(Id nMC = 1, bool verbose = false, Id seed = 13132) const;
  void preparePoly(Chebychev& logPoly) const;
  
  const ProjMatrix* getProjMatrix(Id i = 0) const { return (ProjMatrix*) _multiProjData[i];}
  const PrecisionOp* getMultiPrecisionOp(Id i = 0) const { return _multiPrecisionOp[i]; }

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

  VectorDouble computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const;

protected:
  void _allocate(Id i) const;

#ifndef SWIG
  private:
    void _AtA(const std::vector<std::vector<double>>& inv, std::vector<std::vector<double>>& outv) const;
    void _evalDirect(const std::vector<std::vector<double>>& inv, std::vector<std::vector<double>>& outv) const override;

  public:  
    std::vector<std::vector<double>> computeRhs(const std::vector<double>& datVal) const;
    void computeRhsInPlace(const std::vector<double>& datVal,std::vector<std::vector<double>>& rhs) const;
    void simulateOnMeshings(std::vector<std::vector<double>> &result) const;
    void simulateOnMeshing(std::vector<double>& result,Id icov = 0) const;
    void simulateOnDataPointFromMeshings(const std::vector<std::vector<double>>& simus,std::vector<double>& result) const;
    void evalInvCov(const constvect inv, std::vector<double>& result) const;
    double computeQuadratic(const std::vector<double>& x) const;

#endif

private:
  mutable std::vector<PrecisionOp*>        _multiPrecisionOp; // Pointers are simply stored; do not delete
  std::vector<IProj*>                      _multiProjData; // Pointers are simply stored; do not delete
  VectorDouble                             _varianceData; // Dimension: _ndat
  Id                                      _ndat;
  Id                                      _ncova;
  mutable std::vector<double>              _work1;
  mutable std::vector<double>              _work1bis;
  mutable std::vector<double>              _work1ter;
  mutable std::vector<double>              _workdata;
  mutable std::vector<std::vector<double>> _work2;
  mutable std::vector<std::vector<double>> _work3;
};
}