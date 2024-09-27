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
  virtual int push_back(PrecisionOp *pmatElem, IProjMatrix *projDataElem = nullptr);
  virtual double computeLogDetOp(int nbsimu = 1) const;

  /// Interface for ALinearOpMulti
  int  sizes() const override { return static_cast<int> (_multiPrecisionOp.size()); }
  int  size(int i) const override { return _multiPrecisionOp[i]->getSize(); }

  VectorDouble getAllVarianceData() const {return _varianceData;}
  double getVarianceData(int iech)const {return  _varianceData[iech];}
  void setVarianceData(double nugg){ _varianceData = VectorDouble(_ndat,nugg);}
  void setVarianceDataVector(const VectorDouble& nugg){_varianceData = nugg;}
  
  std::pair<double,double> computeRangeEigenVal() const;
  std::pair<double,double> rangeEigenValQ() const;
  double getMaxEigenValProj() const;
  double sumLogVar() const;

  double computeLogDetQ(int nbsimu = 1) const;
  double computeTotalLogDet(int nbsimu = 1) const;
  void preparePoly(Chebychev& logPoly) const;
  
  const ProjMatrix* getProjMatrix(int i = 0) const { return (ProjMatrix*) _multiProjData[i];}
  const PrecisionOp* getMultiPrecisionOp(int i = 0) const { return _multiPrecisionOp[i]; }

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

  VectorDouble computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const;

protected:
  void _allocate(int i) const;

#ifndef SWIG
  private:
    void _AtA(const std::vector<std::vector<double>>& inv, std::vector<std::vector<double>>& outv) const;
    void _evalDirect(const std::vector<std::vector<double>>& inv, std::vector<std::vector<double>>& outv) const override;

  public:  
    std::vector<std::vector<double>> computeRhs(const std::vector<double>& datVal) const;
    void computeRhsInPlace(const std::vector<double>& datVal,std::vector<std::vector<double>>& rhs) const;
    void simulateOnMeshings(std::vector<std::vector<double>> &result) const;
    void simulateOnMeshing(std::vector<double>& result,int icov = 0) const;
    void simulateOnDataPointFromMeshings(const std::vector<std::vector<double>>& simus,std::vector<double>& result) const;
    void evalInvCov(const constvect inv, std::vector<double>& result) const;
    double computeQuadratic(const std::vector<double>& x) const;

#endif

private:
  mutable std::vector<PrecisionOp*>        _multiPrecisionOp; // Pointers are simply stored; do not delete
  std::vector<IProjMatrix*>                _multiProjData; // Pointers are simply stored; do not delete
  VectorDouble                             _varianceData; // Dimension: _ndat
  int                                      _ndat;
  int                                      _ncova;
  mutable std::vector<double>              _work1;
  mutable std::vector<double>              _work1bis;
  mutable std::vector<double>              _work1ter;
  mutable std::vector<double>              _workdata;
  mutable std::vector<std::vector<double>> _work2;
  mutable std::vector<std::vector<double>> _work3;
};
