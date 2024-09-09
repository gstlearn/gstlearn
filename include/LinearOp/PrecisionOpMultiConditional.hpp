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

#ifndef SWIG
  #include <Eigen/src/Core/Matrix.h>
#endif

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
  
  std::pair<double,double> computeRangeEigenVal() const;
  std::pair<double,double> rangeEigenValQ() const;
  double getMaxEigenValProj() const;
  double sumLogVar() const;

  double computeLogDetQ(int nbsimu = 1, int seed = 123) const;
  double computeTotalLogDet(int nbsimu = 1, int seed = 123) const;
  void preparePoly(Chebychev& logPoly) const;
  
  const ProjMatrix* getProjMatrix(int i = 0) const { return (ProjMatrix*) _multiProjData[i];}
  const PrecisionOp* getMultiPrecisionOp(int i = 0) const { return _multiPrecisionOp[i]; }

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

  VectorDouble computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const;

protected:
  void _allocate(int i) const;

#ifndef SWIG
  private:
    void _AtA(const std::vector<Eigen::VectorXd>& inv, std::vector<Eigen::VectorXd>& outv) const;
    void _evalDirect(const std::vector<Eigen::VectorXd>& inv, std::vector<Eigen::VectorXd>& outv) const override;

  public:  
    std::vector<Eigen::VectorXd> computeRhs(const Eigen::VectorXd& datVal) const;
    void computeRhsInPlace(const Eigen::VectorXd& datVal,std::vector<Eigen::VectorXd>& rhs) const;
    void simulateOnMeshings(std::vector<Eigen::VectorXd> &result) const;
    void simulateOnMeshing(Eigen::VectorXd& result,int icov = 0) const;
    void simulateOnDataPointFromMeshings(const std::vector<Eigen::VectorXd>& simus,Eigen::VectorXd& result) const;
    void evalInvCov(const Eigen::VectorXd& inv, Eigen::VectorXd& result) const;
    double computeQuadratic(const Eigen::VectorXd& x) const;

#endif
public : 
  VectorVectorDouble computeRhs(const VectorDouble &datVal) const;

private:
  mutable std::vector<PrecisionOp*>    _multiPrecisionOp; // Pointers are simply stored; do not delete
  std::vector<IProjMatrix*>            _multiProjData; // Pointers are simply stored; do not delete
  VectorDouble                         _varianceData; // Dimension: _ndat
  int                                  _ndat;
  int                                  _ncova;
  mutable Eigen::VectorXd              _work1;
  mutable Eigen::VectorXd              _work1bis;
  mutable Eigen::VectorXd              _work1ter;
  mutable Eigen::VectorXd              _workdata;
  mutable std::vector<Eigen::VectorXd> _work2;
  mutable std::vector<Eigen::VectorXd> _work3;
};
