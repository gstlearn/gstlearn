/*
 * PrecisionOpMultiConditional.hpp
 *
 *  Created on: 18 déc. 2020
 *      Author: ndesassis
 */

#pragma once

#include "Basic/Vector.hpp"
#include "ALinearOpMulti.hpp"
#include "PrecisionOp.hpp"
#include <vector>
#include "ProjMatrix.hpp"

/**
 * Class to store objects for SPDE
 */
class PrecisionOpMultiConditional : public ALinearOpMulti {

public:
  PrecisionOpMultiConditional();

  virtual ~PrecisionOpMultiConditional();

  void push_back(PrecisionOp*  pmatElem,
                 IProjMatrix* projDataElem = nullptr);
  VectorDouble getVarianceData() const {return _varianceData;}
  double getVarianceData(int iech)const {return  _varianceData[iech];}
  void setVarianceData(double nugg){ _varianceData = VectorDouble(_ndat,nugg);}
  void setVarianceData(const VectorDouble& nugg){_varianceData = nugg;}
  /*!  Returns the dimension of the matrix */
  int  size() const override { return static_cast<int> (_multiPrecisionOp.size()); }
  int  size(int i) const override { return _multiPrecisionOp[i]->getSize(); }
  VectorVectorDouble computeRhs(const VectorDouble& datVal) const;
  void computeRhs(const VectorDouble& datVal,VectorVectorDouble& rhs) const;
  void simulateOnMeshing(const VectorDouble& gauss,VectorVectorDouble& result) const;
  void simulateOnDataPointFromMeshings(const VectorVectorDouble& simus,VectorDouble& result) const;
  void evalInvCov(const VectorDouble& in, VectorDouble& result) const;
  VectorDouble computeCoeffs(const VectorDouble& Y, const VectorVectorDouble& X) const;

protected:
  void _evalDirect(const VectorVectorDouble& in,
                   VectorVectorDouble& out) const override;

private:
  std::vector<PrecisionOp*>  _multiPrecisionOp;
  std::vector<IProjMatrix*>  _multiProjData;
  VectorDouble               _varianceData; //dimension: _ndat
  int                        _ndat;
  int                        _ncova;
  mutable VectorDouble       _work1;
  mutable VectorDouble       _work1bis;
  mutable VectorVectorDouble _work2;
  mutable VectorVectorDouble _work3;
};
