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

class PrecisionOpMultiConditional : public ALinearOpMulti {

public:
  PrecisionOpMultiConditional();

  virtual ~PrecisionOpMultiConditional();

  void push_back(PrecisionOp*  pmatElem,
                 IProjMatrix* projDataElem);
  void setNugget(double nugg){_nugget = nugg;}
  /*!  Returns the dimension of the matrix */
  int  size() const override { return _pMat.size(); }
  int  size(int i) const override { return _pMat[i]->getSize(); }

protected:
  void _evalDirect(const VectorVectorDouble& in,
                   VectorVectorDouble& out) const override;

private:
  std::vector<PrecisionOp*> _pMat;
  std::vector<IProjMatrix*> _projData;
  double _nugget;
  int    _ndat;
  mutable VectorDouble       _work1;
  mutable VectorVectorDouble _work2;
};
