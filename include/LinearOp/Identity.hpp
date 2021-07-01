/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "LinearOp/ALinearOp.hpp"
#include "Basic/Vector.hpp"

class Identity: public ALinearOp
{

public:
  Identity(int n);
  virtual ~Identity();

  void evalInverse(const VectorDouble& in, VectorDouble& out) const override;
  int getSize() const override { return _n; }

protected:
  void _evalDirect(const VectorDouble& in, VectorDouble& out) const override;

private:
  int _n;
};
