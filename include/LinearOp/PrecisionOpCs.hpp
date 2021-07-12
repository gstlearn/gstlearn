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

// This class is just a specialisation of PrecisionOp when the shift Operator is built with cs matrices.
// It allows to return the precision matrix as a cs.

#pragma once
#include "LinearOp/PrecisionOp.hpp"
#include "geoslib_enum.h"

class ShiftOpCs;
class Model;

class PrecisionOpCs : public PrecisionOp
{
public:
  PrecisionOpCs(ShiftOpCs* shiftop = nullptr,
                const Model* model = nullptr,
                int icov = 0,
                ENUM_POPTS power = POPT_UNDEFINED,
                bool verbose = false);
  void evalDeriv(const VectorDouble& in, VectorDouble& out,int iapex,int igparam) override;
  void evalDerivPoly(const VectorDouble& in, VectorDouble& out,int iapex,int igparam) override;
  void gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble result) override;
  virtual ~PrecisionOpCs();
  VectorDouble getCoeffs();
  cs* getQ();
};


