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

#include "gstlearn_export.hpp"
#include "LinearOp/PrecisionOp.hpp"

class AMesh;
class ShiftOpCs;
class CovAniso;
class Model;

/** This class is just a specialization of PrecisionOp when the shift
* Operator is built with sparse (cs) matrices.
* It allows to return the precision matrix as a cs. */

class GSTLEARN_EXPORT PrecisionOpCs : public PrecisionOp
{
public:
  PrecisionOpCs(ShiftOpCs* shiftop = nullptr,
                const CovAniso* cova = nullptr,
                const EPowerPT& power = EPowerPT::fromKey("UNDEFINED"),
                bool verbose = false);

  PrecisionOpCs(AMesh* mesh,
                Model* model,
                int igrf,
                const EPowerPT& power,
                bool verbose);


  void evalDeriv(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void evalDerivOptim(VectorDouble& outv,int iapex,int igparam) override;
  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result) override;
  void gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result) override;
  virtual ~PrecisionOpCs();
  VectorDouble getCoeffs();
  cs* getQ();
};
