/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
  PrecisionOpCs(const AMesh* mesh,
                Model* model,
                int icov = 0,
                const EPowerPT& power = EPowerPT::fromKey("ONE"),
                bool verbose = false);
  virtual ~PrecisionOpCs();

  void evalDeriv(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void evalDerivOptim(VectorDouble& outv,int iapex,int igparam) override;
  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result) override;
  void gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result) override;
  VectorDouble getCoeffs();
  const cs* getQ() const { return _Q; }

private:
  void _buildQ();

private:
  cs* _Q;
};
