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
#include "LinearOp/PrecisionOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Basic/Vector.hpp"
#include "Model/Model.hpp"
#include "geoslib_e.h"
#include "csparse_d.h"
#include "LinearOp/ShiftOpCs.hpp"

PrecisionOpCs::PrecisionOpCs(const ShiftOpCs* shiftop,
                             const Model* model,
                             int icov,
                             ENUM_POPTS power,
                             bool verbose)
    : PrecisionOp(shiftop, model, icov, power, verbose)
{
}

PrecisionOpCs::~PrecisionOpCs()
{
  // TODO Auto-generated destructor stub
}

VectorDouble PrecisionOpCs::getCoeffs()
{
  VectorDouble coeffs = getPoly(POPT_ONE)->getCoeffs();
  return coeffs;
}

cs *PrecisionOpCs::getQ()
{
  VectorDouble blin = getPoly(POPT_ONE)->getCoeffs();
  cs* Q = spde_build_Q(getShiftOp()->getS(), getShiftOp()->getLambda(),
                       blin.size(), blin.data());
  return Q;
}
