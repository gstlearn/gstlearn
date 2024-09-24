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

#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/IOptimCost.hpp"

class PrecisionOp;
class ProjMatrix;

class GSTLEARN_EXPORT OptimCostBinary: public IOptimCost
{

public:
  OptimCostBinary();
  OptimCostBinary(const OptimCostBinary &m);
  OptimCostBinary& operator = (const OptimCostBinary &m);
  virtual ~OptimCostBinary();

  void reset(PrecisionOp* pmat,
             const ProjMatrix* projdata,
             const ProjMatrix* projseis = nullptr,
             const VectorDouble& propseis = VectorDouble(),
             const VectorDouble& varseis = VectorDouble());
  VectorDouble minimize(VectorDouble& indic,
                        bool verbose = false,
                        int maxiter = 100,
                        double eps = 5.e-4);
  void calculateGradient(const VectorDouble& indic,
                         const VectorDouble& lambda,
                         double* out);
  int setMeanProportion(double meanprop);
  /*!  Set the constant parameters for internal Pre-Conditioner */
  static void setPreCondParams(int chebncmax = 10001, double chebtol = 5.e-3)
  {
    DECLARE_UNUSED(chebncmax, chebtol);
  }
  int isInitialized() const { return _isInitialized; }
  int getNPoint() const;
  int getNVertex() const;
  void toggleSeismic(bool status);

private:
  double _evaluateCost(const VectorDouble& indic, const VectorDouble& lambda);
  void _evaluateGrad(const VectorDouble &indic,
                     const VectorDouble &lambda,
                     double *normgrad);
  void _contributeSeismic(const VectorDouble& lambda);
  void _contributeSeismicDerivative(const VectorDouble& lambda);

protected:

private:
  bool               _isInitialized;
  bool               _flagSeismic;
  double             _meanPropRaw;
  double             _meanPropGaus;
  PrecisionOp*       _pMat;
  const ProjMatrix*  _projData;
  const ProjMatrix*  _projSeis;
  VectorDouble       _propSeis;
  VectorDouble       _varSeis;

  mutable VectorDouble _grad;
  mutable VectorDouble _workp;   /* Dimension: Npoint  */
  mutable VectorDouble _workx;   /* Dimension: Npoint  */
  mutable VectorDouble _workv;   /* Dimension: Nvertex */
  mutable VectorDouble _lambdav; /* Dimension: Nvertex */
  mutable VectorDouble _works;   /* Dimension: Nseis   */
};
