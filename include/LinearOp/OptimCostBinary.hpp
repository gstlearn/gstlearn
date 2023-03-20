/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "IOptimCost.hpp"

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
  /*!  Set the constant parameters for internal Conjugate Gradient */
  void setCGParams(int cgmaxiter = 100, double cgeps = 1.e-08)
  {
    _cgMaxIter = cgmaxiter;
    _cgEps = cgeps;
  }
  /*!  Set the constant parameters for internal Pre-Conditioner */
  void setPreCondParams(int chebncmax = 10001, double chebtol = 5.e-3)
  {
    _flagCgPreCond = true;
    _chebNcmax = chebncmax;
    _chebTol = chebtol;
  }
  /*!  Checks if the Cost Function Optimization has been initialized */
  int isInitialized()
  {
    return _isInitialized;
  }
  int getNPoint() const;
  int getNVertex() const;
  void toggleSeismic(bool status);

private:
  double _evaluateCost(const VectorDouble& indic, const VectorDouble& lambda);
  void _evaluateGrad(const VectorDouble& indic,
                     const VectorDouble& lambda,
                     double* normgrad);
  void _contributeSeismic(const VectorDouble& lambda);
  void _contributeSeismicDerivative(const VectorDouble& lambda);

protected:

private:
  bool _isInitialized;
  bool _flagSeismic;
  double _meanPropRaw;
  double _meanPropGaus;
  PrecisionOp* _pMat;
  const ProjMatrix*  _projData;
  const ProjMatrix*  _projSeis;
  VectorDouble _propSeis;
  VectorDouble _varSeis;

  // Parameters for Conjugate Gradient
  int    _cgMaxIter;
  double _cgEps;

  // Parameters for Preconditionner (optional)
  bool   _flagCgPreCond;
  int    _chebNcmax;
  double _chebTol;

  mutable VectorDouble _grad;
  mutable VectorDouble _workp;   /* Dimension: Npoint  */
  mutable VectorDouble _workx;   /* Dimension: Npoint  */
  mutable VectorDouble _workv;   /* Dimension: Nvertex */
  mutable VectorDouble _lambdav; /* Dimension: Nvertex */
  mutable VectorDouble _works;   /* Dimension: Nseis   */
};
