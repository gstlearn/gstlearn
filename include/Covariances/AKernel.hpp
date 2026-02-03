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

#include "Basic/AStringFormat.hpp"
#include "Enum/ESimuType.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

#include "Arrays/Array.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"

namespace gstlrn
{

class TurningBandOperate;

/* Covariance basic function for normalized sill and distance:
 * Positive definite function
 * */
class GSTLEARN_EXPORT AKernel: public AStringable
{
public:
  AKernel(const ECov& type, const CovContext& ctxt);
  AKernel(const AKernel& r);
  AKernel& operator=(const AKernel& r);
  virtual ~AKernel();

  ///////////////////////////////////////////////////
  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  ///////////////////////////////////////////////////
  /// AKernel Interface
  virtual String getFormula() const { return String("Equation is not provided"); }
  virtual double getScadef() const { return 1; }
  virtual double getParMax() const { return 0; }
  virtual bool hasInt1D() const;
  virtual bool hasInt2D() const;
  virtual Id hasRange() const { return 1; } // 0:No; 1:Yes; -1:from Sill
  virtual bool hasParam() const { return false; }
  virtual String getCovName() const = 0;
  virtual bool hasCovDerivative() const { return false; }
  virtual bool hasMarkovCoeffs() const { return false; }
  virtual double normalizeOnSphere(Id n = 50, double scale = 1.) const
  {
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(scale);
    return 1.;
  }
  virtual bool isConsistent() const;
  virtual size_t getMaxNDim() const { return MAX_INT; } // No Space Dimension limit
  virtual Id getMinOrder() const { return -1; }         // Valid starting from FAST
  virtual bool getCompatibleSpaceR() const { return false; }
  virtual bool getCompatibleSpaceS() const { return false; }

  // Specific to Turning Band Simulation Method
  virtual double simulateTurningBand(double t0, TurningBandOperate& operTB) const
  {
    DECLARE_UNUSED(t0);
    DECLARE_UNUSED(operTB);
    return TEST;
  }

  // Specific for Spectral Simulation Method
  virtual bool isValidForSimulation(const ESimuType& simuType) const
  {
    DECLARE_UNUSED(simuType);
    return false;
  }
  virtual MatrixDense simulateSpectralOmega(Id nb) const
  {
    DECLARE_UNUSED(nb);
    return MatrixDense();
  }

  static AKernel* createFromType(const ECov& type, Id ndim = 1);

  virtual void evalCorFuncBatch(constvect h, vect res) const;
  ///////////////////////////////////////////////////

  void setParam(double param, Id ipar = 0);
  void setParams(const VectorDouble& params);
  void setField(double field);
  void setContext(const CovContext& ctxt) { _ctxt = ctxt; }
  double evalCorFunc(double h) const;
  double evalCovDerivative(Id degree, double h) const;
  double evalCovFirstDerivativeOverH(double h) const;

  const ECov& getType() const { return _type; }
  const CovContext& getContext() const { return _ctxt; }
  ESpaceType getSpaceType() const { return _ctxt.getSpace()->getType(); }
  Id getNParams() const { return static_cast<Id>(_params.size()); }
  double getParam(Id ipar = 0) const;
  const VectorDouble& getParams() const { return _params; }

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }
  virtual double evaluateSpectrum(double freq) const;
  virtual VectorDouble getMarkovCoeffs() const;
  virtual void setMarkovCoeffs(const VectorDouble& coeffs);
  virtual void computeMarkovCoeffs(Id dim) { DECLARE_UNUSED(dim); }

  double evaluateCovOnSphere(double alpha, double scale, Id degree = 50) const;
  VectorDouble evaluateSpectrumOnSphere(Id n, double scale, bool flagScale = true) const;

  virtual double getCorrec() const { return 1.; }
  virtual void setCorrec(double val) { DECLARE_UNUSED(val); }
  virtual bool needCorrec() const { return false; }
  virtual void computeCorrec(Id ndim);

  double evalDerivative(double h) const { return _evaluateCovFirstDerivative(h); }

protected:
  /// TODO : Gneiting (spatio-temporal covariance) :
  /// Change argument : double h becomes VectorDouble (number of sub-space)
  virtual void _setParam(double param, Id ipar = 0)
  {
    DECLARE_UNUSED(param);
    DECLARE_UNUSED(ipar);
  }
  virtual double _evaluateCov(double h) const
  {
    DECLARE_UNUSED(h);
    return TEST;
  };
  virtual double _evaluateCovFirstDerivative(double h) const
  {
    double eps = EPSILON4;
    return (_evaluateCov(h + eps) - _evaluateCov(h - eps)) / (2. * eps);
  }
  virtual double _evaluateCovDerivative(Id degree, double h) const;
  virtual double _evaluateCovFirstDerivativeOverH(double h) const
  {
    double eps = EPSILON4;
    return (_evaluateCov(h + eps) - _evaluateCov(h - eps)) / (2. * eps) / h;
  }
  virtual double _evaluateCovOnSphere(double alpha, double scale = 1., Id degree = 50) const;
  virtual VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true) const;

private:
  Array _evalCovFFT(const VectorDouble& hmax, Id N = 128) const;
  ECov _type;           /*! Covariance function type */
  CovContext _ctxt;     /*! Context (space, number of variables, ...) */
  VectorDouble _params; /*! Parameters vector (empty if not used) */
};

GSTLEARN_EXPORT VectorInt getAllMaxDimension();
GSTLEARN_EXPORT VectorBool getAllIsSpatialTemporal();
GSTLEARN_EXPORT VectorString getAllFormula();

} // namespace gstlrn
