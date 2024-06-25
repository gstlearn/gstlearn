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

#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

#include "Basic/AStringable.hpp"
#include "Covariances/CovContext.hpp"
#include "Arrays/Array.hpp"
#include "Matrix/MatrixRectangular.hpp"

class TurningBandOperate;

/* Covariance basic function for normalized sill and distance:
 * Positive definite function
 * */

class GSTLEARN_EXPORT ACovFunc : public AStringable
{
public:
  ACovFunc(const ECov& type, const CovContext& ctxt);
  ACovFunc(const ACovFunc &r);
  ACovFunc& operator= (const ACovFunc &r);
  virtual ~ACovFunc();

  ///////////////////////////////////////////////////
  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ACovFunc Interface
  virtual String       getFormula()   const { return String("Equation not yet implemented"); }
  virtual double       getScadef()    const { return 1; }
  virtual double       getParMax()    const { return 0; }
  virtual bool         hasInt1D()     const;
  virtual bool         hasInt2D()     const;
  virtual int          hasRange()     const { return 1 ; } // 0:No; 1:Yes; -1:from Sill
  virtual bool         hasParam()     const { return false; }
  virtual String       getCovName()   const = 0;
  virtual bool         hasCovDerivative()    const { return false; }
  virtual bool         hasCovOnSphere()      const { return false; }
  virtual bool         hasSpectrumOnSphere() const { return false; }
  virtual bool         hasSpectrum()         const { return false; }
  virtual bool         hasMarkovCoeffs()     const { return false; }

  virtual bool         isConsistent() const;
  virtual unsigned int getMaxNDim()   const { return MAX_INT; } // No Space Dimension limit
  virtual int          getMinOrder()  const { return -1; } // Valid for FAST
  virtual bool         getCompatibleSpaceR() const { return false; }
  virtual bool         getCompatibleSpaceS() const { return false; }

  // Specific to Turning Band Simulation Method
  virtual bool isValidForTurningBand() const { return false; }
  virtual double simulateTurningBand(double t0, TurningBandOperate &operTB) const
  {
    DECLARE_UNUSED(t0);
    DECLARE_UNUSED(operTB);
    return TEST;
  }

  // Specific for Spectral Simulation Method
  virtual bool isValidForSpectral() const { return false; }
  virtual MatrixRectangular simulateSpectralOmega(int nb) const
  {
    DECLARE_UNUSED(nb);
    return MatrixRectangular();
  }

  ///////////////////////////////////////////////////

  void setParam(double param);
  void setField(double field);
  double evalCov(double h) const;
  double evalCovDerivative(int degree, double h) const;
  double evalCovOnSphere(double alpha, double scale = 1., int degree = 50) const;
  VectorDouble evalSpectrumOnSphere(int n, double scale = 1.) const;
  VectorDouble evalCovVec(const VectorDouble& vech) const;
  VectorDouble evalCovDerivativeVec(int degree, const VectorDouble& vech) const;
  const ECov&          getType()    const { return _type; }
  const CovContext&    getContext() const { return _ctxt; }
  double               getParam()   const { return _param; }

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }
  virtual double evaluateSpectrum(double freq, int ndim) const;
  virtual VectorDouble getMarkovCoeffs() const;
  virtual void setMarkovCoeffs(VectorDouble coeffs);
  virtual double getCorrec() const {return 1.;}
  virtual void setCorrec(double val)
  {
    DECLARE_UNUSED(val);
  }
  virtual void computeCorrec(int ndim);
  virtual void computeMarkovCoeffs(int dim)
  {
    DECLARE_UNUSED(dim);
  }

protected:
  /// TODO : Gneiting (spatio-temporal covariance) :
  /// Change argument : double h becomes VectorDouble (number of sub-space)
  virtual double _evaluateCov(double h) const
  {
    DECLARE_UNUSED(h);
    return TEST;
  }
  ;
  virtual double _evaluateCovDerivative(int degree, double h) const;
  virtual double _evaluateCovOnSphere(double alpha,
                                      double scale = 1.,
                                      int degree = 50) const;
  virtual VectorDouble _evaluateSpectrumOnSphere(int n,
                                                 double scale = 1.,
                                                 double param = 1.) const;

private:
  Array _evalCovFFT(const VectorDouble& ext, int N = 128) const;
  ECov        _type;    /*! Covariance function type */
  CovContext  _ctxt;    /*! Context (space, number of variables, ...) */
  double      _param;   /*! Third parameter (TEST if not used) */
};

