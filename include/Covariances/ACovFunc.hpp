/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ECov.hpp"

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
  virtual bool         isConsistent() const;
  virtual double       getScadef()    const { return 1; }
  virtual double       getParMax()    const { return 0; }
  virtual unsigned int getMaxNDim()   const { return MAX_INT; } // No Space Dimension limit
  virtual int          getMinOrder()  const { return -1; } // Valid for FAST
  virtual bool         hasInt1D()     const;
  virtual bool         hasInt2D()     const;
  virtual int          hasRange()     const { return 1 ; } // 0:No; 1:Yes; -1:from Sill
  virtual bool         hasParam()     const { return false; }
  virtual String       getCovName()   const = 0;
  virtual bool         hasCovDerivative() const { return false; }
  virtual bool         hasCovOnSphere() const { return false; }
  virtual bool         hasSpectrum() const { return false; }
  ///////////////////////////////////////////////////

  void setParam(double param);
  void setField(double field);
  double evalCov(double h) const;
  double evalCovDerivative(int degree, double h) const;
  double evalCovOnSphere(double alpha, double scale = 1., int degree = 50) const;
  VectorDouble evalCovVec(const VectorDouble& vech) const;
  VectorDouble evalCovDerivativeVec(int degree, const VectorDouble& vech) const;
  const ECov&          getType()    const { return _type; }
  const CovContext&    getContext() const { return _ctxt; }
  double               getParam()   const { return _param; }

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }

protected:
  /// TODO : Gneiting (spatio-temporal covariance) :
  /// Change argument : double h becomes VectorDouble (number of sub-space)
  virtual double _evaluateCov(double h) const = 0;
  virtual double _evaluateCovDerivate(int degree, double h) const;
  virtual double _evaluateCovOnSphere(double scale = 1., int degree = 50) const;
  virtual double _evaluateSpectrum(double freq, double scale, int ndim) const;

private:
  ECov        _type;    /*! Covariance function type */
  CovContext  _ctxt;    /*! Context (space, irfDegree, field, ...) */
  double      _param;   /*! Third parameter (TEST if not used) */
};

