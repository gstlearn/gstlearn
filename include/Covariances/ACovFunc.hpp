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

#include "Basic/Vector.hpp"

#include "Basic/AStringable.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ECov.hpp"
#include "geoslib_enum.h"

/* Covariance basic function for normalized sill and distance:
 * Positive definite function
 * */

class ACovFunc : public AStringable
{
public:
  ACovFunc(const ECov& type, const CovContext& ctxt);
  ACovFunc(const ACovFunc &r);
  ACovFunc& operator= (const ACovFunc &r);
  virtual ~ACovFunc();

  /// Set the third parameter value
  void setParam(double param);
  /// Scale the context field (used for updating to normalized field)
  void setField(double field);

  double evalCov(double h) const;
  double evalCovDerivative(int degree, double h) const;
  VectorDouble evalCovVec(const VectorDouble& vech) const;
  VectorDouble evalCovDerivativeVec(int degree, const VectorDouble& vech) const;

  virtual std::string toString(int level = 0) const override;

  virtual String getFormula()   const { return String("Equation not yet implemented"); }
  virtual bool   isConsistent() const;

  const ECov&          getType()    const { return _type; }
  const CovContext&    getContext() const { return _ctxt; }
  double               getParam()   const { return _param; }

  virtual double       getScadef()    const { return 1; }
  virtual double       getParMax()    const { return 0; }
  virtual unsigned int getMaxNDim()   const { return MAX_INT; } // No Space Dimension limit
  virtual unsigned int getMinOrder()  const { return -1; } // Valid for FAST
  virtual bool         hasInt1D()     const;
  virtual bool         hasInt2D()     const;
  virtual int          hasRange()     const { return 1 ; } // 0:No; 1:Yes; -1:from Sill
  virtual bool         hasParam()     const { return false; }
  virtual String       getCovName()   const = 0;

protected:
  /// TODO : Gneiting (spatio-temporal covariance) :
  /// Change argument : double h becomes VectorDouble (number of sub-space)
  virtual double _evaluateCov (double h) const = 0;
  virtual double _evaluateCovDerivate(int degree, double h) const;

private:
  ECov        _type;    /*! Covariance function type */
  CovContext  _ctxt;    /*! Context (space, irfDegree, field, ...) */
  double      _param;   /*! Third parameter (TEST if not used) */
};

