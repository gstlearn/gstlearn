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

#include "Basic/Vector.hpp"
#include <vector>
class ALinearOpMulti;

class ALinearOpMulti {

public:
  ALinearOpMulti();
  virtual ~ALinearOpMulti();
  void evalDirect(const VectorVectorDouble& in,
                  VectorVectorDouble& out) const;
  virtual void evalInverse(const VectorVectorDouble& in,
                           VectorVectorDouble& out) const;
  virtual int size() const = 0;
  virtual int size(int) const = 0;

  void setNIterMax(int nitermax) { _nIterMax = nitermax; }
  void setEps(double eps) { _eps = eps; }
  void setPrecond(const ALinearOpMulti* precond, int status);

  /*! Reset the Conjugate Gradient statistics */
  void resetStatCG() const;

  /*! Print out the Conjugate Gradient statistics */
  void printStatCG() const;

protected:
  void _init() const;
  virtual void _evalDirect(const VectorVectorDouble& in,
                           VectorVectorDouble& out) const = 0;
  void _linearComb(double val1,
                   const VectorVectorDouble& in1,
                   double val2,
                   const VectorVectorDouble& in2,
                   VectorVectorDouble& out) const;
  void _copyVals(const VectorVectorDouble& in,
                 VectorVectorDouble& out) const;
  void _updated() const;
  double _prod(const VectorDouble& x,const VectorDouble& y) const;
  double _prod(const VectorVectorDouble& x,
               const VectorVectorDouble& y) const;

private:
  void _diff(const VectorVectorDouble&,
             const VectorVectorDouble&,
             VectorVectorDouble&) const;

  void _fillVal(VectorVectorDouble& vect,double val)const;

private:
  int                       _nIterMax;
  double                    _eps;
  bool                      _precondStatus;
  bool                      _userInitialValue;
  const ALinearOpMulti*     _precond;


  // Work arrays

  mutable bool               _initialized;
  mutable VectorVectorDouble _z;
  mutable VectorVectorDouble _r;
  mutable VectorVectorDouble _temp;
  mutable VectorVectorDouble _p;

  // Environment parameters

  mutable double     _timeCG;
  mutable int        _niterCG;
  mutable int        _numberCG;
};
