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

#include <vector>

class ALinearOpMulti;

class GSTLEARN_EXPORT ALinearOpMulti {

public:
  ALinearOpMulti(int nitermax = 1000, double eps = EPSILON8);
  ALinearOpMulti(const ALinearOpMulti &m);
  ALinearOpMulti& operator=(const ALinearOpMulti &m);
  virtual ~ALinearOpMulti();

  virtual void evalInverse(const VectorVectorDouble &inv,
                           VectorVectorDouble &outv) const;

  void evalDirect(const VectorVectorDouble& inv,
                  VectorVectorDouble& outv) const;
  void initLk(const VectorVectorDouble& inv,
                           VectorVectorDouble& outv) const;
  virtual int sizes() const = 0;
  virtual int size(int) const = 0;

  void setNIterMax(int nitermax) { _nIterMax = nitermax; }
  void setEps(double eps) { _eps = eps; }
  void setPrecond(const ALinearOpMulti* precond, int status);

  /*! Reset the Conjugate Gradient statistics */
  void resetStatCG() const;

  /*! Print out the Conjugate Gradient statistics */
  void printStatCG() const;


  void _linearComb(double val1,
                   const VectorVectorDouble& in1,
                   double val2,
                   const VectorVectorDouble& in2,
                   VectorVectorDouble& outv) const;
  void prodScalar(double val,
                  const VectorVectorDouble& inv,
                  VectorVectorDouble& outv) const;
  void  addProdScalar(double val,
                      const VectorVectorDouble& inv,
                      VectorVectorDouble& outv) const;
  void _copyVals(const VectorVectorDouble& inv,
                 VectorVectorDouble& outv) const;
  void _updated() const;
  double innerProduct(const VectorDouble& x,const VectorDouble& y) const;
  double innerProduct(const VectorVectorDouble& x,
               const VectorVectorDouble& y) const;
  double max(const VectorVectorDouble& vect) const;
  void  fillVal(VectorVectorDouble& vect,double val)const;
  void diff(const VectorVectorDouble&,
              const VectorVectorDouble&,
              VectorVectorDouble&) const;

  void sum(const VectorVectorDouble&,
                const VectorVectorDouble&,
                VectorVectorDouble&) const;
  mutable VectorVectorDouble _temp;
  mutable VectorVectorDouble _p;
  mutable VectorVectorDouble _z;

  void _initPublic() const;

protected:
  void _init() const;
  virtual void _evalDirect(const VectorVectorDouble& inv,
                           VectorVectorDouble& outv) const = 0;

private:
  int                       _nIterMax;
  double                    _eps;
  bool                      _precondStatus;
  bool                      _userInitialValue;
  const ALinearOpMulti*     _precond;

  // Work arrays
  mutable bool               _initialized;
  mutable VectorVectorDouble _r;

  // Environment parameters
  mutable double     _timeCG;
  mutable int        _niterCG;
  mutable int        _numberCG;
};
