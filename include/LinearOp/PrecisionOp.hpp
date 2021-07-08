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
#include "geoslib_enum.h"
#include <map>
#include "LinearOp/ShiftOpCs.hpp"

class APolynomial;
class PrecisionOp {

public:
  PrecisionOp(ShiftOpCs* shiftop = nullptr,
              const Model* model = nullptr,
              int icov = 0,
              ENUM_POPTS power = POPT_UNDEFINED,
              bool verbose = false);
  PrecisionOp(const PrecisionOp &pmat);
  PrecisionOp& operator=(const PrecisionOp &pmat);
  virtual ~PrecisionOp();

  int init(const ShiftOpCs* shiftop,
           const Model* model = nullptr,
           int icov = 0,
           ENUM_POPTS power = POPT_UNDEFINED,
           bool verbose = false);

  void   eval(const VectorDouble& in, VectorDouble& out);
  virtual void   evalDeriv(const VectorDouble& in, VectorDouble& out){};

  int    getSize() const { return _shiftOp->getSize(); }
  double computeLogDet(int nsimus = 1, int seed = 0);

  ShiftOpCs* getShiftOp() const { return _shiftOp; }

protected:
  APolynomial* getPoly(ENUM_POPTS power);
  ENUM_POPTS getPower()const{return _power;}
  const ShiftOpCs* getShiftOpCs() const {return _shiftOp;}

private:
  int  _preparePoly(ENUM_POPTS power);
  int  _prepareChebychev(ENUM_POPTS power);
  int  _preparePrecisionPoly();

  int  _evalPoly(ENUM_POPTS power,const VectorDouble& in, VectorDouble& out);
  void _purge();


private:
  mutable ShiftOpCs*     _shiftOp;
  const CovAniso*      _cova;
  ENUM_POPTS           _power;
  std::map<ENUM_POPTS, APolynomial*> _polynomials;
  bool                 _verbose;

protected :
  mutable VectorDouble _work;

};
