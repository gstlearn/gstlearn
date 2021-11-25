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
#include "Basic/Vector.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/EPowerPT.hpp"

#include <map>

class APolynomial;

class GSTLEARN_EXPORT PrecisionOp {

public:
  PrecisionOp(ShiftOpCs* shiftop = nullptr,
              const CovAniso* cova = nullptr,
              const EPowerPT& power = EPowerPT::UNDEFINED,
              bool verbose = false);
  PrecisionOp(const PrecisionOp &pmat);
  PrecisionOp& operator=(const PrecisionOp &pmat);
  virtual ~PrecisionOp();

  int init(const ShiftOpCs* shiftop,
           const CovAniso* cova = nullptr,
           const EPowerPT& power = EPowerPT::UNDEFINED,
           bool verbose = false);

  void   eval(const VectorDouble& in, VectorDouble& out);
  virtual void gradYQX(const VectorDouble& /*X*/,
                       const VectorDouble& /*Y*/,
                       VectorDouble& /*result*/)
  {
  };

  virtual void gradYQXOptim(const VectorDouble& /*X*/,
                            const VectorDouble& /*Y*/,
                            VectorDouble& /*result*/)
  {
  };

  virtual void evalDeriv(const VectorDouble& /*in*/,
                         VectorDouble& /*out*/,
                         int /*iapex*/,
                         int /*igparam*/)
  {
  };

  virtual void evalDerivOptim(VectorDouble& /*out*/,
                              int /*iapex*/,
                              int /*igparam*/)
  {
  };

//  virtual void evalDerivPoly(const VectorDouble& /*in*/,
//                             VectorDouble& /*out*/,
//                             int /*iapex*/,
//                             int /*igparam*/){};

  int    getSize() const { return _shiftOp->getSize(); }
  double computeLogDet(int nsimus = 1, int seed = 0);
  bool getTraining()const {return _training;}
  void setTraining(bool tr){ _training = tr;}
  ShiftOpCs* getShiftOp() const { return _shiftOp; }

protected:
  APolynomial* getPoly(const EPowerPT& power);
  const EPowerPT&   getPower()const{return _power;}
  const ShiftOpCs* getShiftOpCs() const {return _shiftOp;}

private:
  int  _preparePoly(const EPowerPT& power);
  int  _prepareChebychev(const EPowerPT& power);
  int  _preparePrecisionPoly();

  int  _evalPoly(const EPowerPT& power,const VectorDouble& in, VectorDouble& out);
  void _purge();

private:
  mutable ShiftOpCs*               _shiftOp;
  const CovAniso*                  _cova;
  EPowerPT                         _power;
  std::map<EPowerPT, APolynomial*> _polynomials;
  bool                             _verbose;
  bool                             _training;



protected :
  mutable VectorDouble _work;
  mutable VectorDouble _work2;
  mutable VectorDouble _work3;
  mutable VectorDouble _work4;
  mutable VectorDouble _work5;
  mutable VectorVectorDouble _workPoly;
};
