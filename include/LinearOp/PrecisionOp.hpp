/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EPowerPT.hpp"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ShiftOpCs.hpp"

#include <map>

class APolynomial;
class AMesh;
class Model;

class GSTLEARN_EXPORT PrecisionOp {

public:
  PrecisionOp();
  PrecisionOp(ShiftOpCs* shiftop,
              const CovAniso* cova,
              bool verbose = false);
  PrecisionOp(const AMesh* mesh,
              Model* model,
              int icov = 0,
              bool verbose = false);
  PrecisionOp(const PrecisionOp &m);
  PrecisionOp& operator=(const PrecisionOp &m);
  virtual ~PrecisionOp();

  virtual std::pair<double,double> getRangeEigenVal(int ndiscr = 100);

  static PrecisionOp* createFromShiftOp(ShiftOpCs *shiftop = nullptr,
                                        const CovAniso *cova = nullptr,
                                        bool verbose = false);
  static PrecisionOp* create(const AMesh *mesh,
                             Model *model,
                             int icov = 0,
                             bool verbose = false);

  int reset(const ShiftOpCs* shiftop,
            const CovAniso* cova = nullptr,
            bool verbose = false);

  // Interface functions for using PrecisionOp
  virtual void eval(const VectorDouble &inv, VectorDouble &outv);
  virtual void simulateOneInPlace(VectorDouble& whitenoise, VectorDouble& result);
  virtual void evalInvVect(VectorDouble& in, VectorDouble& result);
  virtual double computeLogDet(int nsimus = 1, int seed = 0);
  virtual void gradYQX(const VectorDouble& /*X*/,
                       const VectorDouble& /*Y*/,
                       VectorDouble& /*result*/,
                       const EPowerPT& /*power*/)
  {
  };
  virtual void gradYQXOptim(const VectorDouble& /*X*/,
                            const VectorDouble& /*Y*/,
                            VectorDouble& /*result*/,
                            const EPowerPT& /*power*/)
  {
  };
  virtual void evalDeriv(const VectorDouble& /*inv*/,
                         VectorDouble& /*outv*/,
                         int /*iapex*/,
                         int /*igparam*/,
                         const EPowerPT& /*power*/)
  {
  };
  virtual void evalDerivOptim(VectorDouble& /*outv*/,
                              int /*iapex*/,
                              int /*igparam*/,
                              const EPowerPT& /*power*/)
  {
  };
//  virtual void evalDerivPoly(const VectorDouble& /*inv*/,
//                             VectorDouble& /*outv*/,
//                             int /*iapex*/,
//                             int /*igparam*/){};

  void evalPower(const VectorDouble &inv, VectorDouble &outv, const EPowerPT& power = EPowerPT::fromKey("ONE"));
  VectorDouble evalCov(int imesh);
  VectorVectorDouble simulate(int nbsimus = 1);
  VectorDouble simulateOne();

  int    getSize() const { return _shiftOp->getSize(); }
  bool getTraining()const {return _training;}
  void setTraining(bool tr){ _training = tr;}
  ShiftOpCs* getShiftOp() const { return _shiftOp; }
  VectorDouble getPolyCoeffs(EPowerPT power);
  void setPolynomialFromPoly(APolynomial* polynomial);
  bool isCovaDefined() const { return _cova != nullptr; }
  VectorDouble getCoeffs();

  // Talking to ShiftOp
  void setNIterMax(int nitermax);
  void setEps(double eps);

protected:
  APolynomial*     getPoly(const EPowerPT& power);
  const ShiftOpCs* getShiftOpCs() const {return _shiftOp;}

private:
  int  _preparePoly(const EPowerPT& power,bool force = false);
  int  _prepareChebychev(const EPowerPT& power);
  int  _preparePrecisionPoly();
  int  _evalPoly(const EPowerPT& power,const VectorDouble& inv, VectorDouble& outv);
  void _purge();

private:
  mutable ShiftOpCs*               _shiftOp;
  const CovAniso*                  _cova;
  std::map<EPowerPT, APolynomial*> _polynomials;
  bool                             _verbose;
  bool                             _training;
  bool                             _destroyShiftOp;
  bool                             _userPoly;

protected :
  mutable VectorDouble _work;
  mutable VectorDouble _work2;
  mutable VectorDouble _work3;
  mutable VectorDouble _work4;
  mutable VectorDouble _work5;
  mutable VectorVectorDouble _workPoly;
};
