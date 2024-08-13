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

#include "Enum/EPowerPT.hpp"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/ALinearOp.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <map>

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#endif

class APolynomial;
class AMesh;
class Model;

class GSTLEARN_EXPORT PrecisionOp : public ALinearOp{

public:
  PrecisionOp();
  PrecisionOp(ShiftOpCs* shiftop,
              const CovAniso* cova,
              bool verbose = false);
  PrecisionOp(const AMesh* mesh,
              Model* model,
              int icov = 0,
              bool verbose = false);
  PrecisionOp(const PrecisionOp &pmat);
  PrecisionOp& operator=(const PrecisionOp &pmat);
  virtual ~PrecisionOp();

  // Interface functions for using PrecisionOp
  //virtual void evalDirect(const Eigen::VectorXd &vecin, Eigen::VectorXd &vecout);
  #ifndef SWIG
    virtual void evalSimulate(const Eigen::VectorXd& whitenoise, Eigen::VectorXd& vecout);
    virtual void evalInverse(const  Eigen::VectorXd& vecin, Eigen::VectorXd& vecout);
  #endif
  virtual void makeReady() {};

  virtual std::pair<double,double> getRangeEigenVal(int ndiscr = 100);

  static PrecisionOp* createFromShiftOp(ShiftOpCs *shiftop = nullptr,
                                        const CovAniso *cova = nullptr,
                                        bool verbose = false);
  static PrecisionOp* create(const AMesh* mesh,
                             Model* model,
                             int icov = 0,
                             bool verbose = false);

  int reset(const ShiftOpCs *shiftop,
            const CovAniso *cova = nullptr,
            bool verbose = false);

  virtual double getLogDeterminant(int nbsimu = 1, int seed = 0);
  #ifndef SWIG
    virtual void gradYQX(const Eigen::VectorXd& /*X*/,
                         const Eigen::VectorXd& /*Y*/,
                         Eigen::VectorXd& /*result*/,
                         const EPowerPT& /*power*/)
    {};
    virtual void gradYQXOptim(const Eigen::VectorXd& /*X*/,
                              const Eigen::VectorXd& /*Y*/,
                              Eigen::VectorXd& /*result*/,
                              const EPowerPT& /*power*/)
    {
    };
  virtual void evalDeriv(const Eigen::VectorXd& /*inv*/,
                         Eigen::VectorXd& /*outv*/,
                         int /*iapex*/,
                         int /*igparam*/,
                         const EPowerPT& /*power*/)
  {
  };
  virtual void evalDerivOptim(Eigen::VectorXd& /*outv*/,
                              int /*iapex*/,
                              int /*igparam*/,
                              const EPowerPT& /*power*/)
  {
  };
  std::vector<Eigen::VectorXd> simulate(int nbsimu = 1);

  #endif
  
//  virtual void evalDerivPoly(const Eigen::VectorXd& /*inv*/,
//                             Eigen::VectorXd& /*outv*/,
//                             int /*iapex*/,
//                             int /*igparam*/){};

  void evalPower(const Eigen::VectorXd &inv, Eigen::VectorXd &outv, const EPowerPT& power = EPowerPT::fromKey("ONE"));
  
  Eigen::VectorXd evalCov(int imesh);
  Eigen::VectorXd simulateOne();

  int  getSize() const { return _shiftOp->getSize(); }
  bool getTraining() const {return _training;}
  void setTraining(bool tr){ _training = tr;}
  ShiftOpCs* getShiftOp() const { return _shiftOp; }
  VectorDouble getPolyCoeffs(const EPowerPT& power);
  void setPolynomialFromPoly(APolynomial* polynomial);
  bool isCovaDefined() const { return _cova != nullptr; }
  VectorDouble getCoeffs();

protected:
  APolynomial*     getPoly(const EPowerPT& power);
  const ShiftOpCs* getShiftOpCs() const {return _shiftOp;}

#ifndef SWIG

public:
void evalPower(const VectorDouble &inv, VectorDouble &outv, const EPowerPT& power = EPowerPT::fromKey("ONE"));

protected:
  virtual void _addToDest(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const;
  void _addEvalPower(const Eigen::VectorXd& inv, Eigen::VectorXd& outv, const EPowerPT& power) const;

#endif

private:
  int  _preparePoly(const EPowerPT& power,bool force = false) const;
  int  _prepareChebychev(const EPowerPT& power) const;
  int  _preparePrecisionPoly() const;
#ifndef SWIG
  int  _evalPoly(const EPowerPT& power,const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
#endif
  void _purge();

private:
  mutable ShiftOpCs*                       _shiftOp;
  const CovAniso*                          _cova; // Not to be deleted
  mutable std::map<EPowerPT, APolynomial*> _polynomials;
  bool                                     _verbose;
  bool                                     _training;
  bool                                     _destroyShiftOp;
  bool                                     _userPoly;

#ifndef SWIG
protected :
  mutable Eigen::VectorXd              _work;
  mutable Eigen::VectorXd              _work2;
  mutable Eigen::VectorXd              _work3;
  mutable Eigen::VectorXd              _work4;
  mutable Eigen::VectorXd              _work5;
  mutable std::vector<Eigen::VectorXd> _workPoly;
#endif

};
