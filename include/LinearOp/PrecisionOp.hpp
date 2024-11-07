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

#include "Covariances/CovAniso.hpp"
#include "LinearOp/ASimulable.hpp"

#include "Enum/EPowerPT.hpp"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include <map>

class APolynomial;
class AMesh;

// This class create a precision operator (matrix-free).
// In general, it is built from a Model and a AMesh
// Note that if the model is multivariate, the precision is built with a constant sill = 1.
// Therefore it has to be used only through the PrecisionOpMulti class
// which handles the sills matrix (possibly non stationary)
class GSTLEARN_EXPORT PrecisionOp : public ASimulable {

public:
  PrecisionOp();
  PrecisionOp(ShiftOpCs* shiftop,
              const CovAniso* cova,
              bool verbose = false);
  PrecisionOp(const AMesh* mesh,
              CovAniso* cova,
              bool verbose = false);
  PrecisionOp(const PrecisionOp &pmat);
  PrecisionOp& operator=(const PrecisionOp &pmat);
  virtual ~PrecisionOp();

  // Interface functions for using PrecisionOp

  #ifndef SWIG
  virtual void evalInverse(const constvect vecin, std::vector<double>& vecout);
#endif

  virtual std::pair<double,double> getRangeEigenVal(int ndiscr = 100);

  static PrecisionOp* createFromShiftOp(ShiftOpCs *shiftop = nullptr,
                                        const CovAniso *cova = nullptr,
                                        bool verbose = false);
  static PrecisionOp* create(const AMesh* mesh,
                             CovAniso* cova,
                             bool verbose = false);

  int reset(const ShiftOpCs *shiftop,
            const CovAniso *cova = nullptr,
            bool verbose = false);

  virtual double getLogDeterminant(int nbsimu = 1);
  #ifndef SWIG
  virtual void gradYQX(const constvect /*X*/,
                       const constvect /*Y*/,
                       vect /*result*/,
                       const EPowerPT& /*power*/) {};
  virtual void gradYQXOptim(const constvect /*X*/,
                            const constvect /*Y*/,
                            vect /*result*/,
                            const EPowerPT& /*power*/) {};
  virtual void evalDeriv(const constvect /*inv*/,
                         vect /*outv*/,
                         int /*iapex*/,
                         int /*igparam*/,
                         const EPowerPT& /*power*/) {};
  virtual void evalDerivOptim(vect /*outv*/,
                              int /*iapex*/,
                              int /*igparam*/,
                              const EPowerPT& /*power*/) {};
  VectorVectorDouble simulate(int nbsimu = 1);

  #endif
  
//  virtual void evalDerivPoly(const Eigen::VectorXd& /*inv*/,
//                             Eigen::VectorXd& /*outv*/,
//                             int /*iapex*/,
//                             int /*igparam*/){};

  #ifndef SWIG
  void evalPower(const constvect inm,
                 vect outm,
                 const EPowerPT& power = EPowerPT::fromKey("ONE"));
#endif
  VectorDouble evalCov(int imesh);
  VectorDouble simulateOne();

  int  getSize() const override { return _shiftOp->getSize(); }
  bool getTraining() const {return _training;}
  void setTraining(bool tr){ _training = tr;}
  ShiftOpCs* getShiftOp() const { return _shiftOp; }
  VectorDouble getPolyCoeffs(const EPowerPT& power);
  void setPolynomialFromPoly(APolynomial* polynomial);
  bool isCovaDefined() const { return _cova != nullptr; }
  VectorDouble getCoeffs();

  virtual VectorDouble extractDiag() const;

protected:
  APolynomial*     getPoly(const EPowerPT& power);
  const ShiftOpCs* getShiftOpCs() const {return _shiftOp;}
  const CovAniso*  getCova() const {return _cova;}

#ifndef SWIG

public:
  void evalPower(const VectorDouble &inv, VectorDouble &outv, const EPowerPT& power = EPowerPT::fromKey("ONE"));

protected:
  virtual int _addToDest(const constvect inv, vect outv) const override;
  virtual int _addSimulateToDest(const constvect whitenoise,
                                 vect outv) const override;

  void _addEvalPower(const constvect inv, vect outv, const EPowerPT& power) const;
#endif

private:
  int  _preparePoly(const EPowerPT& power,bool force = false) const;
  int  _prepareChebychev(const EPowerPT& power) const;
  int  _preparePrecisionPoly() const;
#ifndef SWIG
  int  _evalPoly(const EPowerPT& power, const constvect inv, vect outv) const;
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
  mutable std::vector<double>              _work;
  mutable std::vector<double>              _work2;
  mutable std::vector<double>              _work3;
  mutable std::vector<double>              _work4;
  mutable std::vector<double>              _work5;
  mutable std::vector<std::vector<double>> _workPoly;
#endif

};
