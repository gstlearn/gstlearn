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

#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Enum/EPowerPT.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "LinearOp/IPrecisionOp.hpp"
#include <map>
#include <memory>

namespace gstlrn
{
  class APolynomial;
  class AMesh;

  // This class create a precision operator (matrix-free).
  // In general, it is built from a Model and a AMesh
  // Note that if the model is multivariate, the precision is built with a constant sill = 1.
  // Therefore it has to be used only through the PrecisionOpMulti class
  // which handles the sills matrix (possibly non stationary)
  class GSTLEARN_EXPORT PrecisionOp: virtual public IPrecisionOp
  {
  public:
    PrecisionOp();
    PrecisionOp(AShiftOp* shiftop, const CovAniso* cova, bool verbose = false);
    PrecisionOp(
      const AMesh* mesh,
      const CovAniso* cova,
      bool stencil = false,
      bool verbose = false);
    PrecisionOp(const PrecisionOp& pmat);
    PrecisionOp(PrecisionOp&& pmat) noexcept;
    PrecisionOp& operator=(const PrecisionOp& pmat);
    PrecisionOp& operator=(PrecisionOp&& pmat) noexcept;
    virtual ~PrecisionOp();

    // Interface functions for using PrecisionOp

#ifndef SWIG
    virtual void evalInverse(const constvect vecin, VectorDouble& vecout);
#endif
    VectorDouble evalInverse(const VectorDouble& vecin);
    std::pair<double, double> rangeEigenVal(Id ndiscr = 100) const override;

    static PrecisionOp* createFromShiftOp(
      AShiftOp* shiftop = nullptr,
      const CovAniso* cova = nullptr,
      bool verbose = false);
    static PrecisionOp* create(
      const AMesh* mesh,
      CovAniso* cova,
      bool stencil = false,
      bool verbose = false);

    Id reset(
      const AShiftOp* shiftop,
      const CovAniso* cova = nullptr,
      bool verbose = false);

#ifndef SWIG
    virtual void gradYQX(
      const constvect /*X*/,
      const constvect /*Y*/,
      vect /*result*/,
      const EPowerPT& /*power*/) {};
    virtual void gradYQXOptim(
      const constvect /*X*/,
      const constvect /*Y*/,
      vect /*result*/,
      const EPowerPT& /*power*/) {};
    virtual void evalDeriv(
      const constvect /*inv*/,
      vect /*outv*/,
      Id /*iapex*/,
      Id /*igparam*/,
      const EPowerPT& /*power*/) {};
    virtual void evalDerivOptim(
      vect /*outv*/,
      Id /*iapex*/,
      Id /*igparam*/,
      const EPowerPT& /*power*/) {};
    VectorVectorDouble simulate(Id nbsimu = 1);

#endif

    //  virtual void evalDerivPoly(const Eigen::VectorXd& /*inv*/,
    //                             Eigen::VectorXd& /*outv*/,
    //                             Id /*iapex*/,
    //                             Id /*igparam*/){};

#ifndef SWIG
    void evalPower(
      const constvect inm,
      vect outm,
      const EPowerPT& power = EPowerPT::fromKey("ONE"));
    void evalPower(
      const VectorDouble& inv,
      VectorDouble& outv,
      const EPowerPT& power = EPowerPT::fromKey("ONE"));

#endif
    VectorDouble computeCov(Id imesh);
    VectorDouble simulateOne();

    Id getSize() const override { return _shiftOp->getSize(); }

    bool getTraining() const { return _training; }

    void setTraining(bool tr) { _training = tr; }

    AShiftOp* getShiftOp() const { return _shiftOp; }

    VectorDouble getPolyCoeffs(const EPowerPT& power);
    void setPolynomialFromPoly(APolynomial* polynomial);

    bool isCovaDefined() const { return _cova != nullptr; }

    VectorDouble getCoeffs();
    double computeLogDet(Id nMC = 1) const override;
    virtual VectorDouble extractDiag() const;

  protected:
    APolynomial* getPoly(const EPowerPT& power);

#ifndef SWIG
    Id
      _addEvalPoly(const EPowerPT& power, const constvect inv, vect outv) const;
    Id _addToDest(const constvect inv, vect outv) const override;
    Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;

    void _addEvalPower(const constvect inv, vect outv, const EPowerPT& power)
      const;
#endif

  private:
    Id _preparePoly(const EPowerPT& power, bool force = false) const;
    Id _prepareChebychev(const EPowerPT& power) const;
    Id _preparePrecisionPoly() const;
#ifndef SWIG
    Id _evalPoly(const EPowerPT& power, const constvect inv, vect outv) const;
#endif
    void _purge();

  private:
    mutable AShiftOp* _shiftOp;
    const CovAniso* _cova; // Not to be deleted
    mutable std::map<EPowerPT, std::unique_ptr<APolynomial>> _polynomials;
    bool _verbose;
    bool _training;
    bool _destroyShiftOp;
    bool _userPoly;

#ifndef SWIG

  protected:
    mutable VectorDouble _work;
    mutable VectorDouble _work2;
    mutable VectorVectorDouble _workPoly;
#endif
  };
} // namespace gstlrn
