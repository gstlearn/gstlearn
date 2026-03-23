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
#include "Matrix/MatrixDense.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{

  class ACov;
  class ModelGeneric;
  class ASimuSpectral;

  /**
   * Auxiliary class for the interface between ACov and SimuSpectral for operating the Spectral simulations on RN
   */
  class GSTLEARN_EXPORT SpectrumRN
  {
  public:
    SpectrumRN(const MatrixDense& gamma = MatrixDense(),
               const MatrixDense& omega = MatrixDense());
    SpectrumRN(const MatrixDense& gamma,
               const MatrixDense& omega,
               MatrixDense& omega0,
               VectorDouble& xi);
    SpectrumRN(const SpectrumRN& r);
    SpectrumRN& operator=(const SpectrumRN& r);
    virtual ~SpectrumRN();

    Id getNDim() const { return _omega.getNCols(); };

    Id getNVar() const { return _gamma.getNCols(); };

    Id getNs() const { return _omega.getNRows(); };

    MatrixDense getOmega() { return _omega; };

    MatrixDense getGamma() { return _gamma; };

    MatrixDense getOmega0() const { return _omega0; };

    VectorDouble getXi() const { return _xi; };

  private:
    // spectrum for R^n
    MatrixDense _gamma; // Matrix nrows=_ns, ncols= number of variables of _nvar
    MatrixDense
      _omega; // Matrix nrows=_ns, ncols= number of dimensions of _ndim
    // used for CorGneiting
    MatrixDense
      _omega0; // Matrix nrows=_ns, ncols= number of variables of _nvar (spatial frequencies befor rotation)
    VectorDouble _xi; // Vector of random scales, length =_ns
  };

} // namespace gstlrn
