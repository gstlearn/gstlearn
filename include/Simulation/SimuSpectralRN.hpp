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

#include "Basic/NamingConvention.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/ASimuSpectral.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{

class ACov;
class ModelGeneric;
class ASimuSpectral;

/**
 * Class for operating the Spectral simulations on RN
 */
class GSTLEARN_EXPORT SimuSpectralRN: public ASimuSpectral
{
public:
  SimuSpectralRN(const ACov* cova = nullptr, Id ns = 10000, Id nd = 100);
  SimuSpectralRN(const SimuSpectralRN& r);
  SimuSpectralRN& operator=(const SimuSpectralRN& r);
  virtual ~SimuSpectralRN();

  MatrixDense getOmega() { return _omega; };
  MatrixDense getGamma() { return _gamma; };

protected:
  Id _simulate(const ACov* cov0 = nullptr,
               bool verbose     = false) override;
  Id _compute(Db* dbout,
              Id iuid      = 0,
              bool verbose = false) override;

private:
  // spectrum for R^n
  MatrixDense _gamma; // Matrix nrows=_ns, ncols= number of variables of _cova
  MatrixDense _omega; // Matrix nrows=_ns, ncols= number of dimensions of _cova
};

} // namespace gstlrn