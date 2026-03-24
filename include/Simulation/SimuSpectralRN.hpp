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
#include "Simulation/CalcSimuSpectral.hpp"
#include "Simulation/SpectrumOnRN.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{

class ACov;
class ModelGeneric;

/**
 * Class for operating the Spectral simulations on RN
 */
class GSTLEARN_EXPORT SimuSpectralRN: public CalcSimuSpectral
{
public:
  SimuSpectralRN(Id nbsimu        = 1,
                 Id ns            = 10000,
                 Id nd            = 100,
                 Id seed          = 4324324,
                 const ACov* cov0 = nullptr,
                 bool verbose     = false);
  SimuSpectralRN(const SimuSpectralRN& r)            = delete;
  SimuSpectralRN& operator=(const SimuSpectralRN& r) = delete;
  virtual ~SimuSpectralRN();

  MatrixDense getOmega() { return _omega; };
  MatrixDense getGamma() { return _gamma; };
  SpectrumOnRN* getSpectrum() {return _sp;};

protected:
  bool _check() override;
  Id _simulate(const ACov* cova) override;
  Id _compute(Db* dbout, const VectorBool& activeArray, VectorVectorDouble& tab) override;

private:
  // spectrum for R^n (Not used...)
  MatrixDense _gamma; // Matrix nrows=_ns, ncols= number of variables of _cova
  MatrixDense _omega; // Matrix nrows=_ns, ncols= number of dimensions of _cova
 
  // spectrum on R^n
  SpectrumOnRN* _sp;



  const ACov* _cov0;
};

} // namespace gstlrn