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
#include "Simulation/SpectrumRN.hpp"
#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixDense.hpp"

#include <cmath>

namespace gstlrn
{
  /**
   * ---------------------------------
   * Spectrum for simulation on Rn
   * ---------------------------------
   */
  SpectrumRN::SpectrumRN(const MatrixDense& gamma, const MatrixDense& omega)
    : _gamma(gamma)
    , _omega(omega)
    , _omega0(0, 0)
    , _xi(0)
  {
    if (_gamma.getNRows() != _omega.getNRows())
    {
      messerr("The number of simulated harmonic components is not consistent");
      _gamma.reset(0, 0);
      _omega.reset(0, 0);
      _omega0.reset(0, 0);
      _xi.resize(0);
    }
  }

  SpectrumRN::SpectrumRN(const MatrixDense& gamma,
                         const MatrixDense& omega,
                         MatrixDense& omega0,
                         VectorDouble& xi)
    : _gamma(gamma)
    , _omega(omega)
    , _omega0(omega0)
    , _xi(xi)
  {
    if (_gamma.getNRows() != _omega.getNRows())
    {
      messerr("The number of simulated harmonic components is not consistent");
      _gamma.reset(0, 0);
      _omega.reset(0, 0);
      _omega0.reset(0, 0);
      _xi.resize(0);
    }
  }

  SpectrumRN::SpectrumRN(const SpectrumRN& r)
    : _gamma(r._gamma)
    , _omega(r._omega)
    , _omega0(r._omega0)
    , _xi(r._xi)
  {
  }

  SpectrumRN& SpectrumRN::operator=(const SpectrumRN& r)
  {
    if (this != &r)
    {
      _gamma = r._gamma;
      _omega = r._omega;
      _omega0 = r._omega0;
      _xi = r._xi;
    }
    return *this;
  }

  SpectrumRN::~SpectrumRN() {}

} // namespace gstlrn
