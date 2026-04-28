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

#include "Simulation/ACalcSimuGaussian.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Model;

  class GSTLEARN_EXPORT CalcSimuSpectral: public ACalcSimuGaussian
  {
  public:
    CalcSimuSpectral(
      Id nbsimu = 1,
      Id ns = 10000,
      Id nd = 100,
      Id seed = 4324324,
      bool verbose = false);
    CalcSimuSpectral(const CalcSimuSpectral& r) = delete;
    CalcSimuSpectral& operator=(const CalcSimuSpectral& r) = delete;
    virtual ~CalcSimuSpectral();

    // Perform the task for ONE simulation (debugging option)
    Id simulateSpectralTest();
    Id computeSpectralTest(Db* dbout, Id isimu = 0);

  protected:
    bool _check() override;
    bool _preprocess() override;

    // bool _postprocess() override;

    Id _getNs() const { return _ns; };

    Id _getNd() const { return _nd; };

  private:
    Id _ns; // Number of spectral components
    Id _nd; // Maximum number of spectral orders on
  };

} // namespace gstlrn
