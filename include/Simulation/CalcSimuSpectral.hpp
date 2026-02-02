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

#include "Simulation/ACalcSimulation.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Model;

class GSTLEARN_EXPORT CalcSimuSpectral: public ACalcSimulation
{
public:
  CalcSimuSpectral(Id nbsimu = 1, Id ns = 10000, Id nd = 100, Id seed = 4324324, bool verbose = false);
  CalcSimuSpectral(const CalcSimuSpectral& r)            = delete;
  CalcSimuSpectral& operator=(const CalcSimuSpectral& r) = delete;
  virtual ~CalcSimuSpectral();

  // Perform the task for ONE simulation
  Id simulate();
  Id compute(Db* dbout, Id isimu = 0);

  VectorDouble getPhi() { return _phi; };
  double getPhi(Id i) { return _phi[i]; };
  static bool isValidForSpectral(const ModelGeneric* model);
  bool getVerbose() const { return _verbose; }

protected:
  bool _check() override;
  bool _preprocess() override;
  bool _postprocess() override;
  bool _run() override;

  virtual Id _simulate()                       = 0;
  virtual Id _compute(Db* dbout, Id isimu = 0) = 0;

  Id _getNs() const { return _ns; };
  Id _getNd() const { return _nd; };
  Id _getNDim() const;
  Id _getNVar() const;

  void _setIsPrepared(bool status) { _isPrepared = status; }
  bool _getIsPrepared() const { return _isPrepared; }

private:
  bool _isPrepared;
  bool _verbose;
  Id _iattOut;
  Id _ns;            // Number of spectral components
  Id _nd;            // Maximum number of spectral orders on
  VectorDouble _phi; // Vector length=_ns
};

GSTLEARN_EXPORT Id simuSpectral(Db* dbin,
                                Db* dbout,
                                ModelGeneric* model,
                                ANeigh* neigh                   = nullptr,
                                Id nbsimu                       = 1,
                                Id seed                         = 135672,
                                Id ns                           = 10000,
                                Id nd                           = 100,
                                const ACov* cov0                = nullptr,
                                bool verbose                    = false,
                                const NamingConvention& namconv = NamingConvention("Simu"));

} // namespace gstlrn