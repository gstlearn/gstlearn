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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class ACov;
class ModelGeneric;
/**
 * Abstract Class for operating the Spectral simulations
 */
class GSTLEARN_EXPORT ASimuSpectral
{
public:
  ASimuSpectral(const ACov* cova = nullptr);
  ASimuSpectral(const ASimuSpectral& r);
  ASimuSpectral& operator=(const ASimuSpectral& r);
  virtual ~ASimuSpectral();
  void setCov(const ACov*& cova) { _cova = cova; };
  const ACov* getCov() { return _cova; };
  bool isPrepared() const { return _isPrepared; };
  VectorDouble getPhi() { return _phi; };
  Id getNs() const { return _phi.length(); };

  Id simulate(Id ns,
              Id seed          = 4273,
              bool verbose     = false,
              const ACov* cov0 = nullptr,
              Id nd            = 100);
  Id compute(Db* dbout,
             Id iuid                         = 0,
             bool verbose                    = false,
             const NamingConvention& namconv = NamingConvention("Simu"),
             const String& qualifier         = "simu");

protected:
  virtual Id _simulate(Id ns,
                       Id nd            = 100,
                       const ACov* cov0 = nullptr,
                       bool verbose     = false) = 0;

  virtual Id _compute(Db* dbout,
                      Id iuid      = 0,
                      bool verbose = false) = 0;

protected:
  bool _isPrepared;
  VectorDouble _phi; // Vector length=_ns
  const ACov* _cova; // Storing the poIder (not to be deleted)
};

GSTLEARN_EXPORT Id simuSpectral(Db* dbin,
                                Db* dbout,
                                const ACov* cova,
                                Id nbsimu                       = 1,
                                Id seed                         = 0,
                                Id ns                           = 10000,
                                Id nd                           = 100,
                                const ACov* cov0                = nullptr,
                                bool verbose                    = false,
                                const NamingConvention& namconv = NamingConvention(""));

} // namespace gstlrn