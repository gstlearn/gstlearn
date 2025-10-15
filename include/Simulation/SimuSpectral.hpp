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

#include "Matrix/MatrixDense.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
typedef struct
{
  Id _k;
  Id _countP;
  Id _countM;
  std::map<Id, std::map<Id, Id>> _tab;
} spSim;

class ACov;
class ModelGeneric;
/**
 * Class for operating the Spectral simulations
 */
class GSTLEARN_EXPORT SimuSpectral
{
public:
  SimuSpectral(const ACov *cova = nullptr);
  SimuSpectral(const SimuSpectral &r);
  SimuSpectral& operator=(const SimuSpectral &r);
  virtual ~SimuSpectral();

  Id simulate(Id ns, Id seed = 4273, bool verbose = false, Id nd = 100);
  Id compute(Db *dbout,
              Id iuid = 0,
              bool verbose = false,
              const NamingConvention& namconv = NamingConvention("Simu"));

  static bool isValidForSpectral(const ACov* cova);

  void setCov(const ACov *&cova) { _cova = cova; }
  void setNdim(Id ndim) { _ndim = ndim; }
  void setNs(Id ns) { _ns = ns; }

private:
  void _simulateOnSphere(Id nd = 100, bool verbose = false);
  void _simulateOnRn();
  void _computeOnSphere(Db* dbout, Id iuid, bool verbose = false);
  void _computeOnRn(Db *dbout, Id iuid, bool verbose = false);

  static void _printSpSim(const spSim& spsim, Id status = 0);
  void _printSpSims(Id status = 0);
  static Id _getKey1Maximum(const spSim& spsim);
  static Id _getSumValue(const spSim& spsim);
  static VectorInt _getKeys1(const spSim& spsim);
  static VectorInt _getKeys2(const spSim& spsim, Id key1);
  static VectorInt _getValues2(const spSim& spsim, Id key1);

private:
  Id _ndim;    // Space dimension
  Id _ns;      // Number of simulated harmonic components
  bool _isPrepared;
  VectorDouble _phi;
  VectorDouble _gamma;
  MatrixDense _omega; // Matrix nrows=_ns, ncols=ndim
  std::vector<spSim> _spSims;

  const ACov* _cova; // Storing the pointer (not to be deleted)
};

GSTLEARN_EXPORT Id simuSpectral(Db *dbin = nullptr,
                                 Db *dbout = nullptr,
                                 ModelGeneric* model = nullptr,
                                 Id nbsimu = 1,
                                 Id seed = 43431,
                                 Id ns = 100,
                                 Id nd = 100,
                                 bool verbose = false,
                                 const NamingConvention &namconv = NamingConvention("Simu"));
}