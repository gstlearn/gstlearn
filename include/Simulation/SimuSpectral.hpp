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
  int _k;
  int _countP;
  int _countM;
  std::map<int, std::map<int, int>> _tab;
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

  int simulate(int ns, int seed = 4273, bool verbose = false, int nd = 100);
  int compute(Db *dbout,
              int iuid = 0,
              bool verbose = false,
              const NamingConvention& namconv = NamingConvention("Simu"));

  static bool isValidForSpectral(const ACov* cova);

  void setCov(const ACov *&cova) { _cova = cova; }
  void setNdim(int ndim) { _ndim = ndim; }
  void setNs(int ns) { _ns = ns; }

private:
  void _simulateOnSphere(int nd = 100, bool verbose = false);
  void _simulateOnRn();
  void _computeOnSphere(Db* dbout, int iuid, bool verbose = false);
  void _computeOnRn(Db *dbout, int iuid, bool verbose = false);

  static void _printSpSim(const spSim& spsim, int status = 0);
  void _printSpSims(int status = 0);
  static int _getKey1Maximum(const spSim& spsim);
  static int _getSumValue(const spSim& spsim);
  static VectorInt _getKeys1(const spSim& spsim);
  static VectorInt _getKeys2(const spSim& spsim, int key1);
  static VectorInt _getValues2(const spSim& spsim, int key1);

private:
  int _ndim;    // Space dimension
  int _ns;      // Number of simulated harmonic components
  bool _isPrepared;
  VectorDouble _phi;
  VectorDouble _gamma;
  MatrixDense _omega; // Matrix nrows=_ns, ncols=ndim
  std::vector<spSim> _spSims;

  const ACov* _cova; // Storing the pointer (not to be deleted)
};

GSTLEARN_EXPORT int simuSpectral(Db *dbin = nullptr,
                                 Db *dbout = nullptr,
                                 ModelGeneric* model = nullptr,
                                 int nbsimu = 1,
                                 int seed = 43431,
                                 int ns = 100,
                                 int nd = 100,
                                 bool verbose = false,
                                 const NamingConvention &namconv = NamingConvention("Simu"));
}