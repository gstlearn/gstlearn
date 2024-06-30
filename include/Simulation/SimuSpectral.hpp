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

#include "gstlearn_export.hpp"

#include "ACalcSimulation.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

class MatrixRectangular;

typedef struct
{
  int _k;
  int _countP;
  int _countM;
  std::map<int, std::map<int, int>> _tab;
} spSim;

/**
 * Class for operating the Spectral simulations
 */
class GSTLEARN_EXPORT SimuSpectral
{
public:
  SimuSpectral(const Model *model = nullptr);
  SimuSpectral(const SimuSpectral &r);
  SimuSpectral& operator=(const SimuSpectral &r);
  virtual ~SimuSpectral();

  int simulate(int nb, int seed = 4273);
  int simulateOnSphere(int ns, int nd, int seed = 4273, bool verbose = false);
  int compute(Db *dbout,
              const VectorDouble &xref = VectorDouble(),
              bool verbose = false,
              const NamingConvention& namconv = NamingConvention("Simu"));
  int computeOnSphere(Db *dbout,
                      bool verbose = false,
                      const NamingConvention& namconv = NamingConvention("Simu"));
  void setModel(const Model *&model) { _model = model; }

  static bool isValidForSpectral(const Model *model);

private:
  void _simulateOnSphere(bool verbose = false);

  void _printSpSim(const spSim& spsim, int status = 0) const;
  void _printSpSims(int status = 0);
  int _getKey1Maximum(const spSim& spsim) const;
  int _getSumValue(const spSim& spsim) const;
  VectorInt _getKeys1(const spSim& spsim) const;
  VectorInt _getKeys2(const spSim& spsim, int key1) const;
  VectorInt _getValues2(const spSim& spsim, int key1) const;

private:
  int _ndim;    // Space dimension
  int _nb;      // Number of spectral components
  int _nd;      // Number of degrees considered in the spectrum
  int _ns;      // Number of simulated harmonic components
  bool _isPrepared;
  VectorDouble _phi;
  VectorDouble _gamma;
  MatrixRectangular _omega; // Matrix nrows=nb, ncols=ndim
  std::vector<spSim> _spSims;

  const Model* _model; // Storing the pointer (not to be deleted)
};
