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

/**
 * Class for operating the Spectral simulations
 */
class GSTLEARN_EXPORT  SimuSpectral
{
public:
  SimuSpectral(const Model *model = nullptr);
  SimuSpectral(const SimuSpectral &r);
  SimuSpectral& operator=(const SimuSpectral &r);
  virtual ~SimuSpectral();

  int simulate(int nb, int seed = 4273);
  int simulateOnSphere(int nb, int seed = 4273);
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
  VectorVectorInt _simulateOnSphereV0();
  VectorVectorInt _simulateOnSphere();

private:
  int _ndim;    // Space dimension
  int _nb;      // Number of spectral components
  bool _isPrepared;
  VectorDouble _phi;
  VectorDouble _gamma;
  MatrixRectangular _omega; // Matrix nrows=nb, ncols=ndim

  const Model* _model; // Storing the pointer (not to be deleted)
};
