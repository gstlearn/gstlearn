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

#include "geoslib_define.h"

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuSphericalParam.hpp"

class MeshSpherical;

class GSTLEARN_EXPORT SimuSpherical: public ACalcSimulation
{
public:
  SimuSpherical(int nbsimu = 1, int seed = 4324324);
  SimuSpherical(const SimuSpherical &r) = delete;
  SimuSpherical& operator=(const SimuSpherical &r) = delete;
  virtual ~SimuSpherical();

  int simulate(DbGrid *db,
               Model *model,
               const SimuSphericalParam& sphepar,
               int iptr,
               bool verbose = false);

  VectorDouble simulate_mesh(MeshSpherical *mesh,
                             Model *model,
                             const SimuSphericalParam &sphepar,
                             bool verbose = false);

private:
  virtual bool _run() override;

  static VectorDouble _spectrum_chentsov(const SimuSphericalParam& sphepar);
  static VectorDouble _spectrum_exponential(Model *model, const SimuSphericalParam& sphepar);
  static VectorDouble _spectrum_any(Model *model, const SimuSphericalParam& sphepar);
  static void _spectrum_normalize(int verbose, VectorDouble& freqs);
  static int _gdiscrete(VectorDouble& freqs);
  static int _check_degree_order(const VectorDouble& freqs,
                                 VectorInt& degree,
                                 VectorInt& order,
                                 int verbose);
};

GSTLEARN_EXPORT int simsph(DbGrid *db,
                           Model *model,
                           const SimuSphericalParam& sphepar,
                           int seed,
                           bool verbose,
                           const NamingConvention& namconv = NamingConvention("SimSphe"));
