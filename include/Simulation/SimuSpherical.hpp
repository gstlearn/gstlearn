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

namespace gstlrn
{

class MeshSpherical;

class GSTLEARN_EXPORT SimuSpherical: public ACalcSimulation
{
public:
  SimuSpherical(Id nbsimu = 1, Id seed = 4324324);
  SimuSpherical(const SimuSpherical &r) = delete;
  SimuSpherical& operator=(const SimuSpherical &r) = delete;
  virtual ~SimuSpherical();

  Id simulate(DbGrid *db,
               Model *model,
               const SimuSphericalParam& sphepar,
               Id iptr,
               bool verbose = false);

  VectorDouble simulate_mesh(MeshSpherical *mesh,
                             Model *model,
                             const SimuSphericalParam &sphepar,
                             bool verbose = false);

private:
  bool _run() override;

  static VectorDouble _spectrum_chentsov(const SimuSphericalParam& sphepar);
  static VectorDouble _spectrum_exponential(Model *model, const SimuSphericalParam& sphepar);
  static VectorDouble _spectrum_any(Model *model, const SimuSphericalParam& sphepar);
  static void _spectrum_normalize(Id verbose, VectorDouble& freqs);
  static Id _gdiscrete(VectorDouble& freqs);
  static Id _check_degree_order(const VectorDouble& freqs,
                                 VectorInt& degree,
                                 VectorInt& order,
                                 Id verbose);
};

GSTLEARN_EXPORT Id simsph(DbGrid *db,
                           Model *model,
                           const SimuSphericalParam& sphepar,
                           Id seed,
                           bool verbose,
                           const NamingConvention& namconv = NamingConvention("SimSphe"));
}
