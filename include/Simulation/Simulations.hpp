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

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "Simulation/SimuRefineParam.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class DbGrid;
  class ModelGeneric;
  class Model;
  class SimuFFTParam;
  class SimuPartitionParam;
  class SimuRefineParam;
  class SimuSubstitutionParam;

  GSTLEARN_EXPORT Id fluidPropagation(
    DbGrid* dbgrid,
    const String& name_facies,
    const String& name_fluid,
    const String& name_perm,
    const String& name_poro,
    Id nfacies,
    Id nfluids,
    Id niter = 1,
    const VectorInt& speeds = VectorInt(),
    bool show_fluid = false,
    double number_max = TEST,
    double volume_max = TEST,
    Id seed = 321321,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("Eden"));
  GSTLEARN_EXPORT Id simuFFT(
    DbGrid* db,
    ModelGeneric* model,
    SimuFFTParam& param,
    Id nbsimu = 1,
    Id seed = 432431,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("FFT"));
  GSTLEARN_EXPORT VectorDouble getChangeSupport(
    DbGrid* db,
    ModelGeneric* model,
    const SimuFFTParam& param,
    const VectorDouble& sigma = VectorDouble(),
    Id seed = 14333,
    bool verbose = false);
  GSTLEARN_EXPORT Id tessellationVoronoi(
    DbGrid* dbgrid,
    Model* model,
    const SimuPartitionParam& parparam,
    Id seed = 43243,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("Voronoi"));
  GSTLEARN_EXPORT Id tessellationPoisson(
    DbGrid* dbgrid,
    Model* model,
    const SimuPartitionParam& parparam,
    Id seed = 432432,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("Poisson"));
  GSTLEARN_EXPORT DbGrid* simuRefine(
    DbGrid* dbin,
    Model* model,
    const SimuRefineParam& param,
    Id seed = 432432,
    const NamingConvention& namconv = NamingConvention("Refine"));
  GSTLEARN_EXPORT Id simuSpectral(
    Db* dbin,
    Db* dbout,
    ModelGeneric* model,
    ANeigh* neigh = nullptr,
    Id nbsimu = 1,
    Id seed = 135672,
    Id ns = 10000,
    Id nd = 100,
    const ACov* cov0 = nullptr,
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("Simu"));
  GSTLEARN_EXPORT Id simuSubstitution(
    DbGrid* dbgrid,
    SimuSubstitutionParam& subparam,
    Id seed = 43242,
    Id verbose = false,
    const NamingConvention& namconv = NamingConvention("SimSub"));
  GSTLEARN_EXPORT Id simtub(
    Db* dbin = nullptr,
    Db* dbout = nullptr,
    Model* model = nullptr,
    ANeigh* neigh = nullptr,
    Id nbsimu = 1,
    Id seed = 43431,
    Id nbtuba = 100,
    bool flag_dgm = false,
    const NamingConvention& namconv = NamingConvention("Simu"));
  GSTLEARN_EXPORT Id simbayes(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    Id nbsimu = 1,
    Id seed = 132141,
    Id nbtuba = 100,
    const NamingConvention& namconv = NamingConvention("SimBayes"));
} // namespace gstlrn
