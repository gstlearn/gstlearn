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

#include "Enum/ECalcVario.hpp"

#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
class ECalcVario;
class Db;

/**
 * \brief
 * Class containing the Variogram Map which uses an DbGrid provided by the user
 * This function simply calculate and add the results as new field in this DbGrid.
 *
 */
class GSTLEARN_EXPORT VMap : public AVario
{
public:
  VMap(DbGrid* dbmap);
  VMap(const VMap& r);
  VMap& operator=(const VMap& r);
  virtual ~VMap();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(VMap)

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVario Interface
  double _getIVAR(const Db *db, Id iech, Id ivar) const override;
  void _setResult(Id iech1,
                  Id iech2,
                  Id nvar,
                  Id ilag,
                  Id ivar,
                  Id jvar,
                  Id orient,
                  double ww,
                  double dist,
                  double value) override;

  Id compute(Db *db,
              const ECalcVario &calcul_type,
              Id radius,
              bool flag_FFT = true,
              const NamingConvention &namconv = NamingConvention("VMAP"));

private:
  Id _grid_fft(DbGrid *dbgrid,
                const NamingConvention &namconv);
  static void _extract(const Id* nxmap,
                       const Id* nxgrid,
                       Id* dims,
                       VectorDouble& tabin,
                       VectorDouble& tabout);
  Id _vmap_general(Db *db,
                    Id radius,
                    const NamingConvention &namconv);
  Id _vmap_grid(DbGrid *dbgrid,
                 const NamingConvention &namconv);
  static Id _get_variable_order(Id nvar, Id ivar0, Id jvar0);
  static void _complexArrayAlloc(Id size, VectorVectorDouble& tab);
  static Id _vmap_load_simple(DbGrid* dbgrid,
                               Id ndim,
                               Id sizetot,
                               const Id* dims,
                               Id* dinv,
                               Id ivar,
                               Id jvar,
                               VectorVectorDouble& i1i2,
                               VectorVectorDouble& z1i2,
                               VectorVectorDouble& z2i1,
                               VectorVectorDouble& z2z1);
  static Id _vmap_load_cross(DbGrid* dbgrid,
                              Id ndim,
                              Id sizetot,
                              const Id* dims,
                              Id* dinv,
                              Id ivar,
                              Id jvar,
                              VectorVectorDouble& i1i1,
                              VectorVectorDouble& z1i1,
                              VectorVectorDouble& i2i2,
                              VectorVectorDouble& z2i2);
  static void _vmap_blank(VectorVectorDouble& tab);
  static void _product_conjugate(double coef,
                                 VectorVectorDouble& tab1,
                                 VectorVectorDouble& tab2,
                                 VectorVectorDouble& tab);
  static void _vmap_rescale(double scale, VectorDouble& tab1, VectorDouble& tab2);
  static void _vmap_shift(VectorDouble& tab, VectorDouble& tabm1, VectorDouble& tabm2);
  void _vmap_store(VectorDouble& tab, Id iptr);
  void _vmap_normalize(Id nv2);
  Id _findNeighCell(const VectorInt& indg0, const VectorInt& neigh, Id rank, VectorInt& indg1);

private:
  DbGrid* _dbmap; // Pointer to the already existing output DbGrid (not to be deleted)
};

GSTLEARN_EXPORT DbGrid* db_vmap(Db *db,
                                const ECalcVario &calcul_type = ECalcVario::fromKey("VARIOGRAM"),
                                const VectorInt &nxx = VectorInt(),
                                const VectorDouble &dxx = VectorDouble(),
                                Id radius = 0,
                                bool flag_FFT = true,
                                const NamingConvention &namconv = NamingConvention("VMAP"));
}
