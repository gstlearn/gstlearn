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

class Db;
class ECalcVario;

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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVMap Interface
  double _getIVAR(const Db *db, int iech, int ivar) const override;
  void _setResult(int iech1,
                  int iech2,
                  int nvar,
                  int ipas,
                  int ivar,
                  int jvar,
                  int orient,
                  double ww,
                  double dist,
                  double value) override;

  int compute(Db *db,
              const ECalcVario &calcul_type,
              int radius,
              bool flag_FFT = true,
              const NamingConvention &namconv = NamingConvention("VMAP"));

private:
  int _grid_fft(DbGrid *dbgrid,
                const NamingConvention &namconv);
  static void _extract(const int* nxmap,
                       const int* nxgrid,
                       int* dims,
                       VectorDouble& tabin,
                       VectorDouble& tabout);
  int _vmap_general(Db *db,
                    int radius,
                    const NamingConvention &namconv);
  int _vmap_grid(DbGrid *dbgrid,
                 const NamingConvention &namconv);
  static int _get_variable_order(int nvar, int ivar0, int jvar0);
  static void _complexArrayAlloc(int size, VectorVectorDouble& tab);
  static int _vmap_load_simple(DbGrid* dbgrid,
                               int ndim,
                               int sizetot,
                               const int* dims,
                               int* dinv,
                               int ivar,
                               int jvar,
                               VectorVectorDouble& i1i2,
                               VectorVectorDouble& z1i2,
                               VectorVectorDouble& z2i1,
                               VectorVectorDouble& z2z1);
  static int _vmap_load_cross(DbGrid* dbgrid,
                              int ndim,
                              int sizetot,
                              const int* dims,
                              int* dinv,
                              int ivar,
                              int jvar,
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
  void _vmap_store(VectorDouble& tab, int iptr);
  void _vmap_normalize(int nv2);
  int _findNeighCell(const VectorInt& indg0, const VectorInt& neigh, int rank, VectorInt& indg1);

private:
  DbGrid* _dbmap; // Pointer to the already existing output DbGrid (not to be deleted)
};

GSTLEARN_EXPORT DbGrid* db_vmap(Db *db,
                                const ECalcVario &calcul_type = ECalcVario::fromKey("VARIOGRAM"),
                                const VectorInt &nxx = VectorInt(),
                                const VectorDouble &dxx = VectorDouble(),
                                int radius = 0,
                                bool flag_FFT = true,
                                const NamingConvention &namconv = NamingConvention("VMAP"));
