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
#include "geoslib_d.h"

#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"

class Db;
class ECalcVario;
class Polygons;

/**
 * \brief
 * Class containing the Variogram Cloud which uses an DbGrid provided by the user
 * This function simply calculate and add the results as new field in this DbGrid.
 */
class GSTLEARN_EXPORT VCloud : public AVario
{
public:
  VCloud(DbGrid *dbcloud, const VarioParam* varioparam);
  VCloud(const VCloud& r);
  VCloud& operator=(const VCloud& r);
  virtual ~VCloud();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(VCloud)

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVCloud Interface
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

  int compute(Db *db, const NamingConvention &namconv = NamingConvention("Cloud"));

  int selectFromPolygon(Db *db, Polygons *polygon, int idir = 0);

private:
  void _variogram_cloud(Db *db, int idir);
  void _final_discretization_grid();
  int  _update_discretization_grid(double x, double y);

private:
  DbGrid* _dbcloud; // Pointer to the already existing output DbGrid (not to be deleted)
  const VarioParam* _varioparam; // Pointer (not to be deleted)
};

GSTLEARN_EXPORT DbGrid* db_vcloud(Db *db,
                                  const VarioParam *varioparam,
                                  double lagmax = TEST,
                                  double varmax = TEST,
                                  int lagnb = 100,
                                  int varnb = 100,
                                  const NamingConvention &namconv = NamingConvention("Cloud"));

