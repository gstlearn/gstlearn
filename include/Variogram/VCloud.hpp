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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
  class Db;
  class ECalcVario;
  class Polygons;

  /**
   * \brief
   * Class containing the Variogram Cloud which uses an DbGrid provided by the user
   * This function simply calculates and adds the results as new field in this DbGrid.
   *
   * Some improvements:
   * - you can define a set of 'samples' so that the calculation compare all samples of 'db'
   *   with these specific samples only. Use method: 'setSamples'
   *   This is useful when you want to compare a set of 'reference' samples with all the others.
   * - you can provide a polygon in order to issue the scores of the Samples included
   *   in this polygon
   * - you can provide a maximum distance ('dismax') and a minimum variability ('varmin')
   *   in order to keep the description of the pair of Samples included in this range
   * These two selections are exclusive.
   * For these two cases, a specific printout is provided after the run, summarizing
   * the samples which are most frequently involved in this selective calculations.
   */
  class GSTLEARN_EXPORT VCloud: public AVario
  {
  public:
    VCloud(DbGrid* dbcloud, const VarioParam* varioparam);
    VCloud(const VCloud& r);
    VCloud& operator=(const VCloud& r);
    virtual ~VCloud();

  public:
    /// ICloneable interface
    IMPLEMENT_CLONING(VCloud)

    /// Has a specific implementation in the Target language
    DECLARE_TOTL;

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// AVCloud Interface
    double _getIVAR(const Db* db, Id iech, Id ivar) const override;
    void _setAVarioResult(
      Id iech1,
      Id iech2,
      Id nvar,
      Id idir,
      Id ilag,
      Id ivar,
      Id jvar,
      Id orient,
      double ww,
      double w1,
      double w2,
      double z1,
      double z2,
      double dist,
      double value) override;

    Id compute(
      Db* db,
      const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
      bool flag_ergodic = true,
      const NamingConvention& namconv = NamingConvention("Cloud"));

    void setPolygon(const Polygons* polygon);
    void setIntervals(double distmax, double varmin, Id countmax);
    void setSamples(const VectorInt& samples);

  private:
    void _variogramCloud(Db* db, Id idir);
    void _variogramCloudBySamples(Db* db, Id idir);
    void _finalDiscretizationOnGrid();
    Id _getDiscretizedCellRank(double x, double y);

    bool _flagBySamples() const { return !_samples.empty(); }

    bool _flagByPolygon() const { return _polygon != nullptr; }

    bool _flagByIntervals() const { return !isNA(_distmax) || !isNA(_varmin); }

  private:
    DbGrid* _dbcloud; // Pointer to DbGrid (not to be deleted)
    const VarioParam* _varioparam; // Pointer (not to be deleted)
    Id _IPTR;
    VectorInt _samples;
    bool _flagSelection;
    VectorInt _IDS;
    const Polygons* _polygon; // Pointer (not to be deleted)
    double _distmax;
    double _varmin;
    Id _countmax;
  };

} // namespace gstlrn
