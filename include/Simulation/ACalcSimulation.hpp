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

#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{
  class CovBase;

  class GSTLEARN_EXPORT ACalcSimulation: public ACalcInterpolator
  {
  public:
    ACalcSimulation(Id nbsimu = 1, Id seed = 4324324);
    ACalcSimulation(const ACalcSimulation& r) = delete;
    ACalcSimulation& operator=(const ACalcSimulation& r) = delete;
    virtual ~ACalcSimulation();

    Id getSeed() const { return _seed; }

    Id getNbSimu() const { return _nbsimu; }

    Id getNVar() const;

    void setSeed(Id seed) { _seed = seed; }

    void setNbSimu(Id nbsimu) { _nbsimu = nbsimu; }

    void saveResults(Db* db,
                     Id icase,
                     Id shift,
                     Id isimu,
                     const VectorBool& activeArray,
                     const VectorVectorDouble& tab) const;
    void scaleAndSaveResults(Db* db,
                             const CovBase* cova,
                             Id icase,
                             Id shift,
                             Id isimu,
                             const VectorBool& activeArray,
                             const VectorVectorDouble& tab) const;

  protected:
    bool _check() override;
    bool _preprocess() override;

    virtual Id _computeTB(Db* db, Id icase, Id shift)
    {
      DECLARE_UNUSED(db);
      DECLARE_UNUSED(icase);
      DECLARE_UNUSED(shift);
      return 0;
    };

    void _computeGradient(Db* dbgrd, double delta);
    void _computeTangent(Db* dbtgt, double delta);
    void _correctStationaryMean(Db* dbout, Id icase, bool flagBayes);
    void _difference(Db* dbin,
                     Id icase,
                     bool flag_pgs = false,
                     bool flag_gibbs = false,
                     bool flag_dgm = false);
    void _updateDataToTarget(Db* dbin,
                             Db* dbout,
                             Id icase,
                             bool flag_pgs = false,
                             bool flag_dgm = false) const;
    Id _checkGaussianDataToGrid(Db* dbin, Db* dbout) const;
    Id _conditionalKriging(Db* dbin,
                           Db* dbout,
                           Id icase,
                           bool flag_bayes,
                           bool flag_dgm);

  private:
    Id _nbsimu;
    Id _seed;
  };

} // namespace gstlrn
