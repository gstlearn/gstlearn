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
    ACalcSimulation(Id nbsimu = 1, Id seed = 4324324, bool verbose = false);
    ACalcSimulation(const ACalcSimulation& r) = delete;
    ACalcSimulation& operator=(const ACalcSimulation& r) = delete;
    virtual ~ACalcSimulation();

    Id getSeed() const { return _seed; }

    Id getNbSimu() const { return _nbsimu; }

    Id getNVar() const;

    void setSeed(Id seed) { _seed = seed; }

    void setNbSimu(Id nbsimu) { _nbsimu = nbsimu; }

    void saveResults(Db* db,
                     Id shift,
                     Id isimu,
                     const VectorBool& activeArray,
                     const VectorVectorDouble& tab) const;
    void scaleAndSaveResults(Db* db,
                             const CovBase* cova,
                             Id shift,
                             Id isimu,
                             const VectorBool& activeArray,
                             const VectorVectorDouble& tab) const;

    void setFlagBayes(bool flag_bayes) { _flagBayes = flag_bayes; }

    void setFlagDGM(bool flag_dgm) { _flagDGM = flag_dgm; }

    void setFlagGibbs(bool flag_gibbs) { _flagGibbs = flag_gibbs; }

    void setFlagPGS(bool flag_pgs) { _flagPGS = flag_pgs; }

    Id getSeedPerSimu(Id isimu) const { return _seedPerSimu[isimu]; }

  protected:
    bool _check() override;
    bool _preprocess() override;

    Id _getSeedPerSimu(Id isimu) const;

    bool _getFlagBayes() const { return _flagBayes; }

    bool _getFlagDGM() const { return _flagDGM; }

    bool _getFlagGibbs() const { return _flagGibbs; }

    bool _getFlagPGS() const { return _flagPGS; }

    virtual Id _computeTB(Db* db)
    {
      DECLARE_UNUSED(db);
      return 0;
    };

    virtual Id _getIcase() const { return 0; }

    virtual void _setShift(Id shift) { DECLARE_UNUSED(shift); }

    void _computeGradient(Db* dbgrd, double delta);
    void _computeTangent(Db* dbtgt, double delta);
    void _correctStationaryMean(Db* dbout);
    void _difference(Db* dbin);
    void _updateDataToTarget(Db* dbin, Db* dbout) const;
    Id _checkGaussianDataToGrid(Db* dbin, Db* dbout) const;
    Id _conditionalKriging(Db* dbin, Db* dbout);

  private:
    Id _nbsimu;
    Id _seed;
    bool _flagCheck;
    bool _flagBayes;
    bool _flagPGS;
    bool _flagGibbs;
    bool _flagDGM;
    VectorInt _seedPerSimu;
  };

} // namespace gstlrn