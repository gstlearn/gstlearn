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

#include "Basic/VectorT.hpp"
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

    void setShift(Id shift) { _shift = shift; }

    void setSeed(Id seed) { _seed = seed; }

    void setNbSimu(Id nbsimu) { _nbsimu = nbsimu; }

    void saveResults(
      Db* db,
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

    Id _getShift() const { return _shift; }

    bool _getFlagBayes() const { return _flagBayes; }

    bool _getFlagDGM() const { return _flagDGM; }

    bool _getFlagGibbs() const { return _flagGibbs; }

    bool _getFlagPGS() const { return _flagPGS; }

    virtual Id _simulate() { return 0; }

    virtual Id _compute(
      Db* db,
      Id isimu,
      const VectorBool& activeArray,
      VectorVectorDouble& tab)
    {
      DECLARE_UNUSED(db);
      DECLARE_UNUSED(isimu);
      DECLARE_UNUSED(activeArray);
      DECLARE_UNUSED(tab);
      return 0;
    }

    virtual Id _getIcase() const { return 0; }

    void _computeGradient(
      Db* dbgrd,
      Id isimu,
      double delta,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _computeTangent(
      Db* dbtgt,
      Id isimu,
      double delta,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _correctMean(const VectorBool& activeArray, VectorVectorDouble& tab);
    void _convertToDifference(
      Db* dbin,
      Id isimu,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _updateDataToTarget(Db* dbin, Db* dbout) const;
    static void _allocateForOneSimulation(
      const Db* db,
      Id nvar,
      VectorBool& activeArray,
      VectorVectorDouble& tab,
      bool flagActive = true);
    Id _conditionalKriging(Db* dbin, Db* dbout);
    void _simulateNugget(
      Db* db,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);

  private:
    Id _nbsimu;
    Id _seed;
    Id _shift;
    bool _flagBayes;
    bool _flagPGS;
    bool _flagGibbs;
    bool _flagDGM;
    VectorInt _seedPerSimu;
  };

} // namespace gstlrn
