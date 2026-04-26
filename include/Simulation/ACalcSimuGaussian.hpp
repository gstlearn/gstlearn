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

#include "Simulation/ACalcSimulation.hpp"

namespace gstlrn
{
  class CovBase;

  class GSTLEARN_EXPORT ACalcSimuGaussian: public ACalcSimulation
  {
  public:
    ACalcSimuGaussian(Id nbsimu = 1, Id seed = 4324324, bool verbose = false);
    ACalcSimuGaussian(const ACalcSimuGaussian& r) = delete;
    ACalcSimuGaussian& operator=(const ACalcSimuGaussian& r) = delete;
    virtual ~ACalcSimuGaussian();

    Id getNVar() const;

    void saveResults(
      Db* db,
      Id isimu,
      const VectorBool& activeArray,
      const VectorVectorDouble& tab) const;

    void setFlagBayes(bool flag_bayes) { _flagBayes = flag_bayes; }

    void setFlagDGM(bool flag_dgm) { _flagDGM = flag_dgm; }

    void setFlagGibbs(bool flag_gibbs) { _flagGibbs = flag_gibbs; }

    void setFlagPGS(bool flag_pgs) { _flagPGS = flag_pgs; }

    void setFlagAllocationAlreadyDone(bool flag)
    {
      _flagAllocationAlreadyDone = flag;
    }

  protected:
    bool _check() override;
    bool _preprocess() override;
    bool _postprocess() override;
    bool _run() override;

    void _initializeOutputAttribute();

    bool _getFlagBayes() const { return _flagBayes; }

    bool _getFlagDGM() const { return _flagDGM; }

    bool _getFlagGibbs() const { return _flagGibbs; }

    bool _getFlagPGS() const { return _flagPGS; }

    bool _isConditional() const { return _flagCond; }

    virtual Id _simulate(Id isimu)
    {
      DECLARE_UNUSED(isimu);
      return 0;
    }

    virtual Id
      _compute(Db* db, const VectorBool& activeArray, VectorVectorDouble& tab)
    {
      DECLARE_UNUSED(db);
      DECLARE_UNUSED(activeArray);
      DECLARE_UNUSED(tab);
      return 0;
    }

    bool _isAllocationAlreadyDone() const { return _flagAllocationAlreadyDone; }

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
      Id isimu,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _updateDataToTarget() const;
    static void _allocateForOneSimulation(
      const Db* db,
      Id nvar,
      VectorBool& activeArray,
      VectorVectorDouble& tab,
      bool flagActive = true);
    Id _conditionalKriging();
    void _simulateNugget(
      Db* db,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);

  private:
    Id _iattOut;
    bool _flagCond;
    bool _flagBayes;
    bool _flagPGS;
    bool _flagGibbs;
    bool _flagDGM;
    bool _flagAllocationAlreadyDone;
  };

} // namespace gstlrn
