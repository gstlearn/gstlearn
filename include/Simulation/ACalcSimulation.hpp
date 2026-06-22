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

    bool isConditional() const { return getDbin() != nullptr; }

    Id getSeed() const { return _seed; }

    Id getNbSimu() const { return _nbsimu; }

    Id getNVar() const;

    void setShift(Id shift) { _shift = shift; }

    void setSeed(Id seed) { _seed = seed; }

    void setNbSimu(Id nbsimu) { _nbsimu = nbsimu; }

    Id getSeedPerSimu(Id isimu) const { return _seedPerSimu[isimu]; }

  protected:
    bool _check() override;

    Id _getShift() const { return _shift; }

    virtual Id _getIcase() const { return 0; }

    static void _allocateForOneSimulation(
      const Db* db,
      Id nvar,
      VectorBool& activeArray,
      VectorVectorDouble& tab,
      bool flagActive = true);

  private:
    Id _nbsimu;
    Id _seed;
    Id _shift;
    VectorInt _seedPerSimu;
  };

} // namespace gstlrn
