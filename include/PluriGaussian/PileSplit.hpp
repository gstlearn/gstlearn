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

#include "Basic/VectorNumT.hpp"
#include "PluriGaussian/PileRelem.hpp"
#include <array>
#include <memory>

namespace gstlrn
{
  /**
   * @brief This class is a Utility class used in the Automatic Rule Fitting procedure
   * It is not exported via SWIG.
   * Moreover this class is not supposed to be duplicated.
   */

  class GSTLEARN_EXPORT PileSplit
  {
  public:
    PileSplit();
    PileSplit(const PileSplit& r) = delete;
    PileSplit& operator=(const PileSplit& r) = delete;
    virtual ~PileSplit() = default;

    static std::unique_ptr<PileSplit> alloc();

    Id getOper() const { return _oper; }

    Id getNrule() const { return _nrule; }

    Id getNbyrule() const { return _nbyrule; }

    VectorVectorInt& getSrules() { return _Srules; }

    VectorVectorInt& getSfipos() { return _Sfipos; }

    PileRelem* getRelem(Id i) const { return _relems[i].get(); }

    void setOper(Id oper) { _oper = oper; }

    void setRelem(Id i, std::unique_ptr<PileRelem> relem)
    {
      _relems[i] = std::move(relem);
    }

    void collapse(bool verbose = false);

  private:
    static Id _defineFipos(Id oper, Id side);

    Id _getNRelem() const { return _relems.size(); }

    void _ruleProduct(
      Id nprod,
      Id nrule1,
      Id nbyrule1,
      const VectorVectorInt& rules1,
      const VectorVectorInt& fipos1,
      Id nrule2,
      Id nbyrule2,
      const VectorVectorInt& rules2,
      const VectorVectorInt& fipos2,
      bool verbose = false);

  private:
    Id _oper; // Rank of operator
    Id _nrule; // Number of generated rules
    Id _nbyrule; // Number of symbols in the Rules
    VectorVectorInt _Srules; // List of rules (Dim: [nitem][NRULE])
    VectorVectorInt _Sfipos; // Position of facies (Dim: [nprod][NCOLOR])
    std::array<std::unique_ptr<PileRelem>, 2> _relems;
  };

  GSTLEARN_EXPORT void st_rule_print2(
    Id rank,
    Id nbyrule,
    const VectorVectorInt& rules,
    const VectorVectorInt& fipos,
    bool flag_rank,
    Id flag_similar,
    Id flag_igrf,
    double score);
  GSTLEARN_EXPORT void st_rules_print2(
    const char* title,
    Id nrule,
    Id nbyrule,
    const VectorVectorInt& rules,
    const VectorVectorInt& fipos);

} // namespace gstlrn
