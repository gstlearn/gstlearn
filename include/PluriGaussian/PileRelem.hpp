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
#include <memory>

namespace gstlrn
{
  /**
   * @brief This class is a Utility class used in the Automatic Rule Fitting procedure
   * It is not exported via SWIG
   * Moreover this class is not supposed to be duplicated.
   */
  class PileSplit;

  class GSTLEARN_EXPORT PileRelem
  {
  public:
    PileRelem(PileSplit* old_split = nullptr);
    PileRelem(const PileRelem& r) = delete;
    PileRelem& operator=(const PileRelem& r) = delete;
    virtual ~PileRelem() = default;

    static std::unique_ptr<PileRelem> alloc(PileSplit* old_split = nullptr);

    VectorInt& getFacies() { return _facies; }

    Id getFacies(Id i) const { return _facies[i]; }

    Id getNrule() const { return _nrule; }

    Id getNbyrule() const { return _nbyrule; }

    VectorVectorInt& getRrules() { return _Rrules; }

    VectorInt& getRrule(Id irule) { return _Rrules[irule]; }

    VectorVectorInt& getRfipos() { return _Rfipos; }

    VectorInt& getRfipos(Id iprod) { return _Rfipos[iprod]; }

    Id sameScore(Id ir0, Id igrf_cas, VectorInt& fgrf, VectorInt& fcmp);

    void explore(bool verbose = false);
    void define(
      const VectorInt& facies,
      Id side = ITEST,
      const VectorInt& poss = VectorInt());
    void subdivide(Id half, Id noper);

  private:
    static Id _permut(Id value, Id igrf);
    static Id _fiposEncode(VectorInt& fgrf);
    static void _fiposDecode(Id fipos, VectorInt& fgrf);
    static Id updateOrientation(Id fac0, Id igrf_cas, VectorInt& fgrf);

    Id _getNSplit() const { return _splits.size(); }

    void _glue(
      Id nrule1,
      Id nbyrule1,
      const VectorVectorInt& rules1,
      const VectorVectorInt& fipos1);

  private:
    Id _nrule; // Number of generated rules
    Id _nbyrule; // Number of symbols in the Rules
    VectorInt _facies; // List of facies
    VectorVectorInt _Rrules; // List of rules (Dim: [nitem][nrule])
    VectorVectorInt _Rfipos; // Position of facies (Dim: [nprod][NCOLOR])
    PileSplit* _oldSplit; // Not allocated
    std::vector<std::unique_ptr<PileSplit>> _splits;
  };

} // namespace gstlrn
