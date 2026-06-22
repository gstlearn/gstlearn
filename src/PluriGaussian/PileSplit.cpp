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
#include "PluriGaussian/PileSplit.hpp"
#include "Basic/Message.hpp"
#include "PluriGaussian/CalcModelPGS.hpp"

namespace gstlrn
{
  PileSplit::PileSplit()
    : _oper(0)
    , _nrule(0)
    , _nbyrule(0)
    , _Srules()
    , _Sfipos()
    , _relems{nullptr, nullptr}
  {
  }

  std::unique_ptr<PileSplit> PileSplit::alloc()
  {
    return std::make_unique<PileSplit>();
  }

  void PileSplit::collapse(bool verbose)
  {
    VectorVectorInt ptr[2];
    VectorVectorInt ruleLocal(1);
    ruleLocal[0].resize(1);

    // Explore the two Relems
    for (Id i = 0; i < 2; i++)
    {
      PileRelem* relem = _relems[i].get();
      if (relem == nullptr) continue;
      relem->explore(verbose);
    }

    // Prepare collapsing

    Id num[2] = {0, 0};
    Id nby[2] = {0, 0};
    Id nprod = 0;
    if (_nrule <= 0)
    {
      for (Id i = 0; i < 2; i++)
      {
        PileRelem* relem = _relems[i].get();
        if (relem == nullptr) continue;
        if (relem->getFacies().size() <= 1)
        {
          ruleLocal[0][0] = relem->getFacies(0);
          num[i] = 1;
          nby[i] = 1;
          ptr[i] = ruleLocal;
        }
        else
        {
          num[i] = relem->getNrule();
          nby[i] = relem->getNbyrule();
          ptr[i] = relem->getRrules();
        }
      }

      // Merge the rules of the two Relem (by product)

      nprod = num[0] * num[1];
      _nbyrule = nby[0] + nby[1] + 1;
      if (nprod > 0)
      {
        _ruleProduct(
          nprod, num[0], nby[0], ptr[0], _relems[0]->getRfipos(), num[1],
          nby[1], ptr[1], _relems[1]->getRfipos(), verbose);
        if (verbose)
          st_rules_print2("Split", _nrule, _nbyrule, _Srules, _Sfipos);
      }
    }
  }

  Id PileSplit::_defineFipos(Id oper, Id side)
  {
    Id reponse;

    if (isNA(side))
      reponse = 1;
    else
      reponse = 2 * (oper - 1) + (1 - side) + 1;

    return (reponse);
  }

  void PileSplit::_ruleProduct(
    Id nprod,
    Id nrule1,
    Id nbyrule1,
    const VectorVectorInt& rules1,
    const VectorVectorInt& fipos1,
    Id nrule2,
    Id nbyrule2,
    const VectorVectorInt& rules2,
    const VectorVectorInt& fipos2,
    bool verbose)
  {
    Id NRULE = getRuleAutoNRULE();
    Id NGRF = getRuleAutoNGRF();
    Id NCOLOR = getRuleAutoNCOLOR();
    Id BASE = 2 * NGRF;
    _Srules.resize(nprod);
    for (Id i = 0; i < nprod; i++) _Srules[i].resize(NRULE, 0);
    _Sfipos.resize(nprod);
    for (Id i = 0; i < nprod; i++) _Sfipos[i].resize(NCOLOR, 0);

    Id ir = 0;
    for (Id i1 = 0; i1 < nrule1; i1++)
      for (Id i2 = 0; i2 < nrule2; i2++, ir++)
      {
        Id ic = 0;
        _Srules[ir][ic++] = 1000 + _oper;

        if (verbose)
        {
          message("Rule Product (with operator %d)\n", _oper);
          st_rule_print2(i1, nbyrule1, rules1, fipos1, false, -1, -1, TEST);
          st_rule_print2(i2, nbyrule2, rules2, fipos2, false, -1, -1, TEST);
        }

        for (Id i = 0; i < nbyrule1; i++) _Srules[ir][ic++] = rules1[i1][i];
        for (Id i = 0; i < nbyrule2; i++) _Srules[ir][ic++] = rules2[i2][i];
        for (Id i = 0; i < NCOLOR; i++)
        {
          if (fipos1[i1][i] > 0)
            _Sfipos[ir][i] = fipos1[i1][i] * BASE + _defineFipos(_oper, 1);
          if (fipos2[i2][i] > 0)
            _Sfipos[ir][i] = fipos2[i2][i] * BASE + _defineFipos(_oper, 0);
        }

        if (verbose)
        {
          message("Product result=");
          st_rule_print2(
            ir, nbyrule1 + nbyrule2 + 1, _Srules, _Sfipos, false, -1, -1, TEST);
        }
      }
    _nrule = nprod;
    _nbyrule = nbyrule1 + nbyrule2 + 1;
  }

  void st_rule_print2(
    Id rank,
    Id nbyrule,
    const VectorVectorInt& rules,
    const VectorVectorInt& fipos,
    bool flag_rank,
    Id flag_similar,
    Id flag_igrf,
    double score)
  {
    Id NCOLOR = getRuleAutoNCOLOR();
    Id NGRF = getRuleAutoNGRF();

    // Print the Rank (optional)

    if (flag_rank) message("%4d:", rank + 1);

    // Print the Rule

    for (Id ic = 0; ic < nbyrule; ic++)
    {
      Id value = rules[rank][ic];
      if (value == 1001)
        message("  S");
      else if (value == 1002)
        message("  T");
      else
        message(" %2d", value);
    }

    // Print the score (if available)

    if (!isNA(score)) message(" -> %lf", score);

    // Print the Facies rank

    message(" (");
    for (Id ic = 0; ic < NCOLOR; ic++) message(" %3d", fipos[rank][ic]);
    message(" )");

    // Print the similar score

    if (flag_similar >= 0)
    {
      message(" [Same as: %3d (", flag_similar + 1);

      Id loc0 = flag_igrf;
      for (Id igrf = 0; igrf < NGRF; igrf++)
      {
        Id loc1 = loc0 / 2;
        if (loc0 - 2 * loc1 > 0) message(" by Symmetry around G%d", igrf + 1);
        loc0 = loc1;
      }
      message(" )]");
    }

    // Print the end-of-line

    message("\n");
  }

  void st_rules_print2(
    const char* title,
    Id nrule,
    Id nbyrule,
    const VectorVectorInt& rules,
    const VectorVectorInt& fipos)
  {
    if (nrule <= 0) return;
    message("%s (Nrule=%d, Nbyrule=%d):\n", title, nrule, nbyrule);
    for (Id ir = 0; ir < nrule; ir++)
      st_rule_print2(ir, nbyrule, rules, fipos, false, -1, -1, TEST);
  }
} // namespace gstlrn
