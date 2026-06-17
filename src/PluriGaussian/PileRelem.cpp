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
#include "PluriGaussian/PileRelem.hpp"
#include "PluriGaussian/CalcModelPGS.hpp"
#include "PluriGaussian/PileSplit.hpp"
#include "geoslib_old_f.h"

namespace gstlrn
{
  PileRelem::PileRelem(PileSplit* old_split)
    : _nrule(0)
    , _nbyrule(0)
    , _facies(0)
    , _Rrules()
    , _Rfipos()
    , _oldSplit(old_split)
    , _splits()
  {
  }

  std::unique_ptr<PileRelem> PileRelem::alloc(PileSplit* old_split)
  {
    return std::make_unique<PileRelem>(old_split);
  }

  void PileRelem::explore(bool verbose)
  {
    for (Id is = 0, ns = _getNSplit(); is < ns; is++)
    {
      PileSplit* split = _splits[is].get();
      if (split != nullptr) split->collapse(verbose);
      _glue(
        split->getNrule(), split->getNbyrule(), split->getSrules(),
        split->getSfipos());
      if (verbose) st_rules_print2("Relem", _nrule, _nbyrule, _Rrules, _Rfipos);
    }
  }

  void PileRelem::_glue(
    Id nrule1,
    Id nbyrule1,
    const VectorVectorInt& rules1,
    const VectorVectorInt& fipos1)
  {
    if (nrule1 <= 0) return;
    Id nnew = _nrule + nrule1;
    Id NRULE = getRuleAutoNRULE();
    Id NCOLOR = getRuleAutoNCOLOR();

    _Rrules.resize(nnew);
    for (Id i = _nrule; i < nnew; i++) _Rrules[i].resize(NRULE, 0);
    _Rfipos.resize(nnew);
    for (Id i = _nrule; i < nnew; i++) _Rfipos[i].resize(NCOLOR, 0);

    Id ir = _nrule;
    for (Id i1 = 0; i1 < nrule1; i1++, ir++)
    {
      for (Id ic = 0; ic < nbyrule1; ic++) _Rrules[ir][ic] = rules1[i1][ic];
      for (Id ic = 0; ic < NCOLOR; ic++) _Rfipos[ir][ic] = fipos1[i1][ic];
    }
    _nrule = nnew;
    _nbyrule = nbyrule1;
  }

  void
    PileRelem::define(const VectorInt& facies, Id side, const VectorInt& poss)
  {
    Id NCOLOR = getRuleAutoNCOLOR();
    Id nfacies = facies.size();
    Id number = 0;
    if (poss.empty())
      number = nfacies;
    else
    {
      for (Id i = 0; i < nfacies; i++)
        if (poss[i] == side) number++;
    }

    _facies.resize(number, 0);
    _Rfipos.resize(1);
    _Rfipos[0].resize(NCOLOR, 0);

    Id ecr = 0;
    for (Id i = 0; i < nfacies; i++)
    {
      if (poss.empty() || poss[i] == side) _facies[ecr++] = facies[i];
    }
    if (number == 1) _Rfipos[0][_facies[0] - 1] = 1;
  }

  /****************************************************************************
   **
   ** FUNCTION: subdivide
   **
   ** PURPOSE:  Define and store all the possibilities of splitting
   ** PURPOSE:  the list of facies of the current Relem
   **
   ** IN_ARGS: half    1 to take the symmetry into account
   ** IN_ARGS: noper   Number of underlying GRF
   **
   *****************************************************************************/
  void PileRelem::subdivide(Id half, Id noper)
  {
    bool verbose = false;
    Id NGRF = getRuleAutoNGRF();

    Id ncur = static_cast<Id>(_facies.size());
    if (ncur <= 1) return;

    _splits.clear();
    Id previous_oper = (_oldSplit == nullptr) ? 1 : _oldSplit->getOper();

    for (Id oper = 1; oper <= noper; oper++)
    {
      half = (oper == previous_oper);
      auto divs = ut_split_into_two(ncur, half, verbose);
      for (Id is = 0, ns = divs.size(); is < ns; is++)
      {
        _splits.push_back(std::make_unique<PileSplit>());
        PileSplit* split = _splits.back().get();
        split->setOper(oper);
        for (Id i = 0; i < 2; i++)
        {
          split->setRelem(i, PileRelem::alloc(split));
          split->getRelem(i)->define(_facies, 1 - i, divs[is]);
          split->getRelem(i)->subdivide(0, NGRF);
        }
      }
    }
  }

  void PileRelem::_fiposDecode(Id fipos, VectorInt& fgrf)
  {
    Id div;

    Id NGRF = getRuleAutoNGRF();
    Id nmax = 1 + NGRF;
    Id base = 2 * NGRF;
    for (Id i = 0; i < nmax; i++) fgrf[i] = 0;
    for (Id i = 0; i < nmax; i++)
    {
      fipos = fipos - 1;
      div = fipos / base;
      if (div > 0) fgrf[i] = fipos - base * div + 1;
      fipos = div;
    }
  }

  Id PileRelem::_fiposEncode(VectorInt& fgrf)
  {
    Id found, fipos;

    Id NGRF = getRuleAutoNGRF();
    Id nmax = 1 + NGRF;
    Id base = 2 * NGRF;

    found = fipos = 0;
    for (Id i = nmax - 1; i >= 0; i--)
    {
      if (fgrf[i] > 0) found++;
      if (found == 1) fipos = 1;
      fipos = fipos * base + fgrf[i];
    }
    return (fipos);
  }

  Id PileRelem::_permut(Id value, Id igrf)
  {
    if (igrf == 0)
    {
      if (value == 1) return (2);
      if (value == 2) return (1);
    }
    else if (igrf == 1)
    {
      if (value == 3) return (4);
      if (value == 4) return (3);
    }
    else if (igrf == 2)
    {
      if (value == 5) return (6);
      if (value == 6) return (5);
    }
    else
    {
      messageAbort("Function st_permut has been programmed up to 3 GRFs");
    }
    return (value);
  }

  Id PileRelem::updateOrientation(Id fac0, Id igrf_cas, VectorInt& fgrf)
  {
    Id facp, loc0, loc1;

    Id NGRF = getRuleAutoNGRF();
    Id nmax = 1 + NGRF;
    Id fac = fac0;
    if (fac0 < 0) return (fac0);

    /* Decomposition */

    _fiposDecode(fac, fgrf);

    /* Loop on the GRF permutations */

    loc0 = igrf_cas;
    for (Id igrf = 0; igrf < NGRF; igrf++)
    {
      loc1 = loc0 / 2;
      if (loc0 - 2 * loc1 > 0)
      {

        /* Update the orientation of 'grf' */

        for (Id i = 0; i < nmax; i++) fgrf[i] = _permut(fgrf[i], igrf);
      }
      loc0 = loc1;
    }

    /* Recomposition */

    facp = _fiposEncode(fgrf);

    return (facp);
  }

  Id PileRelem::sameScore(Id ir0, Id igrf_cas, VectorInt& fgrf, VectorInt& fcmp)
  {
    Id NCOLOR = getRuleAutoNCOLOR();
    VectorVectorInt& fipos = _Rfipos;
    if (ir0 <= 0) return (-1);

    // Modify the orientation of 'grf' for the current 'fipos'

    for (Id ic = 0; ic < NCOLOR; ic++)
      fcmp[ic] = updateOrientation(fipos[ir0][ic], igrf_cas, fgrf);

    // Look if the same 'fipos' has already been calculated

    for (Id ir = 0; ir < ir0; ir++)
    {
      bool flag_same = true;
      for (Id ic = 0; ic < NCOLOR && flag_same; ic++)
      {
        if (fipos[ir][ic] != fcmp[ic]) flag_same = false;
      }
      if (flag_same) return (ir);
    }
    return (-1);
  }

} // namespace gstlrn
