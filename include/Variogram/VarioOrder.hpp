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

#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT VarioOrder
  {
  public:
    VarioOrder(Id flag_dist = 0, Id size_aux = 0);
    VarioOrder(const VarioOrder& r) = default;
    VarioOrder& operator=(const VarioOrder& r) = default;
    virtual ~VarioOrder() = default;

    Id final();
    void clear();
    void printout(Id idir_target, Id ipas_target, Id verbose);
    void getBounds(Id idir, Id ilag, Id* ifirst, Id* ilast) const;
    void getIndices(Id ipair, Id* iech, Id* jech, double* dist) const;
    void getAuxiliary(Id ipair, char* aux_iech, char* aux_jech) const;
    Id add(
      Id iech,
      Id jech,
      void* aux_iech,
      void* aux_jech,
      Id ilag,
      Id idir,
      double dist);

    bool empty() const { return (_npair == 0); }

  private:
    Id _nalloc;
    Id _npair;
    Id _sizeAux;
    bool _flagDist;
    VectorInt _tabIech;
    VectorInt _tabJech;
    VectorInt _tabIpas;
    VectorInt _tabSort;
    std::vector<char> _tabAuxIech;
    std::vector<char> _tabAuxJech;
    VectorDouble _tabDist;
  };
} // namespace gstlrn
