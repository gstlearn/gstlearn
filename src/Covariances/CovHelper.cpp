/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Covariances/CovHelper.hpp"
#include "Basic/VectorT.hpp"
#include "Enum/ECov.hpp"

VectorString CovHelper::getAllCovariances()
{
  VectorString vs;
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    vs.push_back(it.getKey());
    it.toNext();
  }
  return vs;
}
