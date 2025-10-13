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
#include "Covariances/NoStatFunctional.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
NoStatFunctional::NoStatFunctional(const AFunctional* func)
  : _func(func)
{
}

NoStatFunctional::~NoStatFunctional()
{
}

void NoStatFunctional::_informField(const VectorVectorDouble& coords,
                                    VectorDouble& tab,
                                    bool verbose)
{
  DECLARE_UNUSED(verbose)
  Id size = static_cast<Id>(coords[0].size());
  Id ndim = static_cast<Id>(coords.size());
  VectorDouble vec(ndim);
  for (Id icoords = 0; icoords < size; icoords++)
  {
    for (Id idim = 0; idim < ndim; idim++)
    {
      vec[idim] = coords[idim][icoords];
    }
    tab[icoords] = _func->getFunctionValue(vec);
  }
}

String NoStatFunctional::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt)
  std::stringstream sstr;
  if (_func == nullptr) return sstr.str();
  sstr << ANoStat::toString(strfmt);
  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
    sstr << "Functional" << std::endl;
  return sstr.str();
}
} // namespace gstlrn