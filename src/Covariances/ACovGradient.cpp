/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_f_private.h"

#include "Covariances/ACovGradient.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Drifts/ADriftElem.hpp"

#include <math.h>

ACovGradient::ACovGradient(const ECov& type, const CovContext& ctxt)
    : CovAniso(type, ctxt)
{
}

ACovGradient::ACovGradient(const ACovGradient &r)
    : CovAniso(r)
{
}

ACovGradient::ACovGradient(const CovAniso &r)
    : CovAniso(r)
{
}

ACovGradient& ACovGradient::operator=(const ACovGradient &r)
{
  if (this != &r)
  {
    CovAniso::operator =(r);
  }
  return *this;
}

ACovGradient::~ACovGradient()
{
}

