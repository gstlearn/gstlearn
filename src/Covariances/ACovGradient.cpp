/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f_private.h"

#include "Covariances/ACovGradient.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
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

