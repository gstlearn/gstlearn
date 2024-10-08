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
#include "Covariances/ACovGradient.hpp"

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

