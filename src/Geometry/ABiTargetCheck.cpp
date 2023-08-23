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
#include "geoslib_f.h"

#include "Geometry/ABiTargetCheck.hpp"

ABiTargetCheck::ABiTargetCheck()
    : AStringable()
{
}

ABiTargetCheck::ABiTargetCheck(const ABiTargetCheck &r)
    : AStringable(r)
{
}

ABiTargetCheck& ABiTargetCheck::operator=(const ABiTargetCheck &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
  }
  return *this;
}

ABiTargetCheck::~ABiTargetCheck()
{
}
