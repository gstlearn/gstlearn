/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include <Geometry/ABiTargetCheck.hpp>
#include "geoslib_f.h"


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
