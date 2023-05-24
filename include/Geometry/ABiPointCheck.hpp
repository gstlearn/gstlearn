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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT ABiPointCheck: public AStringable, public ICloneable
{
public:
  ABiPointCheck();
  ABiPointCheck(const ABiPointCheck& r);
  ABiPointCheck& operator=(const ABiPointCheck& r);
  virtual ~ABiPointCheck();

  virtual bool isOK(const SpacePoint &P1,
                    const SpacePoint &P2,
                    int iech1 = -1,
                    int iech2 = -1) const = 0;
};
