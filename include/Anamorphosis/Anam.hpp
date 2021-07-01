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
#pragma once

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_enum.h"

class Anam : public AStringable
{
private:
  int _type;

public:
  Anam(int type = ANAM_UNDEFINED);
  Anam(const Anam &m);
  Anam& operator= (const Anam &m);
  virtual ~Anam();

  virtual String toString(int level = 0) const override;

  int  getType() const { return _type; }
  void setType(int type) { _type = type; }
};
