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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

//Enums
#include "Anamorphosis/EAnam.hpp"

#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Anam : public AStringable
{
public:
  Anam(const EAnam& type = EAnam::UNDEFINED);
  Anam(const Anam &m);
  Anam& operator= (const Anam &m);
  virtual ~Anam();

  virtual String toString(int level = 0) const override;

  const EAnam&  getType() const { return _type; }

private:
  EAnam _type;
};
