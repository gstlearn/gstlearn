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
#include "Basic/IClonable.hpp"

class GSTLEARN_EXPORT Anam : public AStringable, public IClonable
{
public:
  Anam(const EAnam& type = EAnam::UNDEFINED);
  Anam(const Anam &m);
  Anam& operator= (const Anam &m);
  virtual ~Anam();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual IClonable* clone() const override { return new Anam(*this); }

  const EAnam&  getType() const { return _type; }

private:
  EAnam _type;
};
