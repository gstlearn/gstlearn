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
#include "Anamorphosis/Anam.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"

Anam::Anam(const EAnam& type)
    : AStringable(),
      _type(type)
{
}

Anam::Anam(const Anam &m)
    : _type(m._type)
{

}

Anam& Anam::operator=(const Anam &m)
{
  if (this != &m)
  {
    _type = m._type;
  }
  return *this;
}

Anam::~Anam()
{

}

String Anam::toString(int level) const
{
  std::stringstream sstr;
  sstr << toTitle(1, "Anamorphosis characteristics");
  return sstr.str();
}

