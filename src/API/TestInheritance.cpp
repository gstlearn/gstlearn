/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "API/TestInheritance.hpp"

#include <sstream>

TestInheritance::TestInheritance()
 : _iproj(nullptr)
{
}

TestInheritance::~TestInheritance()
{
}

String TestInheritance::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (nullptr != _iproj)
  {
    sstr << "Number of Apices: " << _iproj->getApexNumber()  << std::endl;
    sstr << "Number of Points: " << _iproj->getPointNumber();
  }
  else
  {
    sstr << "Projection Matrix not yet defined!";
  }
  return sstr.str();
}
