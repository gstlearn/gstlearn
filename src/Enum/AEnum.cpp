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
#include "Enum/AEnum.hpp"
#include "Basic/AStringable.hpp"

#include <iostream>
#include <iomanip>

void AEnum::printEnum() const {
  _printMsg("  %2d - %11s : %s\n", _value, _key.c_str(), _descr.c_str());
}

/**
 * This function is used to call standard 'message' function
 * (in order to route the message to the relevant terminal)
 * @param format Printing format
 * @param args   Variable list of arguments
 */
template<typename ... Args>
void AEnum::_printMsg(const char *format, Args... args)
{
  message(format, args...);
}
