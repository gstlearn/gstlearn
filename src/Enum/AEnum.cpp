/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
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
