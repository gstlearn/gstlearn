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
#include "Enum/AEnum.hpp"
#include "Basic/AStringable.hpp"

void AEnum::printEnum() const {
  _printMsg("  %2d - %11s : %s\n", _value, _key.data(), _descr.data());
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
