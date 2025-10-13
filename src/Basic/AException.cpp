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
#include "Basic/AException.hpp"

#include <iostream>
#include <sstream>

namespace gstlrn
{
AException::AException(const std::string& msg)
: std::exception()
, _msg(msg)
{
}

AException::~AException()
{
}

const char* AException::what() const noexcept
{
  return _msg.c_str();
}

void throw_exp(const std::string& msg,
               const std::string& file,
               Id line)
{
  std::stringstream sstr;
  if (!file.empty())
  {
    sstr << file;
    if (line > 0) sstr << "@" << line;
    sstr << ": ";
  }
  sstr << msg;
  std::cout << "Error: " << sstr.str() << std::endl;
  throw(AException(sstr.str()));
}
}