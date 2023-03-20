/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Basic/AException.hpp"

#include <iostream>
#include <sstream>

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
               int line)
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
