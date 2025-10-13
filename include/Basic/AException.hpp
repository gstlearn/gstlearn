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
#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <exception>
#include <string>

namespace gstlrn
{
 
class GSTLEARN_EXPORT AException : public std::exception
{
public:
  AException(const std::string& msg = "");
  virtual ~AException();

  const char* what() const noexcept;

private:
  std::string _msg;
};

[[noreturn]]
GSTLEARN_EXPORT void throw_exp(const std::string& msg = "",
                               const std::string& file = "",
                               Id line = 0);

#define my_throw(msg) throw_exp(msg, __FILE__, __LINE__)

#define my_throw_impossible(msg) throw_exp(msg, __FILE__, __LINE__)
} // namespace gstlrn
