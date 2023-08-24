/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include <exception>
#include <string>

class GSTLEARN_EXPORT AException : public std::exception
{
public:
  AException(const std::string& msg = "");
  virtual ~AException();

  const char* what() const noexcept;

private:
  std::string _msg;
};

GSTLEARN_EXPORT void throw_exp(const std::string& msg = "",
                               const std::string& file = "",
                               int line = 0);

#define my_throw(msg) throw_exp(msg, __FILE__, __LINE__)

