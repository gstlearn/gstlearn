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

