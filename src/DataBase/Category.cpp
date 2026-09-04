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
#include "DataBase/Category.hpp"

#include <iostream>

namespace gstlrn
{
#ifndef SWIG
  std::ostream& operator<<(std::ostream& os, const Category& cat)
  {
    os << "Category(id=" << cat.getId() << ", name=\"" << cat.getName()
       << "\")";
    return os;
  }
#endif
} // namespace gstlrn
