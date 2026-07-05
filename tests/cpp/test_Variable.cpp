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
#include "Basic/VectorNumT.hpp"
#include "DataBase/DbCol.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main()
{
  DbData data{};

  data.addArray(VectorDouble{1., 2., 3.}, "hello");
  data.addArray(VectorInt{5, 6, 7}, "world");
  data.addArray(VectorString{"foo", "bar", "baz"}, "foobar");
  std::cout << data.getValue<double>(2, 0).value_or(-1.) << " "
            << data.getValue<int>(2, 1).value_or(-1) << " "
            << data.getValue<String>(2, 2).value_or("failed") << "\n";

  data.setValue(2, 0, 4);
  data.setValue(2, 1, 8);
  data.setValue(2, 2, "foobar");

  std::cout << data.getValue<double>(2, 0).value_or(-1.) << " "
            << data.getValue<int>(2, 1).value_or(-1) << " "
            << data.getValue<String>(2, 2).value_or("failed") << "\n";
  return 0;
}
