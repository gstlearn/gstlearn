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
#include "Basic/ASerializable.hpp"
#include "Basic/VectorNumT.hpp"
#include "DataBase/DbData.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_DataBase-");

  DbData data{};

  data.addArray(VectorDouble{1., 2., 3.}, "hello");
  data.addArray(VectorInt{5, 6, 7}, "world");
  data.addArray(VectorString{"foo", "bar", "baz"}, "foobar");

  message("Initial values\n");
  message("%lf\n", data.getValue<double>(2, 0).value_or(-1.));
  message("%d\n", data.getValue<int>(2, 1).value_or(-1));
  message("%s\n", data.getValue<String>(2, 2).value_or("failed").c_str());

  data.setValue(4., 2, 0);
  data.setValue(8, 2, 1);
  data.setValue("foobar", 2, 2);

  message("Values after modification\n");
  message("%lf\n", data.getValue<double>(2, 0).value_or(-1.));
  message("%d\n", data.getValue<int>(2, 1).value_or(-1));
  message("%s\n", data.getValue<String>(2, 2).value_or("failed").c_str());
  return 0;
}
