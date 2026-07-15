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
#include "Basic/File.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/DbCol.hpp"

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

  DbData data{};

  data.AddArray(VectorDouble{1., 2., 3.}, "hello");
  data.AddArray(VectorInt{5, 6, 7}, "world");
  data.AddArray(VectorString{"foo", "bar", "baz"}, "foobar");
  data.AddArray(VectorBool{true, false, true}, "gstlearn");

  Dictionary dict{{{1, "red"}, {3, "blue"}, {7, "green"}}};
  data.AddArray(VectorCategory{3, dict}, "categories");

  std::cout << data.getValue<double>(2, 0).value_or(-1.) << " "
            << data.getValue<int>(2, 1).value_or(-1) << " "
            << data.getValue<String>(2, 2).value_or("failed") << " "
            << data.getValue<bool>(2, 3).value_or(false) << " "
            << data.getValue<Dictionary::Category>(2, 4)
                 .value_or(Dictionary::Category{-1, "NOK"})
                 .second
            << "\n";

  data.setValue(2, 0, 0, 4);
  data.setValue(2, 1, 0, 8);
  data.setValue(2, 2, 0, "foobar");
  data.setValue(2, 3, 0, false);
  data.setValue(2, 4, 0, Dictionary::Category{1, "red"});

  std::cout << data.getValue<double>(2, 0).value_or(-1.) << " "
            << data.getValue<int>(2, 1).value_or(-1) << " "
            << data.getValue<String>(2, 2).value_or("failed") << " "
            << data.getValue<bool>(2, 3).value_or(false) << " "
            << data.getValue<Dictionary::Category>(2, 4)
                 .value_or(Dictionary::Category{-1, "NOK"})
                 .second
            << "\n";
  return 0;
}
