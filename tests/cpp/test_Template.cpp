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
#include "Basic/Message.hpp"
#include "Basic/VectorHelper.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to test the Template of class VectorNumT
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // ===========
  // Real values
  // ===========
  mestitle(1, "Case of VectorDouble");
  VectorDouble rtab1 = {0.3, 9.1, 1.2, TEST, 4.3};
  rtab1.dump("'rtab1'", false);

  message("Count of elements = %d\n", rtab1.count());
  message("Count of elements Defined = %d\n", rtab1.count(1));
  message("Count of elements Undefined = %d\n", rtab1.count(-1));
  message("Sum of elements = %lf\n", rtab1.sum());
  message("Product of elements = %lf\n", rtab1.prod());
  message("Mean = %lf\n", rtab1.mean());
  message("Median = %lf\n", rtab1.median());
  message("Minimum = %lf\n", rtab1.minimum());
  message("Maximum = %lf\n", rtab1.maximum());
  message("Variance = %lf\n", rtab1.variance());
  message("St. Dev. = %lf\n", rtab1.stdv());
  message("Norm = %lf\n", rtab1.norm());

  VectorDouble rtab2 = {5.9, TEST, -1.6, 2.2, -1.};
  rtab2.dump("'rtab2'", false);
  VectorDouble rtab3;

  rtab1 += rtab2;
  rtab1.dump("'rtab1' + 'rtab2'", false);
  rtab1 += 12;
  rtab1.dump("'rtab1' + 12", false);
  rtab3 = rtab1 + rtab2;
  rtab3.dump("'rtab1' + 'rtab2' -> 'rtab3''", false);

  rtab1 -= rtab2;
  rtab1.dump("'rtab1' - 'rtab2'", false);
  rtab1 -= 12;
  rtab1.dump("'rtab1' - 12", false);
  rtab3 = rtab1 - rtab2;
  rtab3.dump("'rtab1' - 'rtab2' -> 'rtab3''", false);

  rtab1 *= rtab2;
  rtab1.dump("'rtab1' * 'rtab2'", false);
  rtab1 *= 12;
  rtab1.dump("'rtab1' * 12", false);
  rtab3 = rtab1 * rtab2;
  rtab3.dump("'rtab1' * 'rtab2' -> 'rtab3''", false);

  rtab1 /= rtab2;
  rtab1.dump("'rtab1' / 'rtab2'", false);
  rtab1 /= 12;
  rtab1.dump("'rtab1' / 12", false);
  rtab3 = rtab1 / rtab2;
  rtab3.dump("'rtab1' / 'rtab2' -> 'rtab3''", false);

  // ==============
  // Integer values
  // ==============
  mestitle(1, "Case of VectorInt");
  VectorInt itab1 = {-2, 9, 1, ITEST, 4};
  itab1.dump("'itab1'", false);

  message("Count of elements = %d\n", itab1.count());
  message("Count of elements Defined = %d\n", itab1.count(1));
  message("Count of elements Undefined = %d\n", itab1.count(-1));
  message("Sum of elements = %lf\n", itab1.sum());
  message("Product of elements = %lf\n", itab1.prod());
  message("Minimum = %lf\n", itab1.minimum());
  message("Maximum = %lf\n", itab1.maximum());

  VectorInt itab2 = {5, ITEST, -1, 2, -1};
  itab2.dump("'itab2'", false);
  VectorInt itab3;

  itab1 += itab2;
  itab1.dump("'itab1' + 'itab2''", false);
  itab1 += 12;
  itab1.dump("'itab1' + 12", false);
  itab3 = itab1 + itab2;
  itab3.dump("'itab1'' + 'itab2' -> 'itab3'", false);

  itab1 -= itab2;
  itab1.dump("'itab1' - 'itab2''", false);
  itab1 -= 12;
  itab1.dump("'itab1' - 12", false);
  itab3 = itab1 - itab2;
  itab3.dump("'itab1'' - 'itab2' -> 'itab3'", false);

  itab1 *= itab2;
  itab1.dump("'itab1' * 'itab2''", false);
  itab1 *= 12;
  itab1.dump("'itab1' * 12", false);
  itab3 = itab1 * itab2;
  itab3.dump("'itab1'' * 'itab2' -> 'itab3'", false);

  itab1 /= itab2;
  itab1.dump("'itab1' / 'itab2''", false);
  itab1 /= 12;
  itab1.dump("'itab1' / 12", false);
  itab3 = itab1 / itab2;
  itab3.dump("'itab1'' / 'itab2' -> 'itab3'", false);

  return 0;
}
