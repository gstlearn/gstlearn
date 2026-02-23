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
 ** This program is meant to test the Template
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

  // ========================
  // Vector of Vector of Real
  // ========================
  mestitle(1, "Case of VectorVectorDouble");
  VectorVectorDouble vrtab1({{1.1, 2.2, 3.3}, {4.4, 5.5, TEST, 7.7}});
  vrtab1.dump("'vrtab1'");

  message("Count of elements = %d\n", vrtab1.count());
  message("Count of elements Defined = %d\n", vrtab1.count(1));
  message("Count of elements Undefined = %d\n", vrtab1.count(-1));
  message("Sum of elements = %lf\n", vrtab1.sum());
  message("Product of elements = %lf\n", vrtab1.prod());
  message("Minimum = %lf\n", vrtab1.minimum());
  message("Maximum = %lf\n", vrtab1.maximum());

  VectorVectorDouble vrtab2({{-2.9, TEST, 3.8}, {4.5, 1.4, 6.3, -5.2}});
  vrtab2.dump("'vrtab2'");
  VectorVectorDouble vrtab3;

  vrtab1.add(vrtab2);
  vrtab1.dump("'vrtab1' + 'vrtab2''");
  vrtab1.addCst(12);
  vrtab1.dump("'vrtab1' + 12");
  vrtab3 = vrtab1.addVec(vrtab2);
  vrtab3.dump("'vrtab1' + 'vrtab2' -> 'vrtab3'");

  vrtab1.subtract(vrtab2);
  vrtab1.dump("'vrtab1' - 'vrtab2''");
  vrtab1.subtractCst(12);
  vrtab1.dump("'vrtab1' - 12");
  vrtab3 = vrtab1.subtractVec(vrtab2);
  vrtab3.dump("'vrtab1' - 'vrtab2' -> 'vrtab3'");

  vrtab1.multiply(vrtab2);
  vrtab1.dump("'vrtab1' * 'vrtab2''");
  vrtab1.multiplyCst(12);
  vrtab1.dump("'vrtab1' * 12");
  vrtab3 = vrtab1.multiplyVec(vrtab2);
  vrtab3.dump("'vrtab1' * 'vrtab2' -> 'vrtab3'");

  vrtab1.divide(vrtab2);
  vrtab1.dump("'vrtab1' / 'vrtab2''");
  vrtab1.divideCst(12);
  vrtab1.dump("'vrtab1' / 12");
  vrtab3 = vrtab1.divideVec(vrtab2);
  vrtab3.dump("'vrtab1' / 'vrtab2' -> 'vrtab3'");

  // ===========================
  // Vector of Vector of Integer
  // ===========================
  mestitle(1, "Case of VectorVectorInt");
  VectorVectorInt vitab1({{1, 2, 3}, {4, 5, ITEST, 7}});
  vitab1.dump("'vitab1'");

  message("Count of elements = %d\n", vitab1.count());
  message("Count of elements Defined = %d\n", vitab1.count(1));
  message("Count of elements Undefined = %d\n", vitab1.count(-1));
  message("Sum of elements = %lf\n", vitab1.sum());
  message("Product of elements = %lf\n", vitab1.prod());
  message("Minimum = %lf\n", vitab1.minimum());
  message("Maximum = %lf\n", vitab1.maximum());

  VectorVectorInt vitab2({{-2, ITEST, 3}, {4, 1, 6, -5}});
  vitab2.dump("'vitab2'");
  VectorVectorInt vitab3;

  vitab1.add(vitab2);
  vitab1.dump("'vitab1' + 'vitab2''");
  vitab1.addCst(12);
  vitab1.dump("'vitab1' + 12");
  vitab3 = vitab1.addVec(vitab2);
  vitab3.dump("'vitab1' + 'vitab2' -> 'vitab3'");

  vitab1.subtract(vitab2);
  vitab1.dump("'vitab1' - 'vitab2''");
  vitab1.subtractCst(12);
  vitab1.dump("'vitab1' - 12");
  vitab3 = vitab1.subtractVec(vitab2);
  vitab3.dump("'vitab1' - 'vitab2' -> 'vitab3'");

  vitab1.multiply(vitab2);
  vitab1.dump("'vitab1' * 'vitab2''");
  vitab1.multiplyCst(12);
  vitab1.dump("'vitab1' * 12");
  vitab3 = vitab1.multiplyVec(vitab2);
  vitab3.dump("'vitab1' * 'vitab2' -> 'vitab3'");

  vitab1.divide(vitab2);
  vitab1.dump("'vitab1' / 'vitab2''");
  vitab1.divideCst(12);
  vitab1.dump("'vitab1' / 12");
  vitab3 = vitab1.divideVec(vitab2);
  vitab3.dump("'vitab1' / 'vitab2' -> 'vitab3'");

  return 0;
}
