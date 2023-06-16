/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT VectorHelper
{
public:
  static VectorInt          initVInt(int nval, int value = 0.);
  static VectorDouble       initVDouble(int nval, double value = 0.);
  static VectorVectorDouble initVVDouble(int nval1, int nval2, double value = 0.);
  static VectorVectorInt    initVVInt(int nval1, int nval2, int value = 0);

  static VectorInt          initVInt(const int* values, int number);
  static VectorDouble       initVDouble(const double* values, int number);
  static VectorVectorDouble initVVDouble(const double* value, int n1, int n2);

  static void display(const String &title, const VectorDouble &vect); // TODO rename
  static void display(const String &title, const VectorVectorDouble &vect);
  static void display(const String &title, const VectorString &vect);
  static void display(const String &title, const VectorInt &vect);

  static String toString(const VectorDouble& vec); // TODO rename
  static String toString(const VectorVectorDouble& vec);
  static String toString(const VectorString& vec);
  static String toString(const VectorInt& vec);

  static void displayStats(const String &title, const VectorDouble &vect);
  static void displayRange(const String &title, const VectorDouble &vect);
  static void displayRange(const String &title, const VectorInt &vect);
  static void displayNNZ(const String &title, const VectorDouble &vect, int nclass = 10);

  static int maximum(const VectorInt &vec, bool flagAbs = false);
  static int minimum(const VectorInt &vec, bool flagAbs = false);
  static double maximum(const VectorDouble &vec, bool flagAbs = false);
  static double minimum(const VectorDouble &vec, bool flagAbs = false);
  static double maximum(const VectorVectorDouble &vec, bool flagAbs = false);
  static double minimum(const VectorVectorDouble &vec, bool flagAbs = false);
  static int product(const VectorInt& nx);
  static double product(const VectorDouble& vec);
  static int countUndefined(const VectorDouble& vec);
  static int countDefined(const VectorDouble& vec);
  static double extensionDiagonal(const VectorDouble& mini, const VectorDouble& maxi);

  static int    cumul(const VectorInt& vec);
  static double cumul(const VectorDouble &vec);
  static double mean(const VectorDouble &vec);
  static double variance(const VectorDouble &vec);
  static double stdv(const VectorDouble &vec);
  static double norm(const VectorDouble &vec);
  static double normOptim(const VectorDouble&vec);

  static double correlation(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble quantiles(const VectorDouble& vec,
                                const VectorDouble& probas);

  static bool isConstant(const VectorDouble& vect, double refval = TEST);
  static bool isConstant(const VectorInt& vect, int refval = ITEST);
  static bool isSame(const VectorDouble &v1,
                     const VectorDouble &v2,
                     double eps = EPSILON10);
  static bool isSame(const VectorInt &v1, const VectorInt &v2);

  static VectorInt sequence(int number, int ideb = 0);
  static VectorDouble sequence(double valFrom,
                               double valTo,
                               double valStep = 1.,
                               double ratio = 1.);
  static void fill(VectorDouble& vec, double v, int size = 0);
  static void fill(VectorInt& vec, int v, int size = 0);
  static void fill(VectorVectorDouble &vec, double value);
  static void fillUndef(VectorDouble& vec, double repl);

  static VectorDouble add(const VectorDouble &veca, const VectorDouble &vecb);
  static void addInPlace(VectorDouble &dest, const VectorDouble &src);
  static void addInPlace(const VectorDouble &veca,
                         const VectorDouble &vecb,
                         VectorDouble &res);
  static VectorDouble subtract(const VectorDouble& veca, const VectorDouble& vecb);
  static void subtractInPlace(VectorDouble &dest, const VectorDouble &src);

  static void multiplyInPlace(VectorDouble& vec, const VectorDouble& v);
  static void divideInPlace(VectorDouble& vec, const VectorDouble& v);

  static void multiplyConstant(VectorDouble& vec, double v);
  static void divideConstant(VectorDouble& vec, double v);
  static void copy(VectorDouble& veca, const VectorDouble& vecb);
  static void addConstant(VectorDouble& vec, double v);
  static void addConstant(VectorInt& vec, int v);

  static void normalize(VectorDouble& vec);
  static void normalizeFromGaussianDistribution(VectorDouble &vec,
                                                double mini = 0.,
                                                double maxi = 1.);
  static VectorDouble normalScore(const VectorDouble& data,
                                  const VectorDouble& wt = VectorDouble());
  static VectorDouble qnormVec(const VectorDouble& vec);
  static VectorDouble pnormVec(const VectorDouble& vec);
  static VectorDouble concatenate(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble power(const VectorDouble& vec, double power);
  static VectorDouble inverse(const VectorDouble& vec);

  static double innerProduct(const VectorDouble &veca, const VectorDouble &vecb);
  static double innerProductSubVec(const VectorDouble &veca,
                                   const VectorDouble &vecb,
                                   int size);
  static VectorDouble crossProduct(const VectorDouble &veca, const VectorDouble &vecb);

  static void cumulate(VectorDouble &veca,
                       const VectorDouble &vecb,
                       double coeff = 1.,
                       double addval = 0.);

  static VectorDouble simulateUniform(int n = 1,
                                      double mini = 0.,
                                      double maxi = 1.);
  static VectorDouble simulateBernoulli(int n = 1,
                                        double proba = 0.5,
                                        double vone = 1.,
                                        double velse = 0.);
  static VectorDouble simulateGaussian(int n = 1,
                                       double mean = 0.,
                                       double sigma = 1.);
  static void simulateGaussianInPlace(VectorDouble &vect,
                                      double mean = 0.,
                                      double sigma = 1.);
  static VectorInt sampleRanks(int ntotal,
                               double proportion = 0.,
                               int number = 0,
                               int seed = 242141,
                               int optSort = 0);

  static VectorInt    sort(const VectorInt& vecin, bool ascending = true);
  static VectorDouble sort(const VectorDouble& vecin, bool ascending = true);
  static VectorInt    orderRanks(const VectorDouble& vecin, bool ascending = true);
  static VectorInt    sortRanks(const VectorDouble& vecin, bool ascending = true);
  static VectorDouble unique(const VectorDouble& vecin);
  static VectorInt    unique(const VectorInt& vecin);
  static VectorInt filter(const VectorInt &vecin,
                          int vmin = ITEST,
                          int vmax = ITEST,
                          bool ascending = true);

  static std::pair<double,double> rangeVals(const VectorDouble& vec);

  static VectorDouble flatten(const VectorVectorDouble& vvd);
};

//typedef VectorHelper VH;
class VH: public VectorHelper {};

