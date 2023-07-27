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

  static void display(const String &title, const VectorDouble &vect, bool skipLine = true); // TODO rename
  static void display(const String &title, const VectorVectorDouble &vect, bool skipLine = true);
  static void display(const String &title, const VectorString &vect, bool skipLine = true);
  static void display(const String &title, const VectorInt &vect, bool skipLine = true);

  static String toStringAsVD(const VectorDouble& vec); // TODO rename
  static String toStringAsVVD(const VectorVectorDouble& vec);
  static String toStringAsVS(const VectorString& vec);
  static String toStringAsVI(const VectorInt& vec);

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
  static double normDistance(const VectorDouble& veca, const VectorDouble& vecb);
  static double correlation(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble quantiles(const VectorDouble& vec,
                                const VectorDouble& probas);

  static bool isConstant(const VectorDouble& vect, double refval = TEST);
  static bool isConstant(const VectorInt& vect, int refval = ITEST);
  static bool isSame(const VectorDouble &v1,
                     const VectorDouble &v2,
                     double eps = EPSILON10);
  static bool isSame(const VectorInt &v1, const VectorInt &v2);

  static VectorInt sequence(int number, int ideb = 0, int step = 1);
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
                         VectorDouble &res,
                         int size = 0);
  static void addInPlace(const double *veca,
                         const double *vecb,
                         double *res,
                         int size);
  static VectorDouble subtract(const VectorDouble& veca, const VectorDouble& vecb);
  static void subtractInPlace(VectorDouble &dest, const VectorDouble &src);
  static void subtractInPlace(VectorInt &dest, const VectorInt &src);

  static void multiplyInPlace(VectorDouble& vec, const VectorDouble& v);
  static void divideInPlace(VectorDouble& vec, const VectorDouble& v);

  static void multiplyConstant(VectorDouble& vec, double v);
  static void multiplyConstantInPlace(const VectorDouble& vec, double v, VectorDouble& vecout);
  static void divideConstant(VectorDouble& vec, double v);
  static void copy(const VectorDouble& vecin, VectorDouble& vecout, int size = -1);
  static void copy(const VectorInt &vecin, VectorInt &vecout, int size = -1);
  static void addConstant(VectorDouble& vec, double v);
  static void addConstant(VectorInt& vec, int v);

  static void normalize(VectorDouble& vec);
  static void normalize(double *tab, int ntab);
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

  static double innerProduct(const VectorDouble &veca, const VectorDouble &vecb, int size = -1);
  static double innerProduct(const double* veca, const double* vecb, int size);

  static VectorDouble crossProduct3D(const VectorDouble &veca, const VectorDouble &vecb);
  static void crossProduct3DInPlace(const double *a, const double *b, double *v);

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

  static VectorInt    sort(const VectorInt& vecin, bool ascending = true, int size = -1);
  static VectorDouble sort(const VectorDouble& vecin, bool ascending = true, int size = -1);
  static void         sortInPlace(VectorInt& vecin, bool ascending = true, int size = -1);
  static void         sortInPlace(VectorDouble& vecin, bool ascending = true, int size = -1);
  static VectorDouble unique(const VectorDouble& vecin, int size = -1);
  static VectorInt    unique(const VectorInt& vecin, int size = -1);
  static VectorInt    orderRanks(const VectorInt& vecin, bool ascending = true, int size = -1);
  static VectorInt    orderRanks(const VectorDouble& vecin, bool ascending = true, int size = -1);
  static VectorInt    sortRanks(const VectorDouble& vecin, bool ascending = true, int size = -1);
  static VectorInt    reorder(const VectorInt& vecin, const VectorInt& order, int size = -1);
  static VectorDouble reorder(const VectorDouble& vecin, const VectorInt& order, int size = -1);
  static void arrangeInPlace(int safe,
                             VectorInt &ranks,
                             VectorDouble &values,
                             bool ascending = true,
                             int size = -1);
  static void arrangeInPlace(int safe,
                             VectorInt &ranks,
                             VectorInt &values,
                             bool ascending = true,
                             int size = -1);
  static VectorInt filter(const VectorInt &vecin,
                          int vmin = ITEST,
                          int vmax = ITEST,
                          bool ascending = true);

  static std::pair<double,double> rangeVals(const VectorDouble& vec);

  static VectorDouble flatten(const VectorVectorDouble& vvd);
};

//typedef VectorHelper VH;
class VH: public VectorHelper {};

