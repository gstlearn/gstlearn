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
#pragma once

#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT VectorHelper
{
public:
  static VectorInt initVInt(Id nval, Id value = 0.);
  static VectorDouble initVDouble(Id nval, double value = 0.);
  static VectorVectorDouble initVVDouble(Id nval1, Id nval2, double value = 0.);
  static VectorVectorInt initVVInt(Id nval1, Id nval2, Id value = 0);

  static VectorInt initVInt(const Id* values, Id number);
  static VectorInt initVInt(const I32* values, Id number);
  static VectorDouble initVDouble(const double* values, Id number);
  static VectorVectorDouble initVVDouble(const double* value, Id n1, Id n2);

  static VectorString initVString(Id ntab, char** names);

#ifndef SWIG
  static void dumpStats(const String& title, constvect vect, Id nmax = -1);
  static void dumpRange(const String& title, constvect vect, Id nmax = -1);
#endif
  static void dumpStats(const String& title, const VectorDouble& vectin, Id nmax = -1);
  static void dumpRange(const String& title, const VectorDouble& vectin, Id nmax = -1);
  static void dumpRange(const String& title, const VectorInt& vect);
  static void dumpNNZ(const String& title, const VectorDouble& vect, Id nclass = 10);

  static void capInPlace(VectorDouble& vec, double vmin = TEST, double vmax = TEST);
  static void capInPlaceVVD(VectorVectorDouble& vec, double vmin = TEST, double vmax = TEST);
  static double extensionDiagonal(const VectorDouble& mini, const VectorDouble& maxi);

  static Id count(const VectorVectorInt& vec);
  static double cumulLog(const VectorDouble& vec);
  static VectorInt cumulIncrement(const VectorVectorInt& vec);
  static VectorDouble quantiles(const VectorDouble& vec, const VectorDouble& probas);

  static bool isEqual(const VectorDouble& v1,
                      const VectorDouble& v2,
                      double eps = EPSILON10);
  static bool isEqualExtended(const VectorDouble& v1,
                              const VectorDouble& v2,
                              double eps           = EPSILON10,
                              bool flagRelative    = true,
                              bool flagAbsolute    = false,
                              const String& string = "");

  static void sequenceInPlace(Id n, VectorInt& vec);
  static VectorInt sequence(Id number, Id ideb = 0, Id step = 1);
  static VectorDouble sequenceVD(double valFrom,
                                 double valTo,
                                 double valStep = 1.,
                                 double ratio   = 1.);
  static void fillUndef(VectorDouble& vec, double repl);

#ifndef SWIG
  static void addMultiplyConstantInPlace(double val1,
                                         const constvect in,
                                         vect out,
                                         Id iad);
  static double innerProductCV(const constvect veca, const constvect vecb);
  static void addMultiplyVectVectInPlace(const constvect in1,
                                         const constvect in2,
                                         vect out,
                                         Id iad);

  static void addInPlace(constvect in, vect dest);
#endif
  static void addInPlace(const double* veca,
                         const double* vecb,
                         double* res,
                         Id size);
  static void addSquareInPlace(VectorDouble& dest, const VectorDouble& src);

#ifndef SWIG
  static VectorDouble subtract(constvect veca, constvect vecb);
  static void subtractInPlace(const constvect in1, const constvect in2, vect outv);
#endif
  static void multiplyComplexInPlace(const VectorDouble& vecaRe,
                                     const VectorDouble& vecaIm,
                                     const VectorDouble& vecbRe,
                                     const VectorDouble& vecbIm,
                                     VectorDouble& resRe,
                                     VectorDouble& resIm);

  static void addMultiplyConstantInPlace(double val1,
                                         const VectorDouble& in1,
                                         VectorDouble& outv,
                                         Id iad);
  static void addMultiplyConstantInPlace(double val1,
                                         const VectorVectorDouble& in1,
                                         VectorVectorDouble& outv);
  static void copy(const VectorDouble& vecin, VectorDouble& vecout, Id size = -1);
  static void copy(const VectorInt& vecin, VectorInt& vecout, Id size = -1);
  static void copy(const VectorVectorDouble& inv, VectorVectorDouble& outv);

  static void mean1AndMean2ToStdev(const VectorDouble& mean1,
                                   const VectorDouble& mean2,
                                   VectorDouble& std,
                                   Id number);

  static void normalize(double* tab, Id ntab);
  static void normalizeFromGaussianDistribution(VectorDouble& vec,
                                                double mini = 0.,
                                                double maxi = 1.);
  static VectorDouble normalScore(const VectorDouble& data,
                                  const VectorDouble& wt = VectorDouble());
  static VectorDouble qnormVec(const VectorDouble& vec);
  static VectorDouble pnormVec(const VectorDouble& vec);
  static VectorDouble concatenate(const VectorDouble& veca, const VectorDouble& vecb);
  static void concatenateInPlace(VectorDouble& veca, const VectorDouble& vecb);
  static VectorDouble power(const VectorDouble& vec, double power);
  static VectorDouble inverse(const VectorDouble& vec);
#ifndef SWIG
  static void power(VectorDouble& res, const constvect vec, double power);
  static void inverse(VectorDouble& res, const constvect vec);
#endif
  // Next line is kept only fr VVD compatibility
  static double innerProduct(const double* veca, const double* vecb, Id size = -1);

  static VectorDouble crossProduct3D(const VectorDouble& veca, const VectorDouble& vecb);
  static void crossProduct3DInPlace(const double* a, const double* b, double* v);

  static VectorDouble cumsum(const VectorDouble& vecin, bool flagAddZero, bool revert = false);
  static void cumulateInPlace(VectorDouble& vec);
  static void cumulate(VectorDouble& veca, const VectorDouble& vecb, double coeff = 1., double addval = 0.);
  static void getMostSignificant(const VectorDouble& vec, double tol = EPSILON6, Id nmax = -1);

  static VectorDouble simulateUniform(Id n        = 1,
                                      double mini = 0.,
                                      double maxi = 1.);
  static VectorDouble simulateBernoulli(Id n         = 1,
                                        double proba = 0.5,
                                        double vone  = 1.,
                                        double velse = 0.);
  static VectorDouble simulateGaussian(Id n         = 1,
                                       double mean  = 0.,
                                       double sigma = 1.);
  static void simulateGaussianInPlace(VectorDouble& vec,
                                      double mean  = 0.,
                                      double sigma = 1.);
  static VectorInt sampleRanks(Id ntotal,
                               double proportion = 0.,
                               Id number         = 0,
                               Id seed           = 242141,
                               Id optSort        = 0);
  static void normalizeCodir(Id ndim, VectorDouble& codir);

  static bool isInList(const VectorInt& vec, Id item);
  static VectorInt sort(const VectorInt& vecin, bool ascending = true, Id size = -1);
  static VectorDouble sort(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static void sortInPlace(VectorInt& vecin, bool ascending = true, Id size = -1);
  static void sortInPlace(VectorDouble& vecin, bool ascending = true, Id size = -1);
  static bool isSorted(const VectorDouble& vec, bool ascending = true);
  static VectorDouble unique(const VectorDouble& vecin, Id size = -1);
  static VectorInt unique(const VectorInt& vecin, Id size = -1);
  static VectorInt orderRanks(const VectorInt& vecin, bool ascending = true, Id size = -1);
  static VectorInt orderRanks(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static VectorInt sortRanks(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static VectorInt reorder(const VectorInt& vecin, const VectorInt& order, Id size = -1);
  static VectorDouble reorder(const VectorDouble& vecin, const VectorInt& order, Id size = -1);
  static VectorDouble revert(const VectorDouble& vecin);
  static VectorInt revert(const VectorInt& vecin);
  static VectorDouble sample(const VectorDouble& vecin, const VectorInt& indKeep);

  static void arrangeInPlace(Id safe,
                             VectorInt& ranks,
                             VectorDouble& values,
                             bool ascending = true,
                             Id size        = -1);
  static void arrangeInPlace(Id safe,
                             VectorInt& ranks,
                             VectorInt& values,
                             bool ascending = true,
                             Id size        = -1);
  static VectorInt filter(const VectorInt& vecin,
                          Id vmin        = ITEST,
                          Id vmax        = ITEST,
                          bool ascending = true);
  static VectorInt complement(const VectorInt& vec, const VectorInt& sel);

  static std::pair<double, double> rangeVals(const VectorDouble& vec);
  static void unflattenInPlace(const VectorDouble& vd, VectorVectorDouble& vvd);
  static void flattenInPlace(const VectorVectorDouble& vvd, VectorDouble& vd);
  static VectorDouble flatten(const VectorVectorDouble& vvd);
  static VectorVectorDouble unflatten(const VectorDouble& vd, const VectorInt& sizes);
  static void linearCombinationInPlace(double val1,
                                       const VectorDouble& vd1,
                                       double val2,
                                       const VectorDouble& vd2,
                                       VectorDouble& outv);
  static void linearCombinationVVDInPlace(double val1,
                                          const VectorVectorDouble& vvd1,
                                          double val2,
                                          const VectorVectorDouble& vvd2,
                                          VectorVectorDouble& outv);

  static VectorDouble suppressTest(const VectorDouble& vecin);
  static void extractInPlace(const VectorDouble& vecin, VectorDouble& vecout, Id start);
  static void mergeInPlace(const VectorDouble& vecin, VectorDouble& vecout, Id start);

  static void transformVD(VectorDouble& tab, Id oper_choice = 1);

  static void squeezeAndStretchInPlaceForward(const VectorDouble& vecin,
                                              VectorDouble& vecout,
                                              double origin,
                                              double mesh,
                                              double top,
                                              double bot);
  static void squeezeAndStretchInPlaceBackward(const VectorDouble& vecin,
                                               VectorDouble& vecout,
                                               double origin,
                                               double mesh,
                                               double top,
                                               double bot);

  static Id whereMinimum(const VectorDouble& tab);
  static Id whereMaximum(const VectorDouble& tab);
  static Id whereElement(const VectorInt& tab, Id target);
  static Id whereElement(const VectorInt& tab, Id target, Id start);
  static bool isIsotropic(const VectorVectorInt& sampleRanks);

  static VectorDouble reduceOne(const VectorDouble& vecin, Id index);
  static VectorDouble reduce(const VectorDouble& vecin, const VectorInt& vindex);
  static VectorDouble compress(const VectorDouble& vecin, const VectorInt& vindex);
  static void truncateDecimalsInPlace(VectorDouble& vec, Id ndec);
  static void truncateDigitsInPlace(VectorDouble& vec, Id ndec);

  /**
   * \defgroup Operators: List of basic operators bewteen numerical Vectors and scalars
   *
   **/

  /** @addtogroup Operate_1
   * \ingroup Operators
   *
   * Action: These operators perform basic operations between two vectors of the same size
   * (or a vector and a scalar) and return a new vector as result (not in place)
   *  @{
   */
  static VectorDouble add(const VectorDouble& v1, const VectorDouble& v2);
  static VectorDouble addCst(const VectorDouble& v1, double v2);

  static VectorDouble subtract(const VectorDouble& v1, const VectorDouble& v2);
  static VectorDouble subtractCst(const VectorDouble& v1, double v2, bool flagOpposite = false);

  static VectorDouble multiply(const VectorDouble& v1, const VectorDouble& v2);
  static VectorDouble multiplyCst(const VectorDouble& v1, double v2);

  static VectorDouble divide(const VectorDouble& v1, const VectorDouble& v2);
  static VectorDouble divideCst(const VectorDouble& v1, double v2, bool flagOpposite = false);
  /**@}*/

  /** @addtogroup Operate_2
   * \ingroup Operators
   *
   * Action: These operators perform basic operations between two vectors of the same size
   * (or a vector and a scalar) and return the result in the first argument.
   * They are not meant to be exposted in foreign language and, therefore, they benefit from
   * the polymorphism.
   *
   * They have been defined to allow efficient calculations without using any allocation
   *  @{
   */
#ifndef SWIG
  static void add(VectorDouble& vecout, const VectorDouble& v1, const VectorDouble& v2);
  static void increment(VectorDouble& vecout, const VectorDouble& v1);
  static void add(VectorDouble& vecout, const VectorDouble& v1, double v2);
  static void increment(VectorDouble& vecout, const VectorDouble& v1, double v2);

  static void subtract(VectorDouble& vecout, const VectorDouble& v1, const VectorDouble& v2);
  static void decrement(VectorDouble& vecout, const VectorDouble& v1);
  static void subtract(VectorDouble& vecout, const VectorDouble& v1, double v2, bool flagOpposite = false);
  static void decrement(VectorDouble& vecout, const VectorDouble& v1, double v2);

  static void multiply(VectorDouble& vecout, const VectorDouble& v1, const VectorDouble& v2);
  static void multiplyAssign(VectorDouble& vecout, const VectorDouble& v1);
  static void multiply(VectorDouble& vecout, const VectorDouble& v1, double v2);
  static void multiplyAssign(VectorDouble& vecout, const VectorDouble& v1, double v2);

  static void divide(VectorDouble& vecout, const VectorDouble& v1, const VectorDouble& v2);
  static void divideAssign(VectorDouble& vecout, const VectorDouble& v1);
  static void divide(VectorDouble& vecout, const VectorDouble& v1, double v2, bool flagOpposite = false);
  static void divideAssign(VectorDouble& vecout, const VectorDouble& v1, double v2);
#endif
  /**@}*/
};

// typedef VectorHelper VH;
class VH: public VectorHelper
{
};

/**
 * \defgroup Symbolic: List of basic operators bewteen numerical Vectors and scalars
 * (expressed using standad operators)
 *
 **/

/** @addtogroup Operate_3
 * \ingroup Symbolic
 *
 * Action: These operators perform basic operations between two vectors of the same size
 * or between a vector and a scalar.
 * They are presented as an overload of basic operators ''+', '-', '*', and '/'.
 * They are not meant to be exposted in foreign language.
 *  @{
 */

/**
 * @brief Operator attached to "+" symbol
 */
inline VectorDouble operator+(const VectorDouble& v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::add(result, v1, v2);
  return result;
}
inline VectorDouble operator+(const VectorDouble& v1, double v2)
{
  VectorDouble result;
  VectorHelper::add(result, v1, v2);
  return result;
}
inline VectorDouble operator+(double v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::add(result, v2, v1);
  return result;
}
inline VectorDouble& operator+=(VectorDouble& v1, const VectorDouble& v2)
{
  VectorHelper::add(v1, v1, v2);
  return v1;
}
inline VectorDouble& operator+=(VectorDouble& v1, double v2)
{
  VectorHelper::add(v1, v1, v2);
  return v1;
}

/**
 * @brief Operator attached to "-" symbol
 */
inline VectorDouble operator-(const VectorDouble& v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::subtract(result, v1, v2);
  return result;
}
inline VectorDouble operator-(const VectorDouble& v1, double v2)
{
  VectorDouble result;
  VectorHelper::subtract(result, v1, v2, false);
  return result;
}
inline VectorDouble operator-(double v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::subtract(result, v2, v1, true);
  return result;
}
inline VectorDouble& operator-=(VectorDouble& v1, const VectorDouble& v2)
{
  VectorHelper::subtract(v1, v1, v2);
  return v1;
}
inline VectorDouble& operator-=(VectorDouble& v1, double v2)
{
  VectorHelper::subtract(v1, v1, v2, false);
  return v1;
}

/**
 * @brief Operator attached to "*" symbol
 */
inline VectorDouble operator*(const VectorDouble& v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::multiply(result, v1, v2);
  return result;
}
inline VectorDouble operator*(const VectorDouble& v1, double v2)
{
  VectorDouble result;
  VectorHelper::multiply(result, v1, v2);
  return result;
}
inline VectorDouble operator*(double v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::multiply(result, v2, v1);
  return result;
}
inline VectorDouble& operator*=(VectorDouble& v1, const VectorDouble& v2)
{
  VectorHelper::multiply(v1, v1, v2);
  return v1;
}
inline VectorDouble& operator*=(VectorDouble& v1, double v2)
{
  VectorHelper::multiply(v1, v1, v2);
  return v1;
}

/**
 * @brief Operator attached to "/" symbol
 */
inline VectorDouble operator/(const VectorDouble& v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::divide(result, v1, v2);
  return result;
}
inline VectorDouble operator/(const VectorDouble& v1, double v2)
{
  VectorDouble result;
  VectorHelper::divide(result, v1, v2, false);
  return result;
}
inline VectorDouble operator/(double v1, const VectorDouble& v2)
{
  VectorDouble result;
  VectorHelper::divide(result, v2, v1, true);
  return result;
}
inline VectorDouble& operator/=(VectorDouble& v1, const VectorDouble& v2)
{
  VectorHelper::divide(v1, v1, v2);
  return v1;
}
inline VectorDouble& operator/=(VectorDouble& v1, double v2)
{
  VectorHelper::divide(v1, v1, v2, false);
  return v1;
}

} // namespace gstlrn
