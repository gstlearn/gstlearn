/******************************************************************************/
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

#include "Basic/VectorNumT.hpp"

/**
 * The Style class is a class which contains all the recommended styles considerations.
 *
 * This includes the documentation. In particular, it demonstrates
 *
 * - the use of a Global class description
 * - the description of each method and its arguments (placed in hpp or cpp section)
 */
class GSTLEARN_EXPORT Style
{
public:
  Style();
  Style(const Style& r);
  Style& operator=(const Style& r);
  virtual ~Style();

  int DocumentedStandard(int myArg) const;
  int DocumentedWithFormula(int myArg) const;
  int UndocumentedArgument(int a);

  static void myFunction(int myArgInt, double myArgDoubleDef = 2.);

private:
  int    _argInt;
  double _argDouble;
  VectorInt _argVectorInt;
  VectorDouble _argVectorDouble;
};

