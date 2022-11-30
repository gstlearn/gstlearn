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

class GSTLEARN_EXPORT Style
{
public:
  Style();
  Style(const Style& r);
  Style& operator=(const Style& r);
  virtual ~Style();

  void DocumentedFunction() const;
  int DocumentedWithFormula(int myArg) const;
  int UndocumentedArgument(int a);

  static void myFunction(int myArgInt, double myArgDoubleDef = 2.);

private:
  int    _argInt;
  double _argDouble;
  VectorInt _argVectorInt;
  VectorDouble _argVectorDouble;
};

