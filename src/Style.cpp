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
#include "Style.hpp"

#include "Basic/AStringable.hpp"

Style::Style()
    : _argInt(0),
      _argDouble(0.),
      _argVectorInt(),
      _argVectorDouble()
{
}

Style::Style(const Style& r)
    : _argInt(r._argInt),
      _argDouble(r._argDouble),
      _argVectorInt(r._argVectorInt),
      _argVectorDouble(r._argVectorDouble)
{
}

Style& Style::operator=(const Style &r)
{
  if (this != &r)
  {
    _argInt = r._argInt;
    _argDouble = r._argDouble;
    _argVectorInt = r._argVectorInt;
    _argVectorDouble = r._argVectorDouble;
  }
  return *this;
}

Style::~Style()
{
}

/**
 * A method with standard argument documentation
 * @param myArg Here should be placed the description of this argument
 * @return Description of the returned argument
 */
int Style::DocumentedStandard(int myArg) const
{
  message("Documented Function\n");
  message("Value of MyArg = %d\n",myArg);
  return 0;
}

/**
 * Documentation with Formula
 *
 * \f$ p1=(x_1,y_1) \f$ and \f$ p2=(x_2,y_2) \f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$
 *
 * @param myArg Here should be placed the description of this argument
 * @return
 */
int Style::DocumentedWithFormula(int myArg) const
{
  message("Input Argument = %d\n", myArg);
  return myArg;
}

/**
 * Description of a static function
 * @param myArgInt Integer argument
 * @param myArgDoubleDef Double argument
 */
void Style::myFunction(int myArgInt, double myArgDoubleDef)
{
  message("Input Argument (int) = %d\n", myArgInt);
  message("Input Argument (double with default value) = %lf\n",myArgDoubleDef);
}

/**
 * Undocumented argument
 * @param a Input argument not used (but documented)
 * @return Error returned code
 */
int Style::UndocumentedArgument(int a)
{
  SYMBOL_UNUSED(a);
  return 0;
}
