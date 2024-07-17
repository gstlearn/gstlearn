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
// Use double quote "" for include files from our source code
#include "API/Style.hpp"
#include "Basic/AStringable.hpp"
// Use corners <> for other include files
#include <iostream>

/**
 * Default constructor
 *
 * Always use the initialization list to set the default value of each class
 * member. Put attributes in the same order than their declaration in the C++
 * header.
 */
Style::Style()
  : _argInt(0)
  , _argDouble(0.)
  , _argVectorInt()
  , _argVectorDouble()
{
}

/**
 * Copy constructor
 *
 * Always use the initialization list to set the value of each class member.
 * Put attributes in the same order than their declaration in the C++ header.
 */
Style::Style(const Style& r)
  : _argInt(r._argInt)
  , _argDouble(r._argDouble)
  , _argVectorInt(r._argVectorInt)
  , _argVectorDouble(r._argVectorDouble)
{
}

/**
 * Assignment operator
 *
 * Always protect from self copy with the if statement
 */
Style& Style::operator=(const Style& r)
{
  if (this != &r)
  {
    _argInt          = r._argInt;
    _argDouble       = r._argDouble;
    _argVectorInt    = r._argVectorInt;
    _argVectorDouble = r._argVectorDouble;
  }
  return *this;
}

/**
 * Destructor
 *
 * Always free pointers and clear the lists members of the class
 */
Style::~Style() {}

/**
 * A method with standard argument documentation
 *
 * \param[in] myArg Here should be placed the description of this argument
 *
 * \return Description of the returned value
 */
int Style::documentedStandard(int myArg) const
{
  std::cout << "Documented Function" << std::endl;
  message("Value of MyArg = %d\n", myArg);
  
  return 0;
}

/**
 * Documentation with Latex formula
 *
 * The distance between \f$ p1=(x_1,y_1) \f$ and \f$ p2=(x_2,y_2) \f$ is
 * \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$ (this formula may need to do this:
 * https://github.com/doxygen/doxygen/issues/7484#issuecomment-572503569)
 *
 * \see SpaceRN::getDistance
 *
 * \param[in] myArg Here should be placed the description of this argument
 *
 * \return The value of the argument + 1
 */
int Style::documentedWithFormula(int myArg) const
{
  message("Input Argument = %d\n", myArg);
  return _increment(myArg);
}

/**
 * Description of a static function
 *
 * \param[in] myArgInt Integer argument
 * \param[in] myArgDoubleDef Double argument
 */
void Style::myFunction(int myArgInt, double myArgDoubleDef)
{
  message("Input Argument (int) = %d\n", myArgInt);
  message("Input Argument (double with default value) = %lf\n", myArgDoubleDef);
}

/**
 * A function where the argument is not used (could be the case in abstract
 * methods)
 *
 * \param[in] a Input argument not used (but documented)
 *
 * \return Error returned code
 */
int Style::unusedArgument(int a)
{
  // Use the DECLARE_UNUSED macro to prevent compiler warning
  // Do not comment the argument (NOT /*a*/)
  DECLARE_UNUSED(a);
  return 0;
}

/**
 *  Example of a private method
 *
 *  \param[in] arg Input integer argument
 *
 *  \return The input argument incremented
 */
int Style::_increment(int arg) const
{
  arg++;
  return arg;
}
