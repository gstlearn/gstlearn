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

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

/**
 * \brief
 * The Style class contains all the recommended styles considerations.
 * It has no special functionalities.<br/>
 * This is only a template class that developers can mimic.
 *
 * This class includes **coding rules** and **doxygen documentation** examples.
 *
 * Since doxygen 1.8.0, Markdown support is available in doxygen comments
 * ([see here](https://www.doxygen.nl/manual/markdown.html)).<br/>
 *
 * Look at *gstlearn* [READMEs](https://github.com/gstlearn/gstlearn)
 * for a quick overview of Markdown syntax.
 *
 * ## Main class features
 *
 * This class demonstrates:
 * - How to use doxygen for documenting the class:
 *    - this global class description
 *    - the private member attributes documentation
 *    - all methods and arguments documentation (at least public ones)
 * - Some mandatory coding rules (see below):
 *    - constructors
 *    - pro active usage of `const`
 *    - error handling with exceptions,
 *    - etc...
 * - How to use some *gstlearn* global features such has:
 *    - Enumerations (see AEnum) - TODO
 *    - `DECLARE_UNUSED` macro (need including geoslib_define.h)
 *    - ...to be expanded
 * - Some coding constraints due to the [customized] SWIG for R version:
 *    - limit the use of function overriding
 *    - do not use namespace or static variables in default argument values
 *    - do not comment the argument name when it is unused (use
 *      `DECLARE_UNUSED`)
 *    - do not use following argument names: 'in', '_*'
 *    - always inherit from pure virtual classes in last position (ex:
 *      ICloneable)
 *
 * ## C++ naming convention
 *
 * Here are some rules about naming C++ types, functions, variables and symbols:
 * - Only usual abbreviations are accepted (ie: min, max, idx, temp, ...)
 * - Constants (macro): all in capitals with underscore (MAX_LENGTH)
 * - Enums: idem
 * - Class or structure names: upper camel case (ie: SparseMatrix, Database)
 * - Aliases (typedef): idem
 * - Class' methods: lower camel case (ie: loadData, getArmonicMean, isEmpty)
 * - Class' attributes (variables): idem
 * - Anything which is *private* or *protected* in a class starts with
 *   underscore (ie: _myMember)
 * - Local variables: lower case separated by underscore (i.e.: temp_value,
 *   color_idx)
 * - Methods visibility;
 *     Method _func() is internal to the class (private or protected)
 *     Method func_() is not exported via SWIG (using #ifndef SWIG)
 * - Use of "InPlace" suffix:
 *     It is used in methods that modify either a member of this class or a
 *     returned argument. Important remark. The member or agument must have been
 *     correctly dimensioned beforehand; no test will be performed (in order to
 *     speed up the processing).
 *
 * ## C++ coding rules
 *
 * Here are some coding rules (good practices) for C++ developers:
 * - Do not use `goto`
 * - Use `break` only in `switch` cases
 * - Limit the number of global variables
 * - Use `const` every where
 * - Do not use old C-style (use C++ STL)
 * - Do not use `using namespace`
 * - Make the includes list as small as possible in C++ header (use forward
 *   declaration)
 * - Expose only what is necessary (all is private by default)
 *
 * More good practices available
 *   [here](www.possibility.com/Cpp/CppCodingStandard.html)
 *
 * Important note: All C++ coding rules are configured in the .clang-format
 *                 file. located in the root folder. Developers can use the
 *                 "Clang format" VS code extension to benefit from the
 *                 automatic formating feature
 *
 * ## Class' constructors/destructors
 *
 * Here are some rules regarding constructors and destructors:
 * - Always add a default constructor (no arguments or arguments with default
 *   values)
 * - Always add a copy constructor (with some exceptions - see Calc* classes)
 * - Always add an assignment operator (with some exceptions - see Calc*
 *   classes)
 * - Always add a virtual destructor
 *
 * ## Documenting Methods
 *
 * Main principles for documenting class methods are the following:
 * - Add one-line comment before each public prototype in the C++ header
 * - Add a doxygen section above each method definition in the C++ body
 * - All this except for trivial methods\tooltip{getters, setters and inline
 *   functions}, constructors and destructors
 * - Define pointer for ENUM. Example to point to EMorpho Enum:  \link
 *   EMorpho.hpp EMorpho \endlink
 */
class GSTLEARN_EXPORT Style
{
public:
  Style();                          // No need for one-line documentation in the header
  Style(const Style& r);            // idem
  Style& operator=(const Style& r); // idem
  virtual ~Style();                 // idem

  // Example of a method with a standard documentation (do not use doxygen
  // here - all the documentation is in the body file!)
  int documentedStandard(int myArg) const;
  // Example of a method with a documentation having a Latex formula (do not use
  // doxygen here!)
  int documentedWithFormula(int myArg) const;
  // Example of a method where argument is not used (do not use doxygen here!)
  int unusedArgument(int a);
  // Just a tentative for formating while typing
  int tentative(int a);

  // Special static function (global) with a default argument (do not use
  // doxygen here!)
  static void myFunction(int myArgInt, double myArgDoubleDef = 2.);

  /**
   * \defgroup Getters Style: Defining the Style
   *
   **/

  /** @addtogroup Getters_0 Accessors to private members value
   * \ingroup Getters
   *
   * This is the demonstration that you can write a unique doxygen documentation
   * for a group of functions. This is useful when all functions handle
   * the same functionality and differ only by their name or arguments.
   *
   * This is also a way to isolate the documentation for trivial functions and
   * make the doxygen page of the class shorter.
   *  @{
   */
  inline double getArgDouble() const { return _argDouble; }
  inline int getArgInt() const { return _argInt; }
  inline const VectorDouble& getArgVectorDouble() const
  {
    return _argVectorDouble;
  }
  inline const VectorInt& getArgVectorInt() const { return _argVectorInt; }
  /**@}*/

  /** @addtogroup Getters_1 Function that update private members value
   * \ingroup Getters
   *
   * Usually, a class should not have these kind of simple accessors (setters).
   * Otherwise, defining member attributes as *private* makes no sense.
   *
   * Mutators\tooltip{non `const` methods that change the state of a class} are
   * usually non trivial and must be defined in the C++ body files.
   * @{
   */
  inline void setArgDouble(double argDouble) { _argDouble = argDouble; }
  inline void setArgInt(int argInt) { _argInt = argInt; }
  inline void setArgVectorDouble(const VectorDouble& argVectorDouble)
  {
    _argVectorDouble = argVectorDouble;
  }
  inline void setArgVectorInt(const VectorInt& argVectorInt)
  {
    _argVectorInt = argVectorInt;
  }
  /**@}*/

private:
  // Example of a private method
  int _increment(int arg) const;

private:
  // Use same line documentation for private members
  int _argInt;                   //!< Private attribute of type `int`
  double _argDouble;             //!< Private attribute of type `double`
  VectorInt _argVectorInt;       //!< Private attribute of type `VectorInt`
  VectorDouble _argVectorDouble; //!< Private attribute of type `VectorDouble`
};
