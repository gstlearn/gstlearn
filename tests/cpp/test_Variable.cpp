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

#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "DataBase/AVariable.hpp"
#include "DataBase/VariableBool.hpp"
#include "DataBase/VariableDouble.hpp"
#include "DataBase/VariableFormat.hpp"
#include "DataBase/VariableInt.hpp"
#include "DataBase/VariableString.hpp"
#include <iostream>

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_Variable-");

  // ================================
  // Testing the VariableDouble class
  // ================================

  // Creating a VariableDouble called "test_double"
  mestitle(1, "Testing class VariableDouble");
  VariableDouble varDbl("test_double");

  // Filling it with a four values (including a NA) and displaying them (old style)
  auto NADbl = getNA<double>();
  varDbl.setValues({1., 3., NADbl, 4.});
  varDbl.display_old();

  // Checking the number of samples
  std::cout << "Number of values = " << varDbl.getNSamples() << std::endl;
  std::cout << "Number of defined values = " << varDbl.getNAllDefined()
            << std::endl;

  // Resizing ti to 7 values, filling it with 5. and displaying it again
  varDbl.resize(7, 5.);
  varDbl.getValues().display();

  // Filling the variable with 1.4 and displaying it again
  varDbl.fill(1.4);
  varDbl.getValues().display();

  // Inquiring the value of the 4th element, changing it to 3.6 and displaying the variable again
  auto iech = 3;
  auto value = varDbl.getValueAsType(iech);
  std::cout << "Value of sample #" << iech << " = " << value << std::endl;
  varDbl.setValue(3.6, iech);
  varDbl.getValues().display();

  // =============================
  // Testing the VariableInt class
  // =============================

  // Creating a VariableInt called "test_int"
  mestitle(1, "Testing class VariableInt");
  VariableInt varInt("test_int");

  // Filling it with a four values (including a NA) and displaying them (old style)
  auto NAInt = getNA<Id>();
  varInt.setValues({1, 3, NAInt, 4});
  varInt.display_old();

  // Checking the number of samples
  std::cout << "Number of values = " << varInt.getNSamples() << std::endl;
  std::cout << "Number of defined values = " << varInt.getNAllDefined()
            << std::endl;

  // ==============================
  // Testing the VariableBool class
  // ==============================

  // Creating a VariableBool called "test_bool"
  mestitle(1, "Testing class VariableBool");
  VariableBool varBool("test_bool");

  // Filling it with a four values (including a NA) and displaying them (old style)
  varBool.setValues({true, false, true, true});
  varBool.display_old();

  // Checking the number of samples
  std::cout << "Number of values = " << varBool.getNSamples() << std::endl;
  std::cout << "Number of true values = " << varBool.getNTrue() << std::endl;
  std::cout << "Number of false values = " << varBool.getNFalse() << std::endl;

  // ================================
  // Testing the VariableString class
  // ================================

  // Creating a VariableString called "test_string"
  mestitle(1, "Testing class VariableString");
  VariableString varString("test_string");

  // Filling it with a four values (including a NA) and displaying them (old style)
  auto NAString = getNA<String>();
  varString.setValues({"hello", "world", NAString, "bar"});
  varString.display_old();

  // ==============================================
  // Testing the VariableDouble class with nCols>1
  // ==============================================

  // Creating a VariableDouble called "test_double_multi"
  mestitle(1, "Testing class VariableDouble with nCols>1 and Unit");
  VariableDouble varDblM("test_double_multi", 3, "m/s");

  // Filling it with values (including a NA) and resetting the number of Cols
  varDblM.setValues({1., 3., 4., 5., NADbl, 4.}, 0);
  varDblM.setValues({3., 4., 5., NADbl, 4., 1.}, 2);
  varDblM.setValues({4., 5., NADbl, 4., 1., 3.}, 1);
  std::cout << "Number of cols = " << varDblM.getNCols() << std::endl;
  std::cout << "Number of samples = " << varDblM.getNSamples() << std::endl;
  std::cout << "Total number of values = " << varDblM.getNAllValues()
            << std::endl;
  std::cout << "Total number of defined values = " << varDblM.getNAllDefined()
            << std::endl;
  for (size_t item = 0; item < varDblM.getNCols(); ++item)
    std::cout << "Number of defined values for item #" << item << " = "
              << varDblM.getNDefined(item) << std::endl;
  varDblM.display_old();

  // ==============================================
  // Testing the VariableString class with nCols>1
  // ==============================================

  // Creating a VariableString called "test_string_multi"
  mestitle(1, "Testing class VariableString with nCols>1");
  VariableString varStringM("test_string_multi", 2);

  // Filling it with a four values (including a NA) and displaying them (old style)
  varStringM.setValues({"hello", "world", NAString, "bar"}, 0);
  varStringM.setValues({"baz", "qux", "foo", "quux"}, 1);

  std::cout << "Number of samples = " << varStringM.getNSamples() << std::endl;
  std::cout << "Number of cols = " << varStringM.getNCols() << std::endl;
  std::cout << "Total number of values = " << varStringM.getNAllValues()
            << std::endl;
  std::cout << "Total number of defined values = "
            << varStringM.getNAllDefined() << std::endl;
  for (size_t item = 0; item < varStringM.getNCols(); ++item)
    std::cout << "Number of defined values for item #" << item << " = "
              << varStringM.getNDefined(item) << std::endl;
  varStringM.display_old();

  // =====================================================
  // Testing the VariableDouble class with Format and Unit
  // =====================================================

  // Creating a VariableDouble called "test_double_multi"
  mestitle(1, "Testing class VariableDouble with Format and Unit");
  VariableDouble varDblF(
    "test_double_format", 1, "m/s", true,
    VariableFormat(
      10, 2, EVariableAlign::RIGHT, EVariableTruncate::RIGHT,
      EVariableNumeric::FIXED));

  // Filling it with values (including a NA) and resetting the number of Cols
  varDblF.setValues({1.311244232, 3.1343, 4.1441, 5.312, NADbl, 4.});
  varDblF.display_old();

  // ================================================
  // Testing the VariableDouble class with Na refused
  // ================================================

  // Creating a VariableDouble called "test_double_multi"
  mestitle(1, "Testing class VariableDouble with NA not allowed");
  VariableDouble varDblNA("test_double_NAnotallowed", 1, "m/s", false);

  varDblNA.setValues({1.311244232, 5.312, NADbl, 4.});
  varDblNA.display_old();

  return 0;
}
