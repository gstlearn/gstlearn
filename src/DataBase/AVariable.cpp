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
#include "DataBase/AVariable.hpp"
#include "DataBase/VariableBool.hpp"
#include "DataBase/VariableDouble.hpp"
#include "DataBase/VariableInt.hpp"
#include "DataBase/VariableString.hpp"

using namespace gstlrn;

static bool _flagVariableCheck = false;

void AVariable::_setFlagVariableCheck(bool flagVariableCheck)
{
  _flagVariableCheck = flagVariableCheck;
}

bool AVariable::_getFlagVariableCheck()
{
  return _flagVariableCheck;
}

/**
 * Create a variable with appropriate type
 *
 * @param[in] name   name of the created variable.
 * @param[in] type   String describing the desired type
 *                   (int, bool, double, string)
 *
 * @remark: if type does not match any existing type, an exception is throw
 */
std::unique_ptr<AVariable>
  AVariable::create(const String& name, const String& type)
{
  if (type == "int")
  {
    return std::make_unique<VariableInt>(name);
  }
  if (type == "bool")
  {
    return std::make_unique<VariableBool>(name);
  }
  if (type == "double")
  {
    return std::make_unique<VariableDouble>(name);
  }
  if (type == "string")
  {
    return std::make_unique<VariableString>(name);
  }

  throw("Wrong variable type");
}

// namespace gstlrn
