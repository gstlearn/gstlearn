#include "Interfaces/AVariable.hpp"
#include "Interfaces/VariableDouble.hpp"
#include "Interfaces/VariableBool.hpp"
#include "Interfaces/VariableInt.hpp"
#include "Interfaces/VariableString.hpp"

AVariable::AVariable()
    : _name(UNDEF_STRING)
{
}

AVariable::AVariable(const String &name)
    : _name(name)
{
}

AVariable::AVariable(const AVariable &ref)
    : _name(ref._name)
{
}

AVariable::~AVariable()
{
}

AVariable& AVariable::operator=(const AVariable &ref)
{
  if (this != &ref)
  {
    _name = ref._name;
  }
  return (*this);
}

void AVariable::setName(const String&name)
{
  _name = name;
}

const String& AVariable::getName() const
{
  return _name;
}

/**
 * Create a variable with appropriate type
 *
 * @param[in] type   String describing the desired type
 *                   (int, bool, double, string)
 * @param[in] name   name of the created variable.
 *
 * @remark: if type does not match any existing type, an exception is throw
 */
AVariable* AVariable::createVariable(const String& type, const String& name)
{
  AVariable* new_var;
  if (!type.compare("int"))
  {
    new_var = new VariableInt(name);
  }
  else if (!type.compare("bool"))
  {
    new_var = new VariableBool(name);
  }
  else if (!type.compare("double"))
  {
    new_var = new VariableDouble(name);
  }
  else if (!type.compare("string"))
  {
    new_var = new VariableString(name);
  }
  else
  {
    throw("Wrong variable type");
  }
  return new_var;
}
