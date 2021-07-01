#include "Interfaces/Category.hpp"

Category::Category(int value, const std::string& label): _value(value),_label(label)
{
}

Category::~Category()
{
}

const std::string& Category::getLabel() const
{
  return(_label);
}

int Category::getValue() const
{
  return(_value);
}

#ifndef SWIG
std::ostream& operator<<(std::ostream &stream, const Category& cat)
{
  stream<< cat.getLabel();
  return(stream);
}
#endif

