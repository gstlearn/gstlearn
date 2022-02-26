#pragma once

#include "gstlearn_export.hpp"
#include "Interfaces/interface_d.hpp"

#include <iostream>
#include <string>

class GSTLEARN_EXPORT Category
{
public:
  Category(int value = UNDEF_CAT_VAL, const String& label = UNDEF_CAT_LABEL);
  virtual ~Category();

  const String& getLabel() const;
  int   getValue() const;

#ifndef SWIG
  //conversion operator
  explicit operator int()
  {
    return(_value);
  }
  explicit operator double()
  {
    return(_value);
  }

  friend std::ostream& operator<<(std::ostream& stream, const Category& cat);
#endif

private:
  int   _value;
  String _label;
};
