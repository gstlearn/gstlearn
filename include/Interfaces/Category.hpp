#ifndef CATEGORY_HPP
#define CATEGORY_HPP

#include <iostream>
#include <string>

#include "Interfaces/interface_d.hpp"  //USE

class Category
{
public:
  Category(int value = UNDEF_CAT_VAL, const std::string& label = UNDEF_CAT_LABEL);
  virtual ~Category();

  const std::string& getLabel() const;
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
  std::string _label;
};

#endif
