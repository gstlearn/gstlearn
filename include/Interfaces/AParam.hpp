#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT AParam
{
public:
  AParam(){}
  virtual ~AParam(){}

  virtual bool checkConsistence() const = 0;
private:
};
