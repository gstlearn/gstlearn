#pragma once

#include "gstlearn_export.hpp"

#include "LinearOp/IProjMatrix.hpp"

class GSTLEARN_EXPORT TestInheritance
{
public:
  TestInheritance();
  TestInheritance(const TestInheritance& r) = delete;
  TestInheritance& operator=(const TestInheritance& r) = delete;
  void setIproj(IProjMatrix* ipr){ _iproj = ipr;}
  virtual ~TestInheritance();

private:
 IProjMatrix* _iproj;
};
