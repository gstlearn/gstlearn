/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT AFunction
{
public:
  AFunction();
  virtual ~AFunction(){};

public :
  virtual double eval(double x) const {return x;};

};
