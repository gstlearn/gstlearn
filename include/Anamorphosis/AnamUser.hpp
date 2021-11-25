/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT AnamUser: public AnamContinuous
{
private:
  double (*_y2z_function)(double);
  double (*_z2y_function)(double);

public:
  AnamUser();
  AnamUser(const AnamUser &m);
  AnamUser& operator= (const AnamUser &m);
  virtual ~AnamUser();

  void calculateMeanAndVariance() override;
  virtual String toString(int level) const override;

  double GaussianToRawValue(double h) const override
  {
    if (_y2z_function == nullptr) return TEST;
    return _y2z_function(h);
  }

  double RawToGaussianValue(double h) const override
  {
    if (_z2y_function == nullptr) return TEST;
    return _z2y_function(h);
  }

  void setY2zFunction(double (*y2z_function)(double))
  {
    _y2z_function = y2z_function;
  }
  void setZ2yFunction(double (*z2y_function)(double))
  {
    _z2y_function = z2y_function;
  }

};
