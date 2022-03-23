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
#include "geoslib_define.h"

#include "Anamorphosis/AnamContinuous.hpp"
#include "Basic/ASerializable.hpp"

class ECalcMember;

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

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AAnam Interface
  const EAnam& getType() const override { return EAnam::EXTERNAL; }

  /// AnamContinuous Interface
  void   calculateMeanAndVariance() override;
  double TransformToRawValue(double h) const override;
  double RawToTransformValue(double h) const override;

  void setY2zFunction(double (*y2z_function)(double)) { _y2z_function = y2z_function; }
  void setZ2yFunction(double (*z2y_function)(double)) { _z2y_function = z2y_function; }

protected:
  virtual int _deserialize(std::istream& is, bool verbose) override;
  virtual int _serialize(std::ostream& os, bool verbose = false) const override;
};
