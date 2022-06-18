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
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Fractures/Family.hpp"
#include "Fractures/Fault.hpp"

class GSTLEARN_EXPORT Description: public AStringable
{
public:
  Description();
  Description(const Description& r);
  Description& operator=(const Description& r);
  virtual ~Description();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNPoint() const { return (int) _x.size(); }

  int getFamily() const { return _family; }
  void setFamily(int family) { _family = family; }
  double getOrient() const { return _orient; }
  void setOrient(double orient) { _orient = orient; }
  double getXXF(int i) const { return _x[i]; }
  double getYYF(int i) const { return _y[i]; }
  void setXXF(int i, double value) { _x[i] = value; }
  void setYYF(int i, double value) { _y[i] = value; }

  void addPoint(double x, double y);

private:
  int _family;
  double _orient;
  VectorDouble _x;
  VectorDouble _y;
};
