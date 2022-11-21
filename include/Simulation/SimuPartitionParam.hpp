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
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT SimuPartitionParam: public AStringable
{
public:
  SimuPartitionParam(int nbtuba = 100,
                     double intensity = 0.1,
                     const VectorDouble& dilate = VectorDouble());
  SimuPartitionParam(const SimuPartitionParam &r);
  SimuPartitionParam& operator=(const SimuPartitionParam &r);
  virtual ~SimuPartitionParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  const VectorDouble getDilate() const { return _dilate; }
  void setDilate(const VectorDouble& dilate) { _dilate = dilate; }
  double getIntensity() const { return _intensity; }
  void setIntensity(double intensity) { _intensity = intensity; }
  int getNbtuba() const { return _nbtuba; }
  void setNbtuba(int nbtuba) { _nbtuba = nbtuba; }
  double getDilate(int idim) const;

private:
  int _nbtuba;
  double _intensity;
  VectorDouble _dilate;
};
