/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT SimuPartitionParam: public AStringable
{
public:
  SimuPartitionParam(Id nbtuba = 100,
                     double intensity = 0.1,
                     const VectorDouble& dilate = VectorDouble());
  SimuPartitionParam(const SimuPartitionParam &r);
  SimuPartitionParam& operator=(const SimuPartitionParam &r);
  virtual ~SimuPartitionParam();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  VectorDouble getDilate() const { return _dilate; }
  void setDilate(const VectorDouble& dilate) { _dilate = dilate; }
  double getIntensity() const { return _intensity; }
  void setIntensity(double intensity) { _intensity = intensity; }
  Id getNbtuba() const { return _nbtuba; }
  void setNbtuba(Id nbtuba) { _nbtuba = nbtuba; }
  double getDilate(Id idim) const;

private:
  Id _nbtuba;
  double _intensity;
  VectorDouble _dilate;
};
}