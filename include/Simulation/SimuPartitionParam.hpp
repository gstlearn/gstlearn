/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
