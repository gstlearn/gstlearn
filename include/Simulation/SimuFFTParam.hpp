/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT SimuFFTParam: public AStringable
{
public:
  SimuFFTParam(bool flag_aliasing = true, double percent = 0.1);
  SimuFFTParam(const SimuFFTParam &r);
  SimuFFTParam& operator=(const SimuFFTParam &r);
  virtual ~SimuFFTParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  bool isFlagAliasing() const { return _flagAliasing; }
  void setFlagAliasing(bool flagAliasing) { _flagAliasing = flagAliasing; }
  double getPercent() const { return _percent; }
  void setPercent(double percent) { _percent = percent; }

private:
  bool _flagAliasing;
  double _percent;
};
