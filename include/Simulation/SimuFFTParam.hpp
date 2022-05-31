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
