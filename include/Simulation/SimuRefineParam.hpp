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

class GSTLEARN_EXPORT SimuRefineParam: public AStringable
{
public:
  SimuRefineParam(int nmult = 1, bool flag_SK = true);
  SimuRefineParam(const SimuRefineParam &r);
  SimuRefineParam& operator=(const SimuRefineParam &r);
  virtual ~SimuRefineParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  bool isFlagSK() const { return _flagSK; }
  void setFlagKs(bool flagKS) { _flagSK = flagKS; }
  int getNmult() const { return _nmult; }
  void setNmult(int nmult) { _nmult = nmult; }

private:
  int _nmult;
  bool _flagSK;
};
