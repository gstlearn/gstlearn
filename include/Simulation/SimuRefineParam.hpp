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

namespace gstlrn
{
class GSTLEARN_EXPORT SimuRefineParam: public AStringable
{
public:
  SimuRefineParam(Id nmult = 1, bool flag_SK = true);
  SimuRefineParam(const SimuRefineParam &r);
  SimuRefineParam& operator=(const SimuRefineParam &r);
  virtual ~SimuRefineParam();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  bool isFlagSK() const { return _flagSK; }
  void setFlagKs(bool flagKS) { _flagSK = flagKS; }
  Id getNmult() const { return _nmult; }
  void setNmult(Id nmult) { _nmult = nmult; }

private:
  Id _nmult;
  bool _flagSK;
};
}