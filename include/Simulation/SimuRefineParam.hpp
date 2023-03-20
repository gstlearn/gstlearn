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
