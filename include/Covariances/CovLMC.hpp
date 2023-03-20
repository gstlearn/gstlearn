/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovAniso;

class GSTLEARN_EXPORT CovLMC : public ACovAnisoList
{
public:
  CovLMC(const ASpace* space = nullptr);
  CovLMC(const CovLMC &r);
  CovLMC(const ACovAnisoList &r);
  CovLMC& operator= (const CovLMC &r);
  virtual ~CovLMC();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMC)
};
