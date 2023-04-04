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
