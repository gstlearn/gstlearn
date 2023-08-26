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
#include "Covariances/CovLMC.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"

CovLMC::CovLMC(const ASpace* space)
    : ACovAnisoList(space)
{
}

CovLMC::CovLMC(const ACovAnisoList &r)
    : ACovAnisoList(r)
{
}

CovLMC::CovLMC(const CovLMC &r)
    : ACovAnisoList(r)
{
}

CovLMC& CovLMC::operator=(const CovLMC &r)
{
  if (this != &r)
  {
    ACovAnisoList::operator=(r);
  }
  return *this;
}

CovLMC::~CovLMC()
{
  /// TODO : Delete pointers ?
}

