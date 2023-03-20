/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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

