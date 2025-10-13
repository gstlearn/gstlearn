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
#include "Basic/GlobalEnvironment.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

namespace gstlrn
{
GlobalEnvironment* GlobalEnvironment::_env = nullptr;

GlobalEnvironment::GlobalEnvironment()
  : _domainReference(0)
{

}

GlobalEnvironment::~GlobalEnvironment()
{
  // TODO delete singleton?
}

GlobalEnvironment* GlobalEnvironment::getEnv()
{
  if (GlobalEnvironment::_env == nullptr)
    GlobalEnvironment::_env = new GlobalEnvironment();
  return _env;
}

void GlobalEnvironment::setDomainReference(Id domainReference, bool verbose)
{
  if (domainReference < 0) domainReference = 0;
  _domainReference = domainReference;
  if (_domainReference == 0) return;
  if (verbose) printDomainReference();
}

void GlobalEnvironment::printDomainReference(void) const
{
  if (_domainReference > 0)
  {
    mestitle(1, "Parameters for Domaining");
    message("Domain Reference value = %d\n", _domainReference);
    message("Use 'domain.define' to modify or cancel the Domaining\n");
  }
}

/****************************************************************************/
/*!
 **  Check if the Domain value matches the Reference value for the Domain
 **
 ** \param[in]  value    Reference Domain value
 **
 *****************************************************************************/
bool GlobalEnvironment::matchDomainReference(double value) const
{
  if (FFFF(value)) return 0;
  if (static_cast<Id>(value) == _domainReference) return true;
  return false;
}
}
