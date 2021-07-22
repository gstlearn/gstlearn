/*
 * GlobalEnvironment.cpp
 *
 *  Created on: 22 juil. 2021
 *      Author: drenard
 */

#include "Basic/GlobalEnvironment.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

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

void GlobalEnvironment::setDomainReference(int value, bool verbose)
{
  if (value < 0) value = 0;
  _domainReference = value;
  if (_domainReference == 0) return;
  if (verbose) printDomainReference();
  return;
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
bool GlobalEnvironment::matchDomainReference(double value)
{
  if (FFFF(value)) return 0;
  if ((int) value == _domainReference) return true;
  return false;
}

