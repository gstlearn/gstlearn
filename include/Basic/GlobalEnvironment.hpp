/*
 * GlobalEnvironment.hpp
 *
 *  Created on: 22 juil. 2021
 *      Author: drenard
 */

#pragma once

class GlobalEnvironment
{
private:
  static GlobalEnvironment* _env;
  GlobalEnvironment();
  virtual ~GlobalEnvironment();

public:
  static GlobalEnvironment* getEnv();

  bool isDomainReference() const { return _domainReference > 0; }
  int  getDomainReference() const { return _domainReference; }
  void setDomainReference(int domainReference, bool verbose = false);
  void printDomainReference(void) const;
  bool matchDomainReference(double value);

private:
  int _domainReference;
};
