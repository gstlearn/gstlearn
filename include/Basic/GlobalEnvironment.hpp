/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT GlobalEnvironment
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
