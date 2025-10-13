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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{ 
class GSTLEARN_EXPORT GlobalEnvironment
{
private:
  static GlobalEnvironment* _env;
  GlobalEnvironment();
  virtual ~GlobalEnvironment();

public:
  static GlobalEnvironment* getEnv();

  bool isDomainReference() const { return _domainReference > 0; }
  Id  getDomainReference() const { return _domainReference; }
  void setDomainReference(Id domainReference, bool verbose = false);
  void printDomainReference(void) const;
  bool matchDomainReference(double value) const;

private:
  Id _domainReference;
};
} // namespace gstlrn
