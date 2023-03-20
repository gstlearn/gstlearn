/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

class CovAniso;
class ACovFunc;
class CovContext;

class GSTLEARN_EXPORT CovFactory
{
public:
  static ACovFunc*    createCovFunc(const ECov& type, const CovContext& ctxt);
  static ACovFunc*    duplicateCovFunc(const ACovFunc& cov);
  static void         displayList(const CovContext& ctxt);
  static VectorString getCovList(const CovContext& ctxt, int order=3);
  static ECov         identifyCovariance(const String& cov_name,
                                         const CovContext& ctxt);
};

