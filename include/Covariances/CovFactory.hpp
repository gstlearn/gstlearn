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

#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

namespace gstlrn
{ 
class CovAniso;
class ACovFunc;
class CovContext;

class GSTLEARN_EXPORT CovFactory
{
public:
  static ACovFunc*    createCovFunc(const ECov& type, const CovContext& ctxt);
  static ACovFunc*    duplicateCovFunc(const ACovFunc& cov);
  static void         displayCovList(const CovContext& ctxt);
  static VectorString getCovList(const CovContext& ctxt, Id order=3);
  static ECov         identifyCovariance(const String& cov_name,
                                         const CovContext& ctxt);
  static double       getScaleFactor(const ECov &type, double param);
};
}