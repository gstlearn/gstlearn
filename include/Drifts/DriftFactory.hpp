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

#include "Covariances/CovContext.hpp"

namespace gstlrn
{
class ADrift;
class DriftList;

class GSTLEARN_EXPORT DriftFactory
{
public:
  static ADrift* createDriftByRank(Id rank, Id rank_fex);
  static ADrift* createDriftBySymbol(const String &symbol);
  static ADrift* createDriftByIdentifier(const String &driftname);
  static DriftList* createDriftListFromIRF(Id order = 0,
                                           Id nfex = 0,
                                           const CovContext &ctxt = CovContext());
  static DriftList* createDriftListForGradients(const DriftList* olddrifts,
                                                const CovContext& ctxt = CovContext());
};
}