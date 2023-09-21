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

#include "Enum/EDrift.hpp"

#include "Covariances/CovContext.hpp"

class ADriftElem;
class DriftList;

class GSTLEARN_EXPORT DriftFactory
{
public:
  static ADriftElem* createDriftByRank(int rank,
                                       int rank_fex,
                                       const CovContext &ctxt);
  static ADriftElem* createDriftBySymbol(const String &symbol,
                                         const CovContext &ctxt = CovContext());
  static ADriftElem* createDriftByIdentifier(const String &driftname,
                                             const CovContext &ctxt = CovContext());
  static DriftList* createDriftListFromIRF(int order = 0,
                                           int nfex = 0,
                                           const CovContext &ctxt = CovContext());
  static DriftList* createDriftListForGradients(const DriftList* inputlist,
                                                const CovContext &ctxt = CovContext());
};
