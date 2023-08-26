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

class GSTLEARN_EXPORT DriftFactory
{
public:
  static ADriftElem* createDriftFunc(const EDrift &type,
                                     const CovContext &ctxt = CovContext(),
                                     int rank_fex = 0);
  static ADriftElem*  duplicateDriftFunc(const ADriftElem& cov);
  static void         displayList(const CovContext& ctxt);
  static EDrift       identifyDrift(const String& symbol, int* rank, const CovContext& ctxt);
 };
