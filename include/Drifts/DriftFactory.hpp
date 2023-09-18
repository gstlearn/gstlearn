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
  static ADriftElem* createDriftByType(const EDrift &type,
                                       int rank_fex = 0,
                                       const CovContext &ctxt = CovContext());
  static ADriftElem* createDriftBySymbol(const String &symbol,
                                         const CovContext &ctxt);
  static void        displayList(const CovContext& ctxt);
 };
