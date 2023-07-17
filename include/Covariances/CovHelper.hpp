/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT CovHelper
{
public:
  CovHelper() {};
  ~CovHelper() {} ;
  CovHelper(const CovHelper&) = delete;
  CovHelper& operator=(const CovHelper&) = delete;

  static VectorString getAllCovariances(int ndim=2, int minorder=-1, bool hasrange=false, bool flagSimtub=true);
};

