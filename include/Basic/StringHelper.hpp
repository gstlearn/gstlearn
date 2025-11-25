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

#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

// TODO : add Namespace
#define SPACES " \t\r\n"

//
// This contains all the "local" methods which are used in the étemplate" methods
// of the class String.hhp
//

#define CASE_DOUBLE 0
#define CASE_REAL   1
#define CASE_INT    2
#define CASE_COL    3
#define CASE_ROW    4

namespace gstlrn
{
#ifndef SWIG
// Functions for returning parameters from OptCst environment
GSTLEARN_EXPORT Id _getColumnRank();
GSTLEARN_EXPORT Id _getColumnName();
GSTLEARN_EXPORT Id _getColumnSize(Id localSize = 0, Id nColumns = 1);
GSTLEARN_EXPORT Id _getMaxNCols();
GSTLEARN_EXPORT Id _getMaxNRows();
GSTLEARN_EXPORT Id _getNBatch();
GSTLEARN_EXPORT Id _getDecimalNumber();
GSTLEARN_EXPORT double _getThresh();

GSTLEARN_EXPORT String _toStrRowColumn(Id icase, Id value, Id flagAdd);

GSTLEARN_EXPORT String _toStrRowHeader(const VectorString& rownames, Id iy, Id rowSize = 0);

GSTLEARN_EXPORT String _toStrColumnHeaders(const VectorString& colnames,
                                           Id colfrom,
                                           Id colto,
                                           Id colSize = 0);

GSTLEARN_EXPORT String _toStrTrailer(Id ncols, Id nrows, Id ncols_util, Id nrows_util);

GSTLEARN_EXPORT Id _convertToId(const String& v, char dec = '.');
GSTLEARN_EXPORT double _convertToDouble(const String& v, char dec = '.');

GSTLEARN_EXPORT Id _askInt(const String& v, Id defval = getNA<Id>(), bool authTest = true);
GSTLEARN_EXPORT double _askDouble(const String& v, double defval = getNA<double>(), bool authTest = true);
GSTLEARN_EXPORT bool _askBool(const String& v, bool defval = false, bool authTest = true);

#endif

} // namespace gstlrn