/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

/**
 * Stores the multivariate statistics
 * - each line corresponds to a sample
 * - each column corresponds to a variable
 * The organization stands as a vector (variables) or samples.
 * This allows adding the statistics for all variables for a new sample
 */
class GSTLEARN_EXPORT Table: public ASerializable ,public AStringable
{
public:
  Table(int nrows = 0, int ncols = 0);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

public:
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromArray(const VectorVectorDouble& table);

  static Table* createFromNF(const String& neutralFileName, bool verbose = false);
  static Table* createFromArray(const VectorVectorDouble& tabin);

  void init(int nrows, int ncols, bool zero = false) { resize(nrows, ncols, zero); }
  bool isEmpty() const { return _stats.empty(); }
  int getRowNumber() const;
  int getColNumber() const;
  VectorDouble getCol(int icol) const;
  VectorDouble getRow(int irow) const;
  void clear() { _stats.clear(); }
  void resize(int irow, int ncols, bool zero=false);
  void update(int irow, int icol, double value)    { _stats[icol][irow] = value; }
  void increment(int irow, int icol, double value) { _stats[icol][irow] += value; }
  double getValue(int irow, int icol) const;
  void setValue(int irow, int icol, double value);
  VectorDouble getRange(int icol) const;
  VectorDouble getAllRange() const;
  void plot(int isimu) const;

private:
  bool _isColValid(int icol) const;
  bool _isRowValid(int irow) const;

private:
  VectorVectorDouble _stats;
};
