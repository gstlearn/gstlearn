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

#include "Basic/Vector.hpp"
#include "Basic/ASerializable.hpp"

/**
 * Stores the multivariate statistics
 * - each line corresponds to a sample
 * - each column corresponds to a variable
 * The organization stands as a vector (variables) or samples.
 * This allows adding the statistics for all variables for a new sample
 */
class Table: public ASerializable
{
public:
  Table(int nrows = 0, int ncols = 0);
  Table(const VectorVectorDouble& table);
  Table(const String& neutralFileName, bool verbose = false);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

public:
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;

  void init(int nrows, int ncols) { resize(nrows, ncols); }
  bool isEmpty() const { return _stats.empty(); }
  int getRowNumber() const;
  int getColNumber() const;
  VectorDouble getCol(int icol) const;
  VectorDouble getRow(int irow) const;
  void clear() { _stats.clear(); }
  void resize(int irow, int ncols);
  void update(int irow, int icol, double value) { _stats[icol][irow] = value; }
  double getValue(int irow, int icol) const;
  VectorDouble getRange(int icol) const;
  VectorDouble getRange() const;
  void display(int isimu) const;
  void plot(int isimu) const;

private:
  bool _isColValid(int icol) const;
  bool _isRowValid(int irow) const;

private:
  VectorVectorDouble _stats;
};
