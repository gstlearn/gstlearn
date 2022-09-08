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

#include "Matrix/MatrixRectangular.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

/**
 * Stores the multivariate statistics
 * - each line corresponds to a sample
 * - each column corresponds to a variable
 * The organization stands as a vector (variables) or samples.
 * This allows adding the statistics for all variables for a new sample
 */
class GSTLEARN_EXPORT Table: public ASerializable, public AStringable
{
public:
  Table(int nrows = 0, int ncols = 0);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

public:
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromArray(const VectorVectorDouble& table, bool flagByRow = true);

  static Table* create(int nrows = 0, int ncols = 0);
  static Table* createFromNF(const String& neutralFilename, bool verbose = true);
  static Table* createFromArray(const VectorVectorDouble& tabin, bool flagByRow = true);

  void init(int nrows, int ncols) { _tab.reset(nrows, ncols); }
  bool isEmpty() const { return _tab.isEmpty(); }
  int getRowNumber() const;
  int getColumnNumber() const;
  VectorDouble getColumn(int icol) const;
  VectorDouble getRow(int irow) const;
  void addRow();
  void update(int irow, int icol, double value);
  void increment(int irow, int icol, double value);
  double getValue(int irow, int icol) const;
  void setValue(int irow, int icol, double value);
  VectorDouble getRange(int icol) const;
  VectorDouble getAllRange() const;
  void plot(int isimu) const;
  void fill(double valinit = 0);

  void setColumnNames(const VectorString &colNames) { _colNames = colNames; }
  void setColumnName(int icol, const String& name);
  void setRowNames(const VectorString &rowNames) { _rowNames = rowNames; }
  void setRowName(int irow, const String& name);

  VectorString getColumnNames() const {  return _colNames; }
  VectorString getRowNames() const {  return _rowNames; }
  String getColumnName(int icol) const {  return _colNames[icol]; }
  String getRowName(int irow) const {  return _rowNames[irow]; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Table"; }

private:
  bool _isColValid(int icol) const;
  bool _isRowValid(int irow) const;

private:
  MatrixRectangular _tab;
  VectorString _rowNames;
  VectorString _colNames;
};
