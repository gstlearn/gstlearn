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
 * Stores an array of values as a Table, i.e. a MatrixRectangular
 * where rows and columns can be optionally decorated
 */

class GSTLEARN_EXPORT Table : public MatrixRectangular, public ASerializable {

public:
  Table(int nrows = 0, int ncols = 0);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixRectangular)

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Table* create(int nrows = 0, int ncols = 0);
  static Table* createFromNF(const String& neutralFilename, bool verbose = true);

  VectorDouble getRange(int icol) const;
  VectorDouble getAllRange() const;
  void plot(int isimu) const;

  void setColumnNames(const VectorString &colNames);
  void setColumnName(int icol, const String& name);
  void setRowNames(const VectorString &rowNames);
  void setRowName(int irow, const String& name);

  VectorString getColumnNames() const {  return _colNames; }
  VectorString getRowNames() const {  return _rowNames; }
  String getColumnName(int icol) const;
  String getRowName(int irow) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Table"; }
  void    _clearContents() override;

private:
  VectorString _rowNames;
  VectorString _colNames;
};
