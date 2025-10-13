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
#include "geoslib_define.h"

#include "Matrix/MatrixDense.hpp"
#include "Basic/ASerializable.hpp"

namespace gstlrn
{

/**
 * Stores an array of values as a Table, i.e. a MatrixDense
 * where rows and columns can be optionally decorated
 */
class GSTLEARN_EXPORT Table : public MatrixDense, public ASerializable {

public:
  Table(Id nrow = 0, Id ncol = 0, bool skip_title = false, bool skip_description = false);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixDense)

  void reset(Id nrows, Id ncols) override;

  String toString(const AStringFormat* strfmt = nullptr) const override;

  static Table* create(Id nrow = 0, Id ncol = 0);
  static Table* createFromNames(const VectorString &rownames,
                                const VectorString &colnames);
  static Table* createFromNF(const String& NFFilename, bool verbose = true);
  static Table* createFromTable(const Table& table);

  VectorDouble getRange(Id icol) const;
  VectorDouble getAllRange() const;
  void plot(Id isimu) const;

  void setColumnNames(const VectorString &colNames);
  void setColumnName(Id icol, const String& name);
  void setRowNames(const VectorString &rowNames);
  void setRowName(Id irow, const String& name);

  VectorString getColumnNames() const {  return _colNames; }
  VectorString getRowNames() const {  return _rowNames; }
  String getColumnName(Id icol) const;
  String getRowName(Id irow) const;

  const String& getTitle() const { return _title; }
  void setTitle(const String &title) { _title = title; }
  void setSkipDescription(bool skipDescription) { _skipDescription = skipDescription; }
  void setSkipTitle(bool skipTitle) { _skipTitle = skipTitle; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
 String _getNFName() const override{return "Table"; }

private:
  void _clearDecoration();

private:
  String _title;
  VectorString _rowNames;
  VectorString _colNames;
  bool _skipTitle;
  bool _skipDescription;
};
}
