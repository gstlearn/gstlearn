/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
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
  Table(int nrow = 0, int ncol = 0);
  Table(const Table &m);
  Table& operator= (const Table &m);
  virtual ~Table();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixRectangular)

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Table* create(int nrow = 0, int ncol = 0);
  static Table* createFromNames(const VectorString &rownames,
                                const VectorString &colnames);
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

  const String& getTitle() const { return _title; }
  void setTitle(const String &title) { _title = title; }
  void setSkipDescription(bool skipDescription) { _skipDescription = skipDescription; }
  void setSkipTitle(bool skipTitle) { _skipTitle = skipTitle; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Table"; }
  void    _clearContents() override;

private:
  String _title;
  VectorString _rowNames;
  VectorString _colNames;
  bool _skipTitle;
  bool _skipDescription;
};
