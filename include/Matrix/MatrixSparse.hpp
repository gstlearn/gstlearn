/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"

class cs;

/**
 * MatrixSparse wrapper
 */
class GSTLEARN_EXPORT MatrixSparse : public AStringable, public ICloneable
{
public:
  MatrixSparse(int nrow = 0, int ncol = 0);
  MatrixSparse(const MatrixSparse &m);
  MatrixSparse& operator= (const MatrixSparse &m);
  virtual ~MatrixSparse();
  //MatrixSparse(const cs* A);

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSparse);

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  cs*  _cs;

};

