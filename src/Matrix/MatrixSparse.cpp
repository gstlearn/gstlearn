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
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <iostream>
#include <iomanip>

// External library
#include "csparse_d.h"
#include "csparse_f.h"

MatrixSparse::MatrixSparse(int nrow, int ncol)
    : AStringable(),
      _cs(nullptr)
{
  cs* Atriplet;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  cs_entry(Atriplet, nrow - 1, ncol - 1, 0.);
  _cs = cs_triplet(Atriplet);
  Atriplet = cs_spfree(Atriplet);
}

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AStringable(m),
      _cs(nullptr)
{
  _cs = cs_duplicate(m._cs);
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _cs = cs_spfree(_cs);
    _cs = cs_duplicate(m._cs);
  }
  return *this;
}

MatrixSparse::~MatrixSparse()
{
  _cs = cs_spfree(_cs);
}

String MatrixSparse::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toMatrix(String(), _cs);

  return sstr.str();
}
