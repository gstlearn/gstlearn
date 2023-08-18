/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT IProjMatrix
{
public:
  IProjMatrix() { }
  virtual ~IProjMatrix() { }
  virtual int point2mesh(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual int mesh2point(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual int getApexNumber() const = 0;
  virtual int getPointNumber() const = 0;
};
