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

#include "Basic/VectorNumT.hpp"
#include "Spatial/Projection.hpp"

class Db;
class MatrixRectangular;

typedef struct{
  VectorDouble center; // Vector for the center of gravity
  VectorDouble mvalues; // Vector of eigen values (normalized)
  MatrixRectangular mvectors; // Array of eigen vectors
  MatrixRectangular axesEllipse;
  MatrixRectangular axesInertia;
  double inertia; // Value of the inertia
  double wztot; // Sum of weights
  double iso; // Iso index
  int nvalid; // Number of valid samples

  double theta;
  double ra;
  double rb;

  VectorDouble xl1;
  VectorDouble yl1;
  VectorDouble xl2;
  VectorDouble yl2;
} cgiOutput;

GSTLEARN_EXPORT void cgiPrint(const cgiOutput& cgiouput);
GSTLEARN_EXPORT cgiOutput cgi(Db *db, const String &name = "");
GSTLEARN_EXPORT int spatial(Db *db,
                            double *totab,
                            double *parea,
                            double *eqarea);
