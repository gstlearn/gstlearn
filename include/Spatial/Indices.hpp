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

#include "Basic/VectorNumT.hpp"
#include "Spatial/Projection.hpp"

class Db;
class MatrixRectangular;

GSTLEARN_EXPORT int cgi(Db *db,
                        int ivar,
                        VectorDouble& center,
                        VectorDouble& mvalue,
                        MatrixRectangular& mvector,
                        double *inertia,
                        double *wztot);
GSTLEARN_EXPORT int spatial(Db *db,
                            double *totab,
                            double *parea,
                            double *eqarea);
