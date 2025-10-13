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

#include "Basic/VectorT.hpp"
#include "Model/GaussianProcess.hpp"
#include "gstlearn_export.hpp"
#include "Db/Db.hpp"
#include "vector"

namespace gstlrn {

class GaussianProcess;
class Db;
class ECov;

using DataFrame = std::map<std::string, std::vector<double>>;

GSTLEARN_EXPORT Db* createDbFromDataFrame(const DataFrame* dat,
                                          const VectorString& coordinates);

GSTLEARN_EXPORT GaussianProcess* createModelFromData(const Db* dat,
                                                     const VectorString& variables,
                                                     const std::vector<ECov>& structs,
                                                     bool addMeasurementError = false);
}