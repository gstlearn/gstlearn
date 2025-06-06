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

class ModelGeneric;

class GSTLEARN_EXPORT AModelOptimNew 
{
public:
  AModelOptimNew(const ModelGeneric* model = nullptr)
    : _model(model)
  {
  }
  AModelOptimNew(const AModelOptimNew& r)
    : _model(r._model)
  {
  }
  AModelOptimNew& operator=(const AModelOptimNew& r)
  {
    if (this != &r)
    {
      _model = r._model;
    }
    return *this;
  }
  virtual ~AModelOptimNew() = default;
  
  virtual double computeCost(bool verbose = false) = 0;

protected:
  const ModelGeneric* _model;
};
