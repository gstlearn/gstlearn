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
#include "Model/ModelGeneric.hpp"

ModelGeneric::ModelGeneric(const CovContext &ctxt)
    : _cova(nullptr),
      _driftList(nullptr),
      _ctxt(ctxt)
{
}

ModelGeneric::~ModelGeneric()
{
}


void ModelGeneric::setField(double field)
{
  _ctxt.setField(field);
  if (_cova != nullptr)
    _cova->setContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

bool ModelGeneric::isValid() const
{
  return _isValid();
}

bool ModelGeneric::_isValid() const
{
  return true;
}