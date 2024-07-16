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
#include "LinearOp/PrecisionOpMulti.hpp"

PrecisionOpMulti::PrecisionOpMulti(Model* model, const std::vector<AMesh*>& meshes)
  : _isValid(false),
    _covList(),
    _nmeshList(),
    _invSills(),
    _model(nullptr),
    _meshes()
{
  if (! _isValidModel(model)) return;

  if (! _isValidMeshes(meshes)) return;

  if (! _matchModelAndMeshes()) return;

  _isValid = true;
}

PrecisionOpMulti::~PrecisionOpMulti()
{
}

/*****************************************************************************/
/*!
** Check that 'model' is a valid one (composed of Bessel_K and Nugget Effect)
**
** \param[in]  model  Model structure
**
** \remarks If the Model is valid, the vector '_covList' is built
**
*******************************************************************************/
bool PrecisionOpMulti::_isValidModel(Model* model)
{
  if (model == nullptr) return false;

  _covList.clear();
  for (int icov = 0, ncov = model->getCovaNumber(); icov < ncov; icov++)
  {
    if (model->getCovaType(icov) == ECov::NUGGET) continue;
    if (model->getCovaType(icov) == ECov::BESSEL_K)
      _covList.push_back(icov);
    else
    {
      messerr("The covariance type %s is not authorized", model->getCovName(icov).c_str());
      return false;
    };
  }

  // When the Model is valid, it is attached as a member
  _model = model;
  return true;
}

/*****************************************************************************/
/*!
** Check that 'meshes' is a valid argument
** (the number of meshes is equal to the number of covariances)
**
** \param[in]  meshes  Vector of Meshes
**
** \remarks When valid, the member '_nmeshList' is constructed
**
*******************************************************************************/
bool PrecisionOpMulti::_isValidMeshes(const std::vector<AMesh*>& meshes)
{
  if (meshes.empty()) return false;
  int nmeshes = (int) meshes.size();

  // Check that all meshes share the same number of vertices
  _nmeshList.clear();
  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    int nmesh = meshes[imesh]->getNMeshes();
    _nmeshList.push_back(nmesh);
  }

    // When the meshes are valid, they are attached as a member
  _meshes = meshes; 
  return true;
}

bool PrecisionOpMulti::_matchModelAndMeshes() const
{
  return (_getNCov() == _getNMesh());
}

int PrecisionOpMulti::_getNVar() const
{
  if (_model == nullptr) return 0;
  return _model->getVariableNumber();
}

int PrecisionOpMulti::_getNCov() const
{
  if (_covList.empty()) return 0;
  return (int) _covList.size();
}

int PrecisionOpMulti::_getNMesh() const
{
  if (_nmeshList.empty()) return 0;
  return (int) _nmeshList.size();
}

int PrecisionOpMulti::_getSize() const
{
  int nvar = _getNVar();
  int ncov = _getNCov();
  if (ncov != _getNMesh()) return 0;

  int size = 0;
  for (int i = 0; i < ncov; i++)
  {
    size += nvar * _nmeshList[i];
  }
  return size;
}

int PrecisionOpMulti::setModel(Model* model)
{
  if (! _isValidModel(model)) return 1;

  _isValid = _matchModelAndMeshes();

  return 0;
}

int PrecisionOpMulti::setMeshes(const std::vector<AMesh*>& meshes)
{
  if (! _isValidMeshes(meshes)) return 1;

  _isValid = _matchModelAndMeshes();
  return 0;
}

int PrecisionOpMulti::_buildInvSills()
{
  if (_model == nullptr) return 1;

  int ncov = _getNCov();
  _invSills.clear();

  for (int icov = 0; icov < ncov; icov++)
  {
    MatrixSquareSymmetric invs = _model->getSillValues(icov);
    if (invs.invert() != 0) return 1;
    _invSills.push_back(invs);
  }

  return 0;
}