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

PrecisionOpMulti::PrecisionOpMulti(Model* model,
                                   const std::vector<AMesh*>& meshes)
  : _isValid(false)
  , _covList()
  , _pops()
  , _invSills()
  , _cholSills()
  , _model(nullptr)
  , _meshes()
{
  if (! _isValidModel(model)) return;

  if (! _isValidMeshes(meshes)) return;

  if (! _matchModelAndMeshes()) return;

  _isValid = true;
}

PrecisionOpMulti::~PrecisionOpMulti()
{
  _popsClear();
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
      messerr("The covariance type %s is not authorized",
              model->getCovName(icov).c_str());
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
*******************************************************************************/
bool PrecisionOpMulti::_isValidMeshes(const std::vector<AMesh*>& meshes)
{
  if (meshes.empty()) return false;

  _meshes = meshes;

  return true;
}

void PrecisionOpMulti::_popsClear()
{
  for (int i = 0, n = (int) _pops.size(); i < n; i++)
    delete _pops[i];
}

bool PrecisionOpMulti::_matchModelAndMeshes()
{
  if (_getNCov() != _getNMesh()) return false;

  // Create the vector of PrecisionOp
  _popsClear();
  for (int i = 0, number = _getNCov(); i < number; i++)
  {
    // Temporarily normalize the covariance (sill) before instantiating the PrecisionOp
    double sill = _model->getSill(_covList[i], 0, 0);
    _model->setSill(_covList[i],0,0,1.);
    _pops.push_back(PrecisionOp::create(_meshes[i], _model, _covList[i]));
    _model->setSill(_covList[i], 0, 0, sill);
  }
  return true;
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
  if (_meshes.empty()) return 0;
  return (int) _meshes.size();
}

int PrecisionOpMulti::sizes() const
{
  int nvar = _getNVar();
  int ncov = _getNCov();
  if (ncov != _getNMesh()) return 0;

  int size = 0;
  for (int i = 0; i < ncov; i++)
  {
    size += nvar * _meshes[i]->getNApices();
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

void PrecisionOpMulti::clearMeshes()
{
  _meshes.clear();
}
  
void PrecisionOpMulti::addMesh(AMesh* mesh)
{
  if (mesh == nullptr) return;
  _meshes.push_back(mesh);

  _isValid = _matchModelAndMeshes();
}

int PrecisionOpMulti::_buildInvSills()
{
  if (_model == nullptr) return 1;

  int ncov = _getNCov();

  // Do nothing if the array has already been calculated (correct dimension)
  if (ncov == (int) _invSills.size()) return 0;

  for (int icov = 0; icov < ncov; icov++)
  {
    MatrixSquareSymmetric invs = _model->getSillValues(icov);
    if (invs.invert() != 0) return 1;
    _invSills.push_back(invs);
  }
  return 0;
}

int PrecisionOpMulti::_buildCholSills()
{
  if (_model == nullptr) return 1;

  int ncov = _getNCov();

  // Do nothing if the array has already been calculated (correct dimension)
  if (ncov == (int)_cholSills.size()) return 0;

  for (int icov = 0; icov < ncov; icov++)
  {
    MatrixSquareSymmetric chols = _model->getSillValues(icov);
    if (chols.computeCholesky() != 0) return 1;
    _cholSills.push_back(chols);
  }
  return 0;
}

int PrecisionOpMulti::size(int imesh) const 
{ 
  if (imesh < 0 || imesh >= _getNMesh()) return 0;
  return _meshes[imesh]->getNApices();
}

String PrecisionOpMulti::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);

  std::stringstream sstr;

  sstr << "Number of Variables   = " << _getNVar()  << std::endl;
  sstr << "Number of Covariances = " << _getNCov()  << std::endl;
  sstr << "Number of Meshes      = " << _getNMesh() << std::endl;
  sstr << "Vector dimension      = " << sizes()     << std::endl;

  sstr << "Indices of Matérn Covariance = " << VH::toStringAsVI(_covList);

  sstr << "Dimensions of the Meshes = ";
  for (int imesh = 0, nmesh = _getNMesh(); imesh < nmesh; imesh++)
    sstr << size(imesh) << " ";
  sstr << std::endl;

  if (_isValid)
   sstr << std::endl << "Class is Valid for operations!" << std::endl;
  return sstr.str();
}

/**
 * Evaluate the product of this phantom matrix by the input vector
 * @param vecin Input array
 */
VectorDouble PrecisionOpMulti::evalDirect(const VectorDouble& vecin)
{
  int totsize = sizes();
  if ((int)vecin.size() != totsize)
  {
    messerr("'vecin' (%d) should be of dimension (%d)", (int)vecin.size(), totsize);
    return VectorDouble();
  }
  // Invert the matrices of sills
  if (_buildInvSills() != 0)
  {
    messerr("Problem when inverting the matrices of sills");
    return VectorDouble();
  }

  // Blank out the output vector
  VectorDouble vecout(totsize, 0.);

  int nvar = _getNVar();
  int ncov = _getNCov();

  int iad_x = 0;
  int iad_y = 0;
  for (int icov = 0; icov < ncov; icov++)
  {
    int napices = size(icov);
    VectorDouble x(napices);
    VectorDouble y(napices);
    VectorDouble ysum(napices);
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int iad = iad_x;
      ysum.fill(0.);
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        VH::extractInPlace(vecin, x, iad);
        _pops[icov]->evalDirect(x, y);
        VH::multiplyConstantSelfInPlace(y, _invSills[icov].getValue(ivar, jvar));
        VH::addInPlace(ysum, y);
        iad += napices;
      }
      VH::mergeInPlace(ysum, vecout, iad_y);
      iad_y += napices;
    }
  }
  return vecout;
}

/**
 * Simulate based on an input random gaussian vector (Matrix free version)
 * @param vecin Input array
 */
VectorDouble PrecisionOpMulti::evalSimulate(const VectorDouble& vecin)
{
  int totsize = sizes();
  if ((int)vecin.size() != totsize)
  {
    messerr("'vecin' (%d) should be of dimension (%d)", (int)vecin.size(), totsize);
    return VectorDouble();
  }
  // Invert the matrices of sills
  if (_buildCholSills() != 0)
  {
    messerr("Problem when inverting the matrices of sills");
    return VectorDouble();
  }

  // Blank out the output vector
  VectorDouble vecout(totsize, 0.);

  int nvar = _getNVar();
  int ncov = _getNCov();

  int iad_x = 0;
  int iad_y = 0;
  for (int icov = 0; icov < ncov; icov++)
  {
    int napices = size(icov);
    VectorDouble x(napices);
    VectorDouble y(napices);
    VectorDouble ysum(napices);
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int iad = iad_x;
      ysum.fill(0.);
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        VH::extractInPlace(vecin, x, iad);
        _pops[icov]->evalSimulate(x, y);
          VH::multiplyConstantSelfInPlace(y, _cholSills[icov].getCholeskyTL(jvar, ivar));
        VH::addInPlace(ysum, y);
        iad += napices;
      }
      VH::mergeInPlace(ysum, vecout, iad_y);
      iad_y += napices;
    }
  }
  return vecout;
}
