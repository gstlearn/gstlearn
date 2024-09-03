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
#include "LinearOp/ALinearOp.hpp"
#include "Matrix/VectorEigen.hpp"
#include <Eigen/src/Core/Matrix.h>

#define EVALOP(prepare,start,tab,getmat,op) \
  if(_prepareOperator(vecin,vecout))\
  {\
    return 1;\
  }; \
  if (prepare() != 0)\
  {\
    messerr("Problem when preparing the matrix");\
    return 1;\
  }\
  int nvar = _getNVar();\
  int ncov = _getNCov();\
  int iad_x = 0;\
  int iad_struct = 0;\
  for (int icov = 0; icov < ncov; icov++)\
  {\
    int napices = size(icov);\
    Eigen::VectorXd y(napices);\
    for (int ivar = 0; ivar < nvar; ivar++)\
    {\
      Eigen::Map<const Eigen::VectorXd>x(vecin.data()+iad_x,napices);\
      _pops[icov]->op(x, y);\
      int iad_y = iad_struct + start * napices;\
      for (int jvar = start; jvar < nvar; jvar++)\
      {\
        VectorEigen::addMultiplyConstantInPlace(tab[icov].getmat(jvar,ivar),y,vecout,iad_y);\
        iad_y += napices;\
      }   \
      iad_x+= napices;  \
    }\
    iad_struct = iad_x;\
  }\
  return 0;


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
** Check that 'model' is a valid one (composed of Matern and Nugget Effect)
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
    if (model->getCovaType(icov) == ECov::MATERN)
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

int PrecisionOpMulti::getSize() const
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

int PrecisionOpMulti::_buildInvSills() const
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

int PrecisionOpMulti::_buildCholSills() const
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
  sstr << "Vector dimension      = " <<  getSize()  << std::endl;

  sstr << "Indices of Matérn Covariance = " << VH::toStringAsVI(_covList);

  sstr << "Dimensions of the Meshes = ";
  for (int imesh = 0, nmesh = _getNMesh(); imesh < nmesh; imesh++)
    sstr << size(imesh) << " ";
  sstr << std::endl;

  if (_isValid)
   sstr << std::endl << "Class is Valid for operations!" << std::endl;
  return sstr.str();
}

int PrecisionOpMulti::_addToDest(const Eigen::VectorXd& vecin,
                                 Eigen::VectorXd& vecout) const {
  EVALOP(_buildInvSills, 0, _invSills, getValue, addToDest)}

/* int PrecisionOpMulti::evalDirectInPlace(const Eigen::VectorXd& vecin,
                                              Eigen::VectorXd& vecout)
{
  EVALOP(_buildInvSills,0,_invSills,getValue,evalDirect)
}
 */
/**
 * Evaluate the product of this phantom matrix by the input vector
 * @param vecin Input array
 */
Eigen::VectorXd PrecisionOpMulti::evalDirect(const Eigen::VectorXd& vecin)
{
  int totsize = getSize();
  
  // Blank out the output vector
  Eigen::VectorXd vecout(totsize);
  VectorEigen::fill(vecout,0.);
  if (ALinearOp::evalDirect(vecin,vecout))
  {
    return Eigen::VectorXd();
  }

  return vecout;

}

/**
 * Simulate based on an input random gaussian vector (Matrix free version)
 * @param vecin  Input array
 * @param vecout Output array
 */

int PrecisionOpMulti::evalSimulateInPlace(const Eigen::VectorXd& vecin,
                                                Eigen::VectorXd& vecout)
{
  EVALOP(_buildCholSills,ivar,_cholSills,getCholeskyTL,evalSimulate)
}

VectorDouble PrecisionOpMulti::evalSimulate(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(),vec.size());
  Eigen::VectorXd out(vec.size());
  evalSimulateInPlace(vecm,out);
  return VectorEigen::copyIntoVD(out);
}


/**
 * Simulate based on an input random gaussian vector (Matrix free version)
 * @param vecin Input array
 */
Eigen::VectorXd PrecisionOpMulti::evalSimulate(const Eigen::VectorXd& vecin)
{
  int totsize = getSize();
  
  Eigen::VectorXd vecout(totsize);
  VectorEigen::fill(vecout,0.);
  if (evalSimulateInPlace(vecin,vecout))
  {
    return Eigen::VectorXd();
  }

  return vecout;
}

int PrecisionOpMulti::_prepareOperator(const Eigen::VectorXd& vecin,
                                             Eigen::VectorXd& vecout) const
{
   int totsize = getSize();
  if ((int)vecin.size() != totsize)
  {
    messerr("'vecin' (%d) should be of dimension (%d)", (int)vecin.size(), totsize);
    return 1;
  }
  if ((int)vecout.size() != totsize)
  {
    vecout.resize(totsize);
  }

  vecout.fill(0.);
  return 0;
}
                                            
