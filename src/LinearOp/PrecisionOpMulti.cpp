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
#include "Model/ANoStat.hpp"
#include <Eigen/src/Core/Matrix.h>

#define EVALOP(prepare,start,tab,getmat,op) \
    int nvar = _getNVar();\
    if (nvar > 1)\
    {\
      if(_prepareOperator(vecin,vecout))\
      {\
        return 1;\
      }; \
      if (prepare() != 0)\
      {\
        messerr("Problem when preparing the matrix");\
        return 1;\
      }\
    }\
    int ncov = _getNCov();\
    int iad_x = 0;\
    int iad_struct = 0;\
    Eigen::VectorXd work;\
    Eigen::VectorXd* y;\
    for (int icov = 0; icov < ncov; icov++)\
    {\
      int napices = size(icov);\
      if (nvar == 1 && ncov == 1)\
      {\
        y = &vecout;\
      }\
      else\
      {\
        y = &_works[icov];\
        y->resize(napices);\
      }\
      int s = 0;\
      for (int jvar = 0; jvar < nvar; jvar++)\
      {\
        Eigen::Map<const Eigen::VectorXd> x(vecin.data() + iad_x, napices);\
        _pops[icov]->op(x, *y);\
        if (nvar == 1 && ncov == 1) break;\
        int iad_y = iad_struct + start * napices;\
        if (nvar == 1)\
        {\
          Eigen::Map<Eigen::VectorXd> outmap(vecout.data() + iad_y, napices);\
          VectorEigen::copy(*y, outmap);\
          iad_x += napices;\
          continue;\
        }\
        for (int ivar = start; ivar < nvar; ivar++)\
        {\
          if (_isNoStatForVariance[icov])\
          {\
            VectorEigen::addMultiplyVectVectInPlace(tab##NoStat[icov][s++],*y,vecout,iad_y);\
          }\
          else\
          {\
            VectorEigen::addMultiplyConstantInPlace(tab[icov].getmat(ivar,jvar),*y,vecout,iad_y);\
          }\
          iad_y += napices;\
        }\
        iad_x+= napices;  \
      }\
      iad_struct = iad_x;\
    }\
    return 0;


PrecisionOpMulti::PrecisionOpMulti(Model* model,
                                   const std::vector<const AMesh*>& meshes)
  : _pops()
  , _invSills()
  , _cholSills()
  , _model(nullptr)
  , _meshes()
  , _isValid(false)
  , _covList()
{
  if (! _isValidModel(model)) return;

  if (! _isValidMeshes(meshes)) return;

  if (! _matchModelAndMeshes()) return;

  _isValid = true;

  int ncov = (int)meshes.size(); 
  _isNoStatForVariance.resize(ncov,false);
  
  _works.resize(ncov);

  bool isnostat = _model->isNoStat();

  if (isnostat)
  {
    const ANoStat* nostat =  _model->getNoStat();

    for (int icov = 0; icov < ncov; icov++)
    {
      if (nostat != nullptr)
      {
      _isNoStatForVariance[icov] = nostat->isDefinedForVariance(icov);
      }
    }
  }
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
bool PrecisionOpMulti::_isValidMeshes(const std::vector<const AMesh*>& meshes)
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
    _pops.push_back(PrecisionOp::create(_meshes[i], _model, _covList[i]));
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

int PrecisionOpMulti::_buildInvSills() const
{
  if (_model == nullptr) return 1;

  if (_getNVar() == 1) return 0;
  int ncov = _getNCov();

  // Do nothing if the array has already been calculated (correct dimension)
  if (ncov == (int) _invSills.size()) return 0;

  _invSills.resize(ncov);
  _invSillsNoStat.resize(ncov);
  for (int icov = 0; icov < ncov; icov++)
  {
    MatrixSquareSymmetric invs = _model->getSillValues(icov);
    if (invs.invert() != 0) return 1;
    _invSills[icov] = invs;
  }
  return 0;
}

int PrecisionOpMulti::_buildCholSills() const
{
  if (_getNVar() == 1) return 0;
  if (_model == nullptr) return 1;

  int ncov = _getNCov();

  // Do nothing if the array has already been calculated (correct dimension)
  if (ncov == (int)_cholSills.size()) return 0;

  _cholSills.resize(ncov);
  _cholSillsNoStat.resize(ncov);

  MatrixSquareSymmetric chols;
  int nvar = _getNVar();

  const ANoStat* nostat =  _model->getNoStat();
 
  for (int icov = 0; icov < ncov; icov++)
  {
   
   if (_isNoStatForVariance[icov])
   {  
    nostat->attachToMesh(_meshes[icov],false,false);
    int nvertex =  (int)_meshes[icov]->getNApices();
    int nterms = nvar * (nvar + 1)/2;
    _cholSillsNoStat[icov].resize(nterms);

    for (int i = 0; i < nterms; i++)
    {
      _cholSillsNoStat[icov][i].resize(nvertex);
    }

    for (int imesh = 0; imesh < nvertex; imesh++)
    {
      _model->updateCovByMesh(imesh);
      chols = _model->getSillValues(icov);
      if (chols.computeCholesky() != 0) return 1;
      int s = 0;
      for (int icol = 0; icol < nvar; icol++)
      {
        for (int irow = icol; irow < nvar; irow++)
        {
          _cholSillsNoStat[icov][s][imesh] = chols.getCholeskyTL(irow,icol);
           s++;
        }
      }
    }
   }
   else
   {
      chols = _model->getSillValues(icov);
      if (chols.computeCholesky() != 0) return 1;
      _cholSills[icov] = chols;
   }
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
                          Eigen::VectorXd& vecout) const
{
  EVALOP(_buildInvSills,0,_invSills,getValue,addToDest)
}

/**
 * Evaluate the product of this phantom matrix by the input vector
 * @param vecin  Input array
 * @param vecout Output array
 */

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
  EVALOP(_buildCholSills,jvar,_cholSills,getCholeskyTL,evalSimulate)
 /*  int nvar = _getNVar();
  if (nvar > 1)
  {
    if (_prepareOperator(vecin, vecout))
    {
      return 1;
    };
    if (_buildCholSills() != 0)
    {
      messerr("Problem when preparing the matrix");
      return 1;
    }
  }
  int ncov       = _getNCov();
  int iad_x      = 0;
  int iad_struct = 0;
  Eigen ::VectorXd* y;
  for (int icov = 0; icov < ncov; icov++)
  {
    int napices = size(icov);
    if (nvar == 1 && ncov == 1)
    {
      y = &vecout;
    }
    else
    {
      y = &(_works[icov]);
      y->resize(napices);
    }
    int s = 0;
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      Eigen ::Map<const Eigen ::VectorXd> x(vecin.data() + iad_x, napices);
      _pops[icov]->evalSimulate(x, *y);
      if (nvar == 1 && ncov == 1) break;
      int iad_y = iad_struct + jvar * napices;
      if (nvar == 1)
      {
        std ::cout << iad_x << " " << napices << std ::endl;
        Eigen ::Map<Eigen ::VectorXd> outmap(vecout.data() + iad_x, napices);
        iad_x += napices;
        VectorEigen ::copy(*y, outmap);
        continue;
      }
      std ::cout << "not here" << std ::endl;
      for (int ivar = jvar; ivar < nvar; ivar++)
      {
        if (_isNoStatForVariance[icov])
        {
          VectorEigen ::addMultiplyVectVectInPlace(_cholSillsNoStat[icov][s++],
                                                 *y, vecout, iad_y);
        }
        else
        {
          VectorEigen ::addMultiplyConstantInPlace(
          _cholSills[icov].getCholeskyTL(ivar, jvar), *y, vecout, iad_y);
        }
        iad_y += napices;
      }
      iad_x += napices;
    }
    iad_struct = iad_x;
  }
  return 0; */
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
                                            
