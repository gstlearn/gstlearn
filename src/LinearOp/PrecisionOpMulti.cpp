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
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/VectorEigen.hpp"
#include "Model/ANoStat.hpp"
#include <Eigen/src/Core/Matrix.h>

#define EVALOP(IN,OUT,TAB,getmat,OP,IY,COMPUTEOP,XORY,START,END,IVAR,JVAR) \
  {\
    int nvar = _getNVar();\
    int ncov = _getNCov();\
    int iad_x = 0;\
    int iad_struct = 0;\
    Eigen::VectorXd* y;\
    for (int icov = 0; icov < ncov; icov++)\
    {\
      int napices = size(icov);\
      if ( (nvar == 1) && (ncov == 1) )\
      {\
        y = &OUT;\
      }\
      else\
      {\
        y = nullptr;\
        if (COMPUTEOP)\
        {\
          y = &_works[icov];\
          y->resize(napices);\
          VectorEigen::fill(*y,0.);\
        }\
      }\
      for (int jvar = 0; jvar < nvar; jvar++)\
      {\
        int iad_y = IY;\
        Eigen::Map<const Eigen::VectorXd> x(IN.data() + iad_x, napices);\
        if (COMPUTEOP) \
          _pops[icov]->OP(x, *y);\
        if ( (nvar == 1) && (ncov == 1) ) break;\
        if (nvar == 1)\
        {\
          Eigen::Map<Eigen::VectorXd> outmap(OUT.data() + iad_y, napices);\
          VectorEigen::copy(*y, outmap);\
          iad_x += napices;\
          continue;\
        }\
        for (int ivar = START; ivar < END; ivar++)\
        {\
          if (_isNoStatForVariance[icov])\
          {\
            VectorEigen::addMultiplyVectVectInPlace(TAB##NoStat[icov][IND(IVAR,JVAR,nvar)],XORY,OUT,iad_y); \
          }\
          else\
          {\
            VectorEigen::addMultiplyConstantInPlace(TAB[icov].getmat(IVAR,JVAR),XORY,OUT,iad_y);\
          }\
          iad_y += napices;\
        }\
        iad_x += napices;  \
      }\
      iad_struct += napices * nvar;\
    }\
    if (COMPUTEOP) return 0;\
  }\

PrecisionOpMulti::PrecisionOpMulti(Model* model,
                                   const std::vector<const AMesh*>& meshes,
                                   bool buildOp)
  : _pops()
  , _invCholSills()
  , _cholSills()
  , _model(nullptr)
  , _meshes()
  , _isValid(false)
  , _covList() 
  , _allStat(true)
  , _ready(false)
{
  if (! _isValidModel(model)) return;

  if (! _isValidMeshes(meshes)) return;

  if (! _matchModelAndMeshes()) return;

  _isValid = true;
  _computeSize();
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
      bool nostaticov = nostat->isDefinedForVariance(icov);
      _isNoStatForVariance[icov] = nostaticov;
      _allStat = _allStat && !nostaticov;
      }
    }
  }
  _buildMatrices();

  if(buildOp)
  {
    buildQop();
  }
}

void PrecisionOpMulti::buildQop()
{
  _buildQop();
  _ready = true;
}

bool PrecisionOpMulti::_checkReady() const
{
  if (!_ready)
  {
    messerr("Operator has not been built. Computation has not been performed.");
    messerr("Call the method buildQop to make the PrecisionOpMulti ready for use");
  }
  return _ready;
}


void PrecisionOpMulti::_buildQop()
{
  for (int i = 0, number = _getNCov(); i < number; i++)
  {
    _pops.push_back(PrecisionOp::create(_meshes[i], _model, _covList[i]));
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

bool PrecisionOpMulti::_matchModelAndMeshes() const 
{
  return _getNCov() == _getNMesh();
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
  return _size;
}
void PrecisionOpMulti::_computeSize()
{
  int nvar = _getNVar();
  int ncov = _getNCov();
  _size = 0;
  for (int i = 0; i < ncov; i++)
  {
    _size += nvar * _meshes[i]->getNApices();
  }
}

int PrecisionOpMulti::_buildGlobalMatricesStationary(int icov)
{
   MatrixSquareSymmetric sill = _model->getSillValues(icov);
  if (sill.computeCholesky() != 0) return 1;
  _cholSills[icov] = sill;
   if (sill.invertCholesky() != 0) return 1;
  _invCholSills[icov] = sill;
  return 0;
}

int PrecisionOpMulti::_buildLocalMatricesNoStat(int icov)
{
  const ANoStat* nostat =  _model->getNoStat();
  int nvar = _getNVar();
  nostat->attachToMesh(_meshes[icov],false,false);
  int nvertex =  (int)_meshes[icov]->getNApices();
  int nterms = nvar * (nvar + 1)/2;
  _invCholSillsNoStat[icov].resize(nterms);
  _cholSillsNoStat[icov].resize(nterms);

  for (int i = 0; i < nterms; i++)
  {
    _invCholSillsNoStat[icov][i].resize(nvertex);
    _cholSillsNoStat[icov][i].resize(nvertex);
  }
  for (int imesh = 0; imesh < nvertex; imesh++)
  {
    _model->updateCovByMesh(imesh);
    MatrixSquareSymmetric sills = _model->getSillValues(icov);
    if (sills.computeCholesky() != 0) return 1;
    if (sills.invertCholesky()  != 0) return 1;

    int s = 0;
    for (int icol = 0; icol < nvar; icol++)
    {
      for (int irow = icol; irow < nvar; irow++)
      {
        _cholSillsNoStat[icov][s][imesh] = sills.getCholeskyTL(irow,icol);
        _invCholSillsNoStat[icov][s][imesh]  = sills.getCholeskyXL(irow,icol);
         s++;
      }
    }
  }
  return 0;
}

int PrecisionOpMulti::_buildMatrices()
{
  if (_model == nullptr) return 1;
  if (_getNVar() == 1) return 0;

  int ncov = _getNCov();

  // Do nothing if the array has already been calculated (correct dimension)
  if (ncov == (int)_cholSills.size()) return 0;
  _cholSills.resize(ncov);
  _invCholSills.resize(ncov);
  _cholSillsNoStat.resize(ncov);
  _invCholSillsNoStat.resize(ncov);
 
  for (int icov = 0; icov < ncov; icov++)
  {
   if (_isNoStatForVariance[icov])
   {  
    if (_buildLocalMatricesNoStat(icov)) return 1;
   }
   else
   {
    if (_buildGlobalMatricesStationary(icov)) return 1;
   }
  }
  return 0;
}

void PrecisionOpMulti::makeReady()
{
  _makeReady();
}

void PrecisionOpMulti::_makeReady()
{
  for (auto &e : _pops)
    e->makeReady();
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

int PrecisionOpMulti::_addToDestImpl(const Eigen::VectorXd& vecin,
                          Eigen::VectorXd& vecout) const
{
  if (!_checkReady()) return 1;
  if (_getNVar() > 1)
  {
    _workTot.resize(vecin.size());
    VectorEigen::fill(_workTot,0.);
  
   EVALOP(vecin,_workTot ,_invCholSills,getCholeskyXL,addToDest,iad_struct + jvar * napices,false,x,jvar,nvar,ivar,jvar)
   EVALOP(_workTot, vecout ,_invCholSills,getCholeskyXL,addToDest,iad_struct,true,*y,0,(jvar+1),jvar,ivar)
  }
  else 
  {
    EVALOP(vecin, vecout ,_invCholSills,getCholeskyXL,addToDest,iad_struct,true,*y,0,1,0,0)
  }
}

int PrecisionOpMulti::_addToDest(const Eigen::VectorXd& vecin,
                          Eigen::VectorXd& vecout) const
{
  return _addToDestImpl(vecin,vecout);
  
}
/**
 * Simulate based on an input random gaussian vector (Matrix free version)
 * @param vecin  Input array
 * @param vecout Output array
 */

int PrecisionOpMulti::_addSimulateInPlace(const Eigen::VectorXd& vecin,
                                                Eigen::VectorXd& vecout)
{
  if (!_checkReady()) return 1;
  EVALOP(vecin,vecout,_cholSills,getCholeskyTL,evalSimulate,iad_struct + jvar * napices,true,*y,jvar,nvar,ivar,jvar)
}

VectorDouble PrecisionOpMulti::evalSimulate(const VectorDouble& vec)
{
  if (!_checkReady()) return VectorDouble(); 
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(),vec.size());
  Eigen::VectorXd out(vec.size());
  VectorEigen::fill(out,0.);
  _addSimulateInPlace(vecm,out);
  return VectorEigen::copyIntoVD(out);
}

                                            
