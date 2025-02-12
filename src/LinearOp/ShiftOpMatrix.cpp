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
#include "LinearOp/ShiftOpMatrix.hpp"

#include "Enum/EConsElem.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovAniso.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Space/SpaceSN.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"

#include <math.h>
#include <memory>

ShiftOpMatrix::ShiftOpMatrix()
  : AShiftOp()
  , _TildeC()
  , _S(nullptr)
  , _nCovAnisoGradParam(0)
  , _SGrad()
  , _TildeCGrad()
  , _LambdaGrad()
  , _flagNoStatByHH(false)
  , _ndim(0)
{
}

ShiftOpMatrix::ShiftOpMatrix(const AMesh* amesh,
                             const CovAniso* cova,
                             const Db* dbout,
                             bool verbose)
  : AShiftOp()
  , _TildeC()
  , _S(nullptr)
  , _nCovAnisoGradParam(0)
  , _SGrad()
  , _TildeCGrad()
  , _LambdaGrad()
  , _flagNoStatByHH(false)
  , _ndim(amesh->getEmbeddedNDim())
{
  if (amesh != nullptr)
    _napices = amesh->getNApices();
  (void)initFromMesh(amesh, cova, dbout, verbose);
}

ShiftOpMatrix::ShiftOpMatrix(const MatrixSparse* S,
                             const VectorDouble& TildeC,
                             const VectorDouble& Lambda,
                             const CovAniso* cova,
                             bool verbose)
  : AShiftOp()
  , _TildeC()
  , _S(nullptr)
  , _nCovAnisoGradParam(0)
  , _SGrad()
  , _TildeCGrad()
  , _LambdaGrad()
  , _flagNoStatByHH(false)
  , _ndim(0)
{
  if (S != nullptr)
    _napices = S->getNCols();
  (void)initFromCS(S, TildeC, Lambda, cova, verbose);
  (void)initFromCS(S, TildeC, Lambda, cova, verbose);
}

ShiftOpMatrix::ShiftOpMatrix(const ShiftOpMatrix& shift)
  : AShiftOp(shift)
  , _TildeC()
  , _S(nullptr)
  , _nCovAnisoGradParam(0)
  , _SGrad()
  , _TildeCGrad()
  , _LambdaGrad()
  , _flagNoStatByHH(false)
  , _ndim(0)
{
  _reallocate(shift);
}

ShiftOpMatrix& ShiftOpMatrix::operator=(const ShiftOpMatrix &shift)
{
  if (this != &shift)
  {
    _reset();
    _reallocate(shift);
  }
  return *this;
}

ShiftOpMatrix::~ShiftOpMatrix()
{
  _reset();
}

ShiftOpMatrix* ShiftOpMatrix::create(const AMesh *amesh,
                             const CovAniso *cova,
                             const Db *dbout,
                             bool verbose)
{
  return new ShiftOpMatrix(amesh, cova, dbout, verbose);
}

ShiftOpMatrix* ShiftOpMatrix::createFromSparse(const MatrixSparse *S,
                                       const VectorDouble &TildeC,
                                       const VectorDouble &Lambda,
                                       const CovAniso *cova,
                                       bool verbose)
{
  return new ShiftOpMatrix(S, TildeC, Lambda, cova, verbose);
}

/**
 *
 * @param amesh Meshing description (New format)
 * @param cova Pointer to the CovAniso structure
 * @param dbout Pointer to the Db structure
 * @param flagAdvection When TRUE, S is replaced by G
 * @param verbose Verbose flag
 * @return Error return code
 */
int ShiftOpMatrix::initFromMesh(const AMesh* amesh,
                            const CovAniso* cova,
                            const Db* /*dbout*/,
                            bool flagAdvection,
                            bool verbose)
{
  DECLARE_UNUSED(flagAdvection,verbose);

  // Initializations

  _setCovAniso(cova);
  _napices = amesh->getNApices();
  try
  {
    _ndim = amesh->getEmbeddedNDim();

    // Attach the Non-stationary to Mesh and Db (optional)

    _cova->informMeshByMeshForAnisotropy(amesh);

    // Calculating and storing the mesh sizes
    VectorDouble units = amesh->getMeshSizes();

    // Define if parameterization is in HH or in range/angle
    _determineFlagNoStatByHH();

    // Identify the covariance

    // Construct S sparse Matrix
    if (_buildS(amesh))
      my_throw("Problem when buildS");

    // Construct the Lambda vector

    _buildLambda(amesh);
  }

  catch(const AException& e)
  {
    messerr("initFromMesh has failed: %s",e.what());
    _reset();
    return 1;
  }
  catch(const std::exception& e)
  {
    messerr("initFromMesh has failed: %s",e.what());
    _reset();
    return 1;
  }
  return 0;
}

/**
 * Initialize the environment for calculation of derivatives of S
 * @param amesh Meshing description (New format)
 * @param cova Pointer to the CovAniso structure
 * @param verbose Verbose flag
 * @param tol Smallest value below which the value is not stored in sparse matrix
 * @return Error return code
 */
int ShiftOpMatrix::initGradFromMesh(const AMesh* amesh,
                                const CovAniso* cova,
                                bool verbose,
                                double tol)
{
  // Initializations
  DECLARE_UNUSED(verbose)
  _setCovAniso(cova);

  try
  {
    _cova->informMeshByMesh(amesh);

    // Construct S sparse Matrix
    if (_buildSGrad(amesh, tol))
      my_throw("Problem when buildSGrad");

    if (_buildLambdaGrad(amesh))
      my_throw("Problem when buildLambdaGrad");
  }

  catch(const AException& e)
  {
    messerr("initGradFromMesh has failed: %s",e.what());
    _reset();
    return 1;
  }
  catch(const std::exception& e)
  {
    messerr("initGradFromMesh has failed: %s",e.what());
    _reset();
    return 1;
  }

  return 0;
}

/**
 *
 * @param S Sparse matrix describing the S information
 * @param TildeC Diagonal array containing TildeC
 * @param Lambda Normalization vector
 * @param cova Pointer to the CovAniso structure
 * @param verbose Verbose flag
 * @return
 */
int ShiftOpMatrix::initFromCS(const MatrixSparse* S,
                              const VectorDouble& TildeC,
                              const VectorDouble& Lambda,
                              const CovAniso* cova,
                              bool verbose)
{
  DECLARE_UNUSED(verbose);

  // Initializations

  _setCovAniso(cova);

  try
  {
    // Store the TildeC & Lambda vectors

    _TildeC = TildeC;
    _Lambda = Lambda;
    _ndim   = cova->getNDim();

    // Duplicate the Shift Operator sparse matrix

    _S = S->clone();
    _napices = S->getNCols();
  }

  catch(const AException& e)
  {
    messerr("initFromCS has failed: %s",e.what());
    return 1;
  }
  catch(const std::exception& e)
  {
    messerr("initFromCS has failed: %s",e.what());
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Perform the operation: y = x * C^power
 **
 ** \param[in] x       Input vector
 ** \param[in] y       Output vector
 ** \param[in] power   Value of the exponent
 **
 ** \remarks 'C' is a member (_TildeC) that stands as a vector
 ** \remarks Specific coding has been realized for the cases
 ** \remarks where 'power' is equal to 1, -1, 0.5 and -0.5
 **
 *****************************************************************************/
void ShiftOpMatrix::prodTildeC(const VectorDouble& x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
  if (power == EPowerPT::ONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * _TildeC[i];
  }
  else if (power == EPowerPT::MINUSONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / _TildeC[i];
  }
  else if (power == EPowerPT::HALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * sqrt(_TildeC[i]);
  }
  else if (power == EPowerPT::MINUSHALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / sqrt(_TildeC[i]);
  }
  else if (power == EPowerPT::LOG)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i];
  }
  else
  {
    my_throw("Unexpected value for argument 'power'");
  }
}

void ShiftOpMatrix::normalizeLambdaBySills(const AMesh* mesh)
{
  VectorDouble tab;

  if (_cova->isNoStatForVariance())
  {
    _cova->informMeshByApexForSills(mesh);
    int number = (int) _Lambda.size();
                       
    
    for (int imesh = 0; imesh < number; imesh++)
    {
      _cova->updateCovByMesh(imesh,false);
      double sill = _cova->getSill(0,0);
      double invsillsq = 1. / sqrt(sill);
      _Lambda[imesh] *= invsillsq;
    }
  }
  else 
  {
    double invsillsq = 1. / sqrt(_getCovAniso()->getSill(0,0));
    for (auto &e:_Lambda)
    {
      e *= invsillsq;
    }
  }
}





void ShiftOpMatrix::prodLambdaOnSqrtTildeC(const VectorDouble& inv,
                                     VectorDouble& outv,
                                     double puis) const
{
  for (int i = 0; i < getSize(); i++)
    outv[i] = inv[i] * pow(_Lambda[i] / sqrt(_TildeC[i]), puis);
}

/****************************************************************************/
/*!
 **  Perform the operation: y = S * x
 **
 ** \param[in] inv     Input vector
 ** \param[in] outv    Output vector
 **
 ** \remarks 'S' is a member that stands as a sparse matrix
 **
 *****************************************************************************/
int ShiftOpMatrix::_addToDest(const constvect inv, vect outv) const
{
  _S->addProdMatVecInPlaceToDest(inv, outv);
  return 0;
}

int ShiftOpMatrix::_resetGrad()
{
  if (_SGrad.empty()) return 1;
  for (int i = 0; i < (int) _SGrad.size(); i++)
  {
     delete _SGrad[i];
     _SGrad[i] = nullptr;
  }
  return 0;
}

void ShiftOpMatrix::_reset()
{
  delete _S;
  _S = nullptr;
  _resetGrad();
}

void ShiftOpMatrix::_reallocate(const ShiftOpMatrix& shift)
{
  _TildeC = shift._TildeC;
  _Lambda = shift._Lambda;
  _S = shift._S->clone();
  _nCovAnisoGradParam = shift._nCovAnisoGradParam;
  for (int i = 0; i < (int) _SGrad.size(); i++)
    _SGrad[i] = shift._SGrad[i]->clone();
  for (int i = 0; i < (int) _LambdaGrad.size(); i++)
    _LambdaGrad[i] = shift._LambdaGrad[i];
  _flagNoStatByHH = shift._flagNoStatByHH;

  _cova = shift._cova;
  _ndim = shift._ndim;
  _napices = shift._napices;
}

MatrixSparse* ShiftOpMatrix::getTildeCGrad(int iapex, int igparam) const
{

  if (_TildeCGrad.empty())
    {
      messerr("You must initialize the Gradients with 'initGradFromMesh' beforehand");
      return nullptr;
    }
    int iad = getSGradAddress(iapex, igparam);
    if (iad < 0) return nullptr;

    return _TildeCGrad[iad];
}

MatrixSparse* ShiftOpMatrix::getSGrad(int iapex, int igparam) const
{
  if (_SGrad.empty())
  {
    messerr("You must initialize the Gradients with 'initGradFromMesh' beforehand");
    return nullptr;
  }
  int iad = getSGradAddress(iapex, igparam);
  if (iad < 0) return nullptr;

  return _SGrad[iad];
}

/**
 * Locally update the covariance (does nothing in Stationary Case)
 * @param cova Local CovAniso structure (updated here)
 * @param imesh Rank of the active mesh
 */
void ShiftOpMatrix::_updateCova(std::shared_ptr<CovAniso> &cova, int imesh)
{
  cova->updateCovByMesh(imesh,true);

}

void ShiftOpMatrix::_updateHH(MatrixSquareSymmetric& hh, int imesh)
{
  DECLARE_UNUSED(imesh)
  if (! _isNoStat()) return;
  int ndim = getNDim();

  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = idim; jdim < ndim; jdim++)
    {
      double value = 0.; //TODO repare
    //  double value = nostat->getValue(EConsElem::TENSOR, 0, imesh, idim, jdim);
      hh.setValue(idim, jdim, value);
    }
}

/**
 * Calculate HH matrix from parameters
 * Note that this function is also called in Stationary case...
 * So it must be workable without any updates
 * @param hh Output Array
 * @param imesh Rank of the active mesh
 */
void ShiftOpMatrix::_loadHHRegular(MatrixSquareSymmetric &hh, int imesh)
{
  int ndim = getNDim();
  // Locally update the covariance for non-stationarity (if necessary)
  _updateCova(_getCovAniso(), imesh);

  // Calculate the current HH matrix (using local covariance parameters)
  const MatrixSquareGeneral &rotmat = _getCovAniso()->getAnisoInvMat();

  VectorDouble diag = VH::power(_getCovAniso()->getScales(), 2.);
  MatrixSquareSymmetric temp(ndim);
  temp.setDiagonal(diag);
  hh.normMatrix(rotmat, temp);

}

void ShiftOpMatrix::_loadHHVariety(MatrixSquareSymmetric& hh, int imesh)
{
  int ndim = getNDim();

  // Locally update the covariance for non-stationarity (if necessary)
  _updateCova(_getCovAniso(), imesh);

  // Calculate the current HH matrix (using local covariance parameters)
  VectorDouble diag = VH::power(_getCovAniso()->getScales(), 2.);

  hh.fill(0.);
  for (int idim = 0; idim < ndim; idim++)
    hh.setValue(idim, idim, diag[0]);
}

/**
 * Calculate HH Gradient matrix from one of the CovAniso parameters for the given Apex.
 * @param amesh Meshing description (New format)
 * @param hh         Output Array (updated here)
 * @param igparam    Rank of the parameter for derivation
 * @param ipref      Rank of the point
 *
 * @details: The parameters 'igparam' are sorted as follows:
 * @details: - 0:(ndim-1)   : ranges in each Space direction
 * @details: - ndim:ngparam : rotation angles (=ndim or 1 in 2-D)
 */

void ShiftOpMatrix::_loadHHGrad(const AMesh *amesh,
                            MatrixSquareSymmetric &hh,
                            int igparam,
                            int ipref)
{
  int ndim = getNDim();

  if (_flagNoStatByHH)
  {
    // Case where the derivation must be performed on the HH terms

    hh.fill(0.);
    int ecr = 0;
    for (int idim = 0; idim < ndim; idim++)
      for (int jdim = idim; jdim < ndim; jdim++)
      {
        if (ecr == igparam) hh.setValue(idim, jdim, 1.);
        ecr++;
      }
  }
  else
  {
    // Case where the derivation is performed on ranges and angles


    // Locally update the covariance for non-stationarity (if necessary)

    _updateCova(_getCovAniso(), ipref);
    const MatrixSquareGeneral &rotmat = _getCovAniso()->getAnisoRotMat();
    VectorDouble diag = VH::power(_getCovAniso()->getScales(), 2.);

    MatrixSquareSymmetric temp(ndim);
    if (igparam < ndim)
    {
      // Derivation with respect to the Range 'igparam'
      temp.fill(0.);
      temp.setValue(igparam, igparam, 2. * _getCovAniso()->getScale(igparam));
      hh.normMatrix(rotmat, temp);
    }
    else
    {
      // Derivation with respect to the Angle 'igparam'-ndim
      int ir = igparam - ndim;
      auto covini = _getCovAniso();
      auto covaderiv = cloneAndCast(covini);
      _updateCova(covaderiv, ipref);
      _getCovAniso()->setAnisoAngle(ir, covaderiv->getAnisoAngles(ir) + 90.);
      const MatrixSquareGeneral &drotmat = covaderiv->getAnisoRotMat();

      VH::divideConstant(diag, 180. / GV_PI); // Necessary as angles are provided in degrees. Factor 2 is for derivative
      temp.setDiagonal(diag);
      hh.innerMatrix(temp, drotmat, rotmat);
    }
  }

  int number = amesh->getNApexPerMesh();
  hh.prodScalar(1. / number);
}

double ShiftOpMatrix::_computeGradLogDetHH(const AMesh *amesh,
                                       int igparam,
                                       int ipref,
                                       const MatrixSquareSymmetric &invHH,
                                       MatrixSquareSymmetric &work,
                                       MatrixSquareSymmetric &work2)
{
  int ndim = getNDim();
  int number = amesh->getNApexPerMesh();

  if (igparam < ndim)
  {
    _updateCova(_getCovAniso(), ipref);
    const MatrixSquareGeneral &rotmat = _getCovAniso()->getAnisoRotMat();
    MatrixSquareSymmetric temp(ndim);
    temp.setDiagonal(_getCovAniso()->getScales());
    for (int idim = 0; idim < ndim; idim++)
    {
      if (idim != igparam)
      {
        temp.setValue(idim, idim, 0.);
      }
      else
      {
        temp.setValue(idim, idim, 2. * _getCovAniso()->getScale(idim) / number);
      }
    }

    work.normMatrix(rotmat, temp);
    work2.prodMatMatInPlace(&work, &invHH);
    double result = work2.trace();
    return result;
  }
  return 0.;
}

/**
 * Constitute HH (only in non-stationary case)
 * @param hh   Returned HH symmetric matrix
 * @param amesh Pointer to the meshing
 * @param imesh Rank of the mesh
 */
void ShiftOpMatrix::_loadHH(const AMesh *amesh,
                        MatrixSquareSymmetric &hh,
                        int imesh)
{
  if (_flagNoStatByHH)
  {
    _updateHH(hh, imesh);
  }
  else
  {
    if (amesh->getVariety() == 0)
      _loadHHRegular(hh, imesh);
    else
      _loadHHVariety(hh, imesh);
  }
}

void ShiftOpMatrix::_loadAux(VectorDouble &tab, const EConsElem &type, int imesh)
{
  DECLARE_UNUSED(tab,type,imesh)
  // TODO Repare
  // if (! _isNoStat()) return;
  // if (tab.empty()) return;
  // int size = static_cast<int> (tab.size());

  // VectorDouble tabloc(size, 0.);
  // for (int i = 0; i < size; i++) tab[i] = 0;

  // const ANoStatCov* nostat = _getCovAniso()->getNoStat();
  // for (int i = 0; i < (int) tab.size(); i++)
  //   if (nostat->isDefined(type, i, -1))
  //     tab[i] = nostat->getValue(type, 0, imesh, i, -1);
}

int ShiftOpMatrix::_preparMatrices(const AMesh *amesh,
                               int imesh,
                               MatrixSquareGeneral& matu,
                               MatrixRectangular& matw) const
{
  int ndim = _ndim;
  int ncorner = amesh->getNApexPerMesh();

  for (int icorn = 0; icorn < ncorner; icorn++)
  {
    for (int idim = 0; idim < ndim; idim++)
      matu.setValue(idim, icorn, amesh->getCoor(imesh, icorn, idim));
    matu.setValue(ncorner - 1, icorn, 1.);
  }

  if (matu.invert())
  {
    messerr("Problem for Mesh #%d", imesh + 1);
    amesh->printMesh(imesh);
    return 1;
  }

  for (int icorn = 0; icorn < ncorner; icorn++)
    for (int idim = 0; idim < ndim; idim++)
      matw.setValue(idim, icorn, matu.getValue(icorn, idim));

  return 0;
}

MatrixSparse* ShiftOpMatrix::_BuildTildeCGradfromMap(std::map< int, double> &tab) const
{
  std::map<int, double>::iterator it;

  NF_Triplet NF_T;
  int ip1_max = -1;

  it = tab.begin();
  while (it != tab.end())
  {
    int ip1 = it->first;
    NF_T.add(ip1, ip1, it->second);
    if (ip1 > ip1_max) ip1_max = ip1;
    it++;
  }

  // Add the fictitious value at maximum sparse matrix dimension

  NF_T.force(getSize(), getSize());

  /* Optional printout */

  return MatrixSparse::createFromTriplet(NF_T, getSize(), getSize(), -1);
}

int ShiftOpMatrix::_prepareMatricesSVariety(const AMesh *amesh,
                                        int imesh,
                                        VectorVectorDouble &coords,
                                        MatrixRectangular& matM,
                                        MatrixSquareSymmetric &matMtM,
                                        AMatrix &matP,
                                        double *deter) const
{
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();

  amesh->getEmbeddedCoordinatesPerMeshInPlace(imesh, coords);

  for (int icorn = 0; icorn < ncorner - 1; icorn++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      double val = coords[icorn][idim] - coords[ncorner - 1][idim];
      matM.setValue(idim, icorn, val);
    }
  }

  // Calculate M^t %*% M

  matMtM.normMatrix(matM);
  *deter = matMtM.determinant();

  // Calculate (M^t %*% M)^{-1}

  if (matMtM.invert())
  {
    messerr("Problem for Mesh #%d", imesh + 1);
    amesh->printMesh(imesh);
    return 1;
  }

  // Calculate P = (M^t %*% M)^{-1} %*% M^t
  matP.prodMatMatInPlace(&matMtM, &matM, false, true);
  return 0;
}

int ShiftOpMatrix::_prepareMatricesSphere(const AMesh *amesh,
                                      int imesh,
                                      VectorVectorDouble &coords,
                                      AMatrixSquare &matMs,
                                      double *deter) const
{
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();

  amesh->getEmbeddedCoordinatesPerMeshInPlace(imesh, coords);

  for (int icorn = 0; icorn < ncorner - 1; icorn++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      double val = coords[icorn][idim] - coords[ncorner - 1][idim];
      matMs.setValue(idim, icorn, val);
    }
  }

  double detM = matMs.determinant();
  *deter = detM * detM;
  if (matMs.invert())
  {
    messerr("Problem for Mesh #%d", imesh + 1);
    amesh->printMesh(imesh);
    return 1;
  }
  return 0;
}

/**
 * Calculate the private member "_S" directly from the Mesh on a Variety in 3-D space
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 *
 * @remark TildeC is calculated at the same time
 */
int ShiftOpMatrix::_buildS(const AMesh *amesh, double tol)
{
  DECLARE_UNUSED(tol);
  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int napices = amesh->getNApices();
  int nmeshes = amesh->getNMeshes();

  _TildeC.clear();
  _TildeC.resize(napices, 0.);

  // Initialize the arrays

  VectorDouble srot(2);
  MatrixSquareSymmetric hh(ndim);
  MatrixSquareSymmetric matMtM(ncorner-1);
  MatrixRectangular matP(ncorner-1,ndim);
  MatrixRectangular matM(ndim, ncorner - 1);
  MatrixSquareGeneral matMs(ndim);
  MatrixSquareSymmetric matPinvHPt(ncorner-1);
  VectorVectorDouble coords = amesh->getEmbeddedCoordinatesPerMesh();
  double detMtM = 0.;
  double dethh = 0.;

  // Define the global matrices


  if (_isGlobalHH())
  {
    _loadHH(amesh, hh, 0);
    dethh = 1. / hh.determinant();
  }
  
  //TODO : repare shpere
  /* if (! _isNoStat())
    _loadAux(srot, EConsElem::SPHEROT, 0); */

  _S = _prepareSparse(amesh);
  if (_S == nullptr) goto label_end;

  /* Loop on the active meshes */
  
  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    OptDbg::setCurrentIndex(imesh + 1);

    // Non stationary case

    if (_cova->isNoStatForAnisotropy())
    {
        _loadHH(amesh, hh, imesh);
        dethh = 1. / hh.determinant();
    }
     // if (nostat->isDefined(EConsElem::SPHEROT, -1, -1))
     //   _loadAux(srot, EConsElem::SPHEROT, imesh);

    // Prepare M matrix

    bool flagSphere = (amesh->getVariety() != 0);
    flagSphere = false; // Modif DR a valider
    if (! flagSphere)
    {
      if (_prepareMatricesSVariety(amesh, imesh, coords, matM, matMtM, matP, &detMtM))
        my_throw("Matrix inversion");
      matPinvHPt.normMatrix(matP, hh, true);
    }
    else
    {
      if (_prepareMatricesSphere(amesh, imesh, coords, matMs, &detMtM))
        my_throw("Matrix inversion");
      matPinvHPt.normMatrix(matMs, hh, true);
    }

    // Storing in the Element of the Sparse Matrix

    double ratio = sqrt(dethh * detMtM);
    double S = 0.;
    for (int j0 = 0; j0 < ncorner-1; j0++)
    {
      int ip0 = amesh->getApex(imesh, j0);
      _TildeC[ip0] += ratio / 6.;

      double s = 0.;
      for (int j1 = 0; j1 < ncorner-1; j1++)
      {
        int ip1 = amesh->getApex(imesh, j1);
        double vald = matPinvHPt.getValue(j0, j1) * ratio / 2.;
        s += vald;
        _S->addValue(ip0, ip1, vald);
      }
      int ip1 = amesh->getApex(imesh, ncorner - 1);
      _S->addValue(ip0, ip1, -s);
      _S->addValue(ip1, ip0, -s);
      S += s;
    }
    int ip0 = amesh->getApex(imesh, ncorner - 1);
    _TildeC[ip0] += ratio / 6.;
    _S->addValue(ip0, ip0, S);
  }

  // Ending S construction

  _S->prodNormDiagVecInPlace(_TildeC, -3);

  /* Set the error return code */

  error = 0;

  label_end:
  if (error)
  {
    delete _S;
    _S = nullptr;
  }
  return error;
}

MatrixSparse* ShiftOpMatrix::_prepareSparse(const AMesh *amesh)
{
  MatrixSparse* Sret = nullptr;
  MatrixSparse* Sl = nullptr;
  int nmeshes = amesh->getNMeshes();
  int ncorner = amesh->getNApexPerMesh();

  // Define Sl as the sparse matrix giving the clutter of apices among vertices
  NF_Triplet NF_T;
  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    for (int ic = 0; ic < ncorner; ic++)
    {
      int iapex = amesh->getApex(imesh, ic);
      NF_T.add(iapex, imesh, 1.);
    }
  }
  Sl = MatrixSparse::createFromTriplet(NF_T);

  // Operate the product Sl * t(Sl) to get the final matrix Sret
  Sret = prodNormMat(Sl, VectorDouble(), false);
  delete Sl;

  // Blank out the contents of the sparse matrix
  Sret->setConstant(0.);

  return Sret;
}

bool ShiftOpMatrix::_cond(int indref, int igparam, int ipref)
{
  return ipref == indref && igparam == 0;
}

/**
 * Calculate the private member "_SGrad" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpMatrix::_buildSGrad(const AMesh *amesh, double tol)
{
  auto cova = _getCovAniso();
  _nCovAnisoGradParam = cova->getNGradParam();
  int number = _nCovAnisoGradParam * getSize();
  VectorT<std::map<int, double> > tab(number);
  std::vector<std::map<std::pair<int, int>, double> > Mtab(number);

  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int ngparam = _nCovAnisoGradParam;

  // Initialize the arrays

  MatrixSquareSymmetric hh(ndim);
  MatrixSquareSymmetric work(ndim);
  MatrixSquareSymmetric work2(ndim);
  MatrixSquareSymmetric hhGrad(ndim);
  MatrixRectangular matM(ndim, ncorner - 1);
  MatrixSquareSymmetric matMtM(ncorner-1);
  MatrixRectangular matP(ncorner-1,ndim);
  MatrixSquareGeneral matMs(ndim);
  MatrixSquareSymmetric matPGradHPt(ncorner-1);
  MatrixSquareSymmetric matPHHPt(ncorner-1);

  double detMtM = 0.;
  double dethh = 0.;
  double gradLogDetHH = 0.;
  // Define the global matrices

  if (_isGlobalHH())
    _loadHH(amesh, hhGrad, 0);

  /* Loop on the meshes */

  VectorVectorDouble coords = amesh->getEmbeddedCoordinatesPerMesh();
  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setCurrentIndex(imesh + 1);

    // Prepare M matrix
    if (_cova->isNoStatForAnisotropy())
      _loadHH(amesh, hh, imesh);
    bool flagSphere = (amesh->getVariety() == 1);
    flagSphere = false; // Modif DR a valider
    if (! flagSphere)
    {
      if (_prepareMatricesSVariety(amesh, imesh, coords, matM, matMtM, matP, &detMtM))
        my_throw("Matrix inversion");
      matPHHPt.normMatrix(matP, hh, true);
    }

    else
    {
      if (_prepareMatricesSphere(amesh, imesh, coords, matMs, &detMtM))
        my_throw("Matrix inversion");
      matPHHPt.normMatrix(matMs, hh, true);
    }

    dethh = 1. / hh.determinant();
    hh.invert();

    // Loop on the derivative terms

    for (int igparam = 0; igparam < ngparam; igparam++)
    {
      for (int jref = 0; jref < ncorner; jref++)
      {
        int ipref = amesh->getApex(imesh, jref);
        int iad = getSGradAddress(ipref, igparam);

        // Update HH matrix

        _loadHHGrad(amesh, hhGrad, igparam, ipref);
        gradLogDetHH = _computeGradLogDetHH(amesh,igparam,ipref,hh,work,work2);

        if (amesh->getVariety() == 0)
          matPGradHPt.normMatrix(matP, hhGrad, true);
        else
          matPGradHPt.normMatrix(matMs, hhGrad, true);

        // Storing in the Map

        double ratio = sqrt(dethh * detMtM);
        double dratio = - 0.5 * gradLogDetHH * ratio;
        double S = 0.;
        for (int j0 = 0; j0 < ncorner-1; j0++)
        {
          int ip0 = amesh->getApex(imesh, j0);
          _mapTildeCUpdate(tab[iad],ip0,dratio / 6.);
          double s = 0.;
          for (int j1 = 0; j1 < ncorner-1; j1++)
          {
            int ip1 = amesh->getApex(imesh, j1);
            double vald  = matPGradHPt.getValue(j0, j1)  / 2.;
            double valdS = matPHHPt.getValue(j0, j1)  / 2.;
            vald = ratio * vald + dratio * valdS;
            s += vald;
            _mapGradUpdate(Mtab[iad], ip0, ip1, vald, tol);
          }
          int ip1 = amesh->getApex(imesh, ncorner - 1);
          _mapGradUpdate(Mtab[iad],ip0, ip1, -s, tol);
          _mapGradUpdate(Mtab[iad],ip1, ip0, -s, tol);
          S += s;
        }
        int ip0 = amesh->getApex(imesh, ncorner - 1);
        _mapTildeCUpdate(tab[iad],ip0,dratio / 6.);
        _mapGradUpdate(Mtab[iad], ip0, ip0, S, tol);
      }
    }
  }

  // Construct the SGrad member
  _resetGrad();
  _SGrad.resize(number);
  _TildeCGrad.resize(number);

  for (int i = 0; i < number; i++)
  {
    _SGrad[i] = _BuildSGradfromMap(Mtab[i]);
    if (_SGrad[i] == nullptr) return 1;
    _TildeCGrad[i] = _BuildTildeCGradfromMap(tab[i]);
    if (_TildeCGrad[i] == nullptr) return 1;
  }

  VectorDouble sqrtTildeC = VH::power(_TildeC, 0.5);
  VectorDouble invSqrtTildeC = VH::power(_TildeC, -0.5);
  VectorDouble tempVec = VH::inverse(_TildeC);
  VH::multiplyConstant(tempVec, -0.5);

  int ind = 0;
  MatrixSparse* tildeCGradMat = nullptr;
  MatrixSparse* A = nullptr;
  MatrixSparse* At = nullptr;
  for (int ipar = 0; ipar < _nCovAnisoGradParam; ipar++)
  {
    for (int iap = 0; iap < getSize(); iap++)
    {
      VectorDouble tildeCGrad = _TildeCGrad[ind]->getDiagonal();

      VH::multiplyInPlace(tildeCGrad, tempVec);
      _SGrad[ind]->prodNormDiagVecInPlace(invSqrtTildeC, 1);

      tildeCGradMat = MatrixSparse::diagVec(tildeCGrad);
      A = MatrixFactory::prodMatMat<MatrixSparse>(_S, tildeCGradMat);
      delete tildeCGradMat;
      At = A->transpose();
      A->addMatInPlace(*At);
      _SGrad[ind]->addMatInPlace(*A);
      delete At;
      delete A;
      ind++;
    }
  }

  return 0;
}

void ShiftOpMatrix::_mapTildeCUpdate(std::map<int, double> &tab,
                                 int ip0,
                                 double value,
                                 double tol)
{
  std::pair<std::map<int,double>::iterator, bool> ret;

  if (ABS(value) < tol) return;
  ret = tab.insert(std::pair<int, double>(ip0, value));
  if (!ret.second) ret.first->second += value;
}

VectorT<std::map<int,double>> ShiftOpMatrix::_mapTildeCCreate()const
{
  int number = _ndim * getSize();
  VectorT<std::map<int,double>> tab(number);
  for (int i = 0; i < number; i++)
  {
    tab.push_back(std::map<int,double>());
  }
  return tab;
}

VectorT<std::map<int, double>> ShiftOpMatrix::_mapCreate() const
{
  int size = getSize();
  if (size <= 0) my_throw("_mapCreate");
  VectorT<std::map<int, double>> tab(size);
  return tab;
}

VectorT<VectorT<std::map<int, double>>> ShiftOpMatrix::_mapVectorCreate() const
{
  int number = _nCovAnisoGradParam * getSize();
  if (number <= 0) my_throw("_mapVectorCreate");
  VectorT<VectorT<std::map<int, double>>> Mtab(number);
  for (int i = 0; i < number; i++)
  {
    Mtab[i] = _mapCreate();
  }
  return Mtab;
}

/**
 * Construct the _Lambda vector (Dimension: _napices)
 * @param amesh Description of the Mesh
 */
void ShiftOpMatrix::_buildLambda(const AMesh *amesh)
{
  int ndim = getNDim();
  int nvertex = amesh->getNApices();
  //int nmeshes = amesh->getNMeshes();
  auto cova = _getCovAniso();

  /* Load global matrices */

  _Lambda.clear();
  _Lambda.resize(nvertex, 0.);

  MatrixSquareSymmetric hh(ndim);
  //double param = cova->getParam();
  bool flagSphere = (amesh->getVariety() == 1);

  double correc = cova->getCorrec();
  //double sqdethh = 0.;

 if (flagSphere)
  {
    const ASpace *space = getDefaultSpaceSh().get();
    const SpaceSN *spaceSn = dynamic_cast<const SpaceSN*>(space);
    double r = 1.;
    if (spaceSn != nullptr) r = spaceSn->getRadius();
    double normalizing = cova->normalizeOnSphere(50); //useful only for Markov
    correc = pow(r, -2.) * normalizing;
    if (_isGlobalHH())
    {
      _loadHH(amesh, hh, 0);
      //factor = sqrt(hh.determinant());
    }
  }

  for (int ip = 0; ip < nvertex; ip++)
  {
    _Lambda[ip] = sqrt(_TildeC[ip] * correc);
  }
}

/**
 * Construct the _Lambda vector (Dimension: _napices)
 * @param amesh Description of the Mesh (New class)
 * @return
 */
bool ShiftOpMatrix::_buildLambdaGrad(const AMesh *amesh)
{
  int nvertex = amesh->getNApices();
  auto cova = cloneAndCast(_cova);

  /* Core allocation */

  int number = getLambdaGradSize();
  if (_LambdaGrad.empty())
  {
    _LambdaGrad.clear();
    VectorDouble temp(nvertex);
    for(int i = 0; i< number; i++)
      _LambdaGrad.push_back(temp);
  }

  // Filling by range / angle

  for (int ip = 0; ip < nvertex; ip++)
  {
    // Locally update the covariance for non-stationarity (if necessary)
    _updateCova(cova, ip);

    for (int idim = 0; idim < number; idim++)
    {
      _LambdaGrad[idim][ip] = -_Lambda[ip] / (2. * cova->getScale(idim));
    }
  }
  return false;
}

/**
 * Project the coordinates of the mesh vertices on the sphere
 * @param amesh Mesh structure
 * @param srot  Rotation parameters
 * @param imesh Rank of the mesh of interest
 * @param coeff Coordinates of the projected vertices
 */
void ShiftOpMatrix::_projectMesh(const AMesh *amesh,
                             const VectorDouble& srot,
                             int imesh,
                             double coeff[3][2])
{
  double xyz[3][3];

  // Calculate the Mesh Center

  VectorDouble center(3, 0.);
  for (int icorn = 0; icorn < (int) amesh->getNApexPerMesh(); icorn++)
  {
    GH::convertSph2Cart(amesh->getCoor(imesh, icorn, 0),
                        amesh->getCoor(imesh, icorn, 1), &xyz[icorn][0],
                        &xyz[icorn][1], &xyz[icorn][2]);
    for (int i = 0; i < 3; i++)
      center[i] += xyz[icorn][i];
  }
  double ratio = VH::norm(center);
  VH::divideConstant(center, sqrt(ratio));

  // Center gives the vector joining the origin to the center of triangle
  double phi    = srot[1] * GV_PI / 180.;
  double theta  = srot[0] * GV_PI / 180.;
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  double sintet = sin(theta);
  double costet = cos(theta);

  // W is the Pole vector
  VectorDouble w;
  w.push_back(sinphi * costet);
  w.push_back(sinphi * sintet);
  w.push_back(cosphi);

  // V1 = Center ^ w: first axis
  VectorDouble v1 = VH::crossProduct3D(center, w);
  VH::normalize(v1);

  // V2 = Center ^ V1: second axis
  VectorDouble v2 = VH::crossProduct3D(center, v1);
  VH::normalize(v2);

  // Get the end points from Unit vectors
  VectorDouble axe1 = VH::add(center, v1);
  VectorDouble axe2 = VH::add(center, v2);

  /* Projection */

  for (int icorn = 0; icorn < 3; icorn++)
  {
    coeff[icorn][0] = coeff[icorn][1] = 0.;
    for (int i = 0; i < 3; i++)
      coeff[icorn][0] += (axe1[i] - center[i]) * (xyz[icorn][i] - center[i]);
    for (int i = 0; i < 3; i++)
      coeff[icorn][1] += (axe2[i] - center[i]) * (xyz[icorn][i] - center[i]);
  }
}

/**
 * Returns the internal address for a given vertex and a given parameter
 * It returns -1 if the address is ivalid
 * @param iapex  Rank of the target apex
 * @param igparam  Rank of the target parameter
 * @return
 */
int ShiftOpMatrix::getSGradAddress(int iapex, int igparam) const
{
  int ngparam = _nCovAnisoGradParam;
  int napices = getSize();
  if (!checkArg("Mesh Apex index", iapex, napices)) return -1;
  if (!checkArg("Rank of the CovAniso parameter", igparam, ngparam)) return -1;
  return napices * igparam + iapex;
}

double ShiftOpMatrix::getMaxEigenValue() const
{
  return getS()->L1Norm();
}


int ShiftOpMatrix::getLambdaGradSize() const
{
  return _ndim;
}


void ShiftOpMatrix::_mapGradUpdate(std::map<std::pair<int, int>, double>& tab,
                               int ip0,
                               int ip1,
                               double value,
                               double tol)
{
  std::pair<std::map<std::pair<int, int>, double>::iterator, bool> ret;

  if (ABS(value) < tol) return;
  std::pair<int, int> key(ip0, ip1);
  ret = tab.insert(std::pair<std::pair<int, int>, double>(key, value)); // ret.second = false if key is already in the map
  if (!ret.second) ret.first->second += value;
}

/**
 * Transform the Map into a square cparse matrix
 * @param tab   Input Map
 * @return
 */
MatrixSparse* ShiftOpMatrix::_BuildSGradfromMap(std::map<std::pair<int, int>, double> &tab) const
{
  std::map<std::pair<int, int>, double>::iterator it;

  NF_Triplet NF_T;
  int ip0_max = -1;
  int ip1_max = -1;

  it = tab.begin();
  while (it != tab.end())
  {
    int ip0 = it->first.first;
    int ip1 = it->first.second;
    NF_T.add(ip0, ip1, it->second);
    if (ip0 > ip0_max) ip0_max = ip0;
    if (ip1 > ip1_max) ip1_max = ip1;
    it++;
  }

  // Add the fictitious value at maximum sparse matrix dimension (if 'nmax' provided)

  NF_T.force(getSize(), getSize());

  /* Optional printout */

  return MatrixSparse::createFromTriplet(NF_T, getSize(), getSize());
}

void ShiftOpMatrix::_determineFlagNoStatByHH()
{
  _flagNoStatByHH = false;
  if (! _isNoStat()) return;
  _flagNoStatByHH = _cova->isNoStatForTensor();
}

// void ShiftOpMatrix::multiplyByValueAndAddDiagonal(double v1,double v2)
// {
//   MatrixSparse* T1 = MatrixSparse::diagConstant(getSize(), 1.);
//   if (T1 == nullptr) my_throw("Problem in cs_eye");
//   _S->addMatInPlace(*T1, v1, v2);
//   delete T1;  
// }