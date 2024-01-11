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
#include "LinearOp/ShiftOpCs.hpp"
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Timer.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/ANoStat.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/Model.hpp"
#include "Space/SpaceSN.hpp"
#include "Space/ASpaceObject.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <math.h>

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

ShiftOpCs::ShiftOpCs(const CGParam params)
    : ALinearOp(params),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _TildeCGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _model(nullptr),
      _igrf(0),
      _icov(0),
      _ndim(0),
      _napices(0)
{
}

ShiftOpCs::ShiftOpCs(const AMesh* amesh,
                     Model* model,
                     const Db* dbout,
                     int igrf,
                     int icov,
                     const CGParam params,
                     bool verbose)
    : ALinearOp(params),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _TildeCGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _model(model),
      _igrf(0),
      _icov(0),
      _ndim(amesh->getEmbeddedNDim()),
      _napices(amesh->getNApices())
{
  (void) initFromMesh(amesh, model, dbout, igrf, icov, verbose);
}

#ifndef SWIG
ShiftOpCs::ShiftOpCs(const cs* S,
                     const VectorDouble& TildeC,
                     const VectorDouble& Lambda,
                     Model* model,
                     const CGParam params,
                     bool verbose)
    : ALinearOp(params),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _TildeCGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _model(model),
      _igrf(0),
      _icov(0),
      _ndim(0),
      _napices(cs_getncol(S))
{
  (void) initFromCS(S, TildeC, Lambda, model, verbose);
}
#endif

ShiftOpCs::ShiftOpCs(const ShiftOpCs &shift)
    : ALinearOp(shift),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _model(nullptr),
      _igrf(0),
      _icov(0),
      _ndim(0),
      _napices(0)
{
  _reallocate(shift);
}

ShiftOpCs& ShiftOpCs::operator=(const ShiftOpCs &shift)
{
  _reset();
  _reallocate(shift);
  return *this;
}

ShiftOpCs::~ShiftOpCs()
{
  _reset();
}

ShiftOpCs* ShiftOpCs::create(const AMesh *amesh,
                             Model *model,
                             const Db *dbout,
                             int igrf,
                             int icov,
                             const CGParam params,
                             bool verbose)
{
  return new ShiftOpCs(amesh, model, dbout, igrf, icov, params, verbose);
}

#ifndef SWIG
ShiftOpCs* ShiftOpCs::createFromSparse(const cs *S,
                                       const VectorDouble &TildeC,
                                       const VectorDouble &Lambda,
                                       Model *model,
                                       const CGParam params,
                                       bool verbose)
{
  return new ShiftOpCs(S, TildeC, Lambda, model, params, verbose);
}
#endif

/**
 *
 * @param amesh Meshing description (New format)
 * @param model Pointer to the Model structure
 * @param dbout Pointer to the Db structure
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance within the Model
 * @param flagAdvection When TRUE, S is replaced by G
 * @param verbose Verbose flag
 * @return Error return code
 */
int ShiftOpCs::initFromMesh(const AMesh* amesh,
                            Model* model,
                            const Db* /*dbout*/,
                            int igrf,
                            int icov,
                            bool flagAdvection,
                            bool verbose)
{
  DECLARE_UNUSED(flagAdvection);

  // Initializations

  _setModel(model);
  _setIgrf(igrf);
  _setIcov(icov);
  _napices = amesh->getNApices();
  try
  {
    _ndim = amesh->getEmbeddedNDim();

    // Attach the Non-stationary to Mesh and Db (optional)

    if (model->isNoStat())
    {
      if (model->getNoStat()->attachToMesh(amesh, verbose))
      {
        messerr("Problem when attaching 'mesh' to Non_stationary Parameters");
        return 1;
      }
    }

    // Calculating and storing the mesh sizes
    VectorDouble units = amesh->getMeshSizes();

    // Define if parameterization is in HH or in range/angle
    _determineFlagNoStatByHH();

    // Identify the covariance
    CovAniso cova = *model->getCova(icov);

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
 * @param model Pointer to the Model structure
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance within the Model
 * @param verbose Verbose flag
 * @param tol Smallest value below which the value is not stored in sparse matrix
 * @return Error return code
 */
int ShiftOpCs::initGradFromMesh(const AMesh* amesh,
                                Model* model,
                                int igrf,
                                int icov,
                                bool verbose,
                                double tol)
{
  // Initializations

  _setModel(model);
  _setIgrf(igrf);
  _setIcov(icov);

  try
  {

    // Attach the Non-stationary to Mesh and Db (optional)

    if (model->isNoStat())
    {
      if (model->getNoStat()->attachToMesh(amesh, verbose))
      {
        messerr("Problem when attaching 'mesh' to Non_stationary Parameters");
        return 1;
      }
    }

    // Identify the covariance
    CovAniso cova = *model->getCova(icov);

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
 * @param model Pointer to the Model structure
 * @param verbose Verbose flag
 * @return
 */
int ShiftOpCs::initFromCS(const cs* S,
                          const VectorDouble& TildeC,
                          const VectorDouble& Lambda,
                          Model* model,
                          bool verbose)
{
  DECLARE_UNUSED(verbose);

  // Initializations

  _setModel(model);

  try
  {
    // Store the TildeC & Lambda vectors

    _TildeC = TildeC;
    _Lambda = Lambda;
    _ndim = model->getDimensionNumber();

    // Duplicate the Shift Operator sparse matrix

    _S = cs_duplicate(S);
    if (_S == nullptr) my_throw("Problem when duplicating S sparse matrix");
    _napices = cs_getncol(S);
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
void ShiftOpCs::prodTildeC(const VectorDouble& x,
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

void ShiftOpCs::prodLambda(const VectorDouble& x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
  if (power == EPowerPT::ONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * _Lambda[i];
  }
  else if (power == EPowerPT::MINUSONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / _Lambda[i];
  }
  else if (power == EPowerPT::HALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * sqrt(_Lambda[i]);
  }
  else if (power == EPowerPT::MINUSHALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / sqrt(_Lambda[i]);
  }
  else
  {
    my_throw("Unexpected value for argument 'power'");
  }
}

void ShiftOpCs::prodLambdaOnSqrtTildeC(const VectorDouble& inv,
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
 ** \param[in] x       Input vector
 ** \param[in] y       Output vector
 **
 ** \remarks 'S' is a member that stands as a sparse matrix
 **
 *****************************************************************************/
void ShiftOpCs::_evalDirect(const VectorDouble& x, VectorDouble& y) const
{
  Timer time;
  int n = (int) x.size();
  cs_vecmult(_S, n, x.data(), y.data());
  getLogStats().incrementStatsDirect(time.getIntervalSeconds());
}

void ShiftOpCs::_resetGrad()
{
  if (_SGrad.empty()) return;
  for (int i = 0; i < (int) _SGrad.size(); i++)
     _SGrad[i] = cs_spfree(_SGrad[i]);
}

void ShiftOpCs::_reset()
{
  _S = cs_spfree(_S);
  _resetGrad();
}

void ShiftOpCs::_reallocate(const ShiftOpCs& shift)
{
  _TildeC = shift._TildeC;
  _Lambda = shift._Lambda;
  _S = cs_duplicate(shift._S);
  _nModelGradParam = shift._nModelGradParam;
  for (int i = 0; i < (int) _SGrad.size(); i++)
    _SGrad[i] = cs_duplicate(shift._SGrad[i]);
  for (int i = 0; i < (int) _LambdaGrad.size(); i++)
    _LambdaGrad[i] = shift._LambdaGrad[i];
  _flagNoStatByHH = shift._flagNoStatByHH;

  _model = shift._model;
  _igrf = shift._igrf;
  _icov = shift._icov;
  _ndim = shift._ndim;
  _napices = shift._napices;
}

cs* ShiftOpCs::getTildeCGrad(int iapex, int igparam) const
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

cs* ShiftOpCs::getSGrad(int iapex, int igparam) const
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

Triplet ShiftOpCs::getSToTriplet(bool flag_from_1) const
{
  return csToTriplet(getS(), flag_from_1);
}
Triplet ShiftOpCs::getTildeCGradToTriplet(int iapex, int igparam, bool flag_from_1) const
{
  return csToTriplet(getTildeCGrad(iapex, igparam), flag_from_1);
}
Triplet ShiftOpCs::getSGradToTriplet(int iapex, int igparam, bool flag_from_1) const
{
  return csToTriplet(getSGrad(iapex, igparam), flag_from_1);
}

/**
 * Locally update the covariance (does nothing in Stationary Case)
 * @param cova Local CovAniso structure (updated here)
 * @param imesh Rank of the active mesh
 */
void ShiftOpCs::_updateCova(CovAniso* cova, int imesh)
{
  if (! _isNoStat()) return;
  int igrf = _getIgrf();
  int icov = _getIcov();
  int ndim = getNDim();
  const ANoStat* nostat = _getModel()->getNoStat();

  // Third parameter
  if (nostat->isDefined(EConsElem::PARAM, icov, -1, -1, igrf))
  {
    double param = nostat->getValue(EConsElem::PARAM, 0, imesh, igrf, -1, -1, igrf);
    cova->setParam(param);
  }

  // Anisotropy coefficients
  if (nostat->isDefinedforRotation(icov, igrf))
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      if (nostat->isDefined(EConsElem::RANGE, icov, idim, -1, igrf))
      {
        double range = nostat->getValue(EConsElem::RANGE, 0, imesh, icov, idim, -1, igrf);
        cova->setRange(idim, range);
      }
      if (nostat->isDefined(EConsElem::SCALE, icov, idim, -1, igrf))
      {
        double scale = nostat->getValue(EConsElem::SCALE, 0, imesh, icov,  idim, -1, igrf);
        cova->setScale(idim, scale);
      }
    }

    // Anisotropy Rotation
    for (int idim = 0; idim < ndim; idim++)
    {
      if (nostat->isDefined(EConsElem::ANGLE, icov, idim, -1, igrf))
      {
        double anisoAngle = nostat->getValue(EConsElem::ANGLE, 0, imesh, icov, idim, -1, igrf);
        cova->setAnisoAngle(idim, anisoAngle);
      }
    }
  }
}

void ShiftOpCs::_updateHH(MatrixSquareSymmetric& hh, int imesh)
{
  if (! _isNoStat()) return;
  int igrf = _getIgrf();
  int icov = _getIcov();
  int ndim = getNDim();
  const ANoStat* nostat = _getModel()->getNoStat();

  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = idim; jdim < ndim; jdim++)
    {
      double value = nostat->getValue(EConsElem::TENSOR, 0, imesh, icov, idim, jdim, igrf);
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
void ShiftOpCs::_loadHHRegular(MatrixSquareSymmetric &hh, int imesh)
{
  int ndim = getNDim();
  const CovAniso *covini = _getCova();
  CovAniso *cova = covini->clone();

  // Locally update the covariance for non-stationarity (if necessary)
  _updateCova(cova, imesh);

  // Calculate the current HH matrix (using local covariance parameters)
  const MatrixSquareGeneral &rotmat = cova->getAnisoInvMat();

  VectorDouble diag = VH::power(cova->getScales(), 2.);
  MatrixSquareSymmetric temp(ndim);
  temp.setDiagonal(diag);
  hh.normMatrix(temp, rotmat);

  delete cova;
}

void ShiftOpCs::_loadHHVariety(MatrixSquareSymmetric& hh, int imesh)
{
  int ndim = getNDim();
  const CovAniso *covini = _getCova();
  CovAniso *cova = covini->clone();

  // Locally update the covariance for non-stationarity (if necessary)
  _updateCova(cova, imesh);

  // Calculate the current HH matrix (using local covariance parameters)
  VectorDouble diag = VH::power(cova->getScales(), 2.);

  hh.fill(0.);
  for (int idim = 0; idim < ndim; idim++)
    hh.setValue(idim, idim, diag[0]);

  delete cova;
}

/**
 * Calculate HH Gradient matrix from one of the Model parameters for the given Apex.
 * @param amesh Meshing description (New format)
 * @param hh         Output Array (updated here)
 * @param igparam    Rank of the parameter for derivation
 * @param ipref      Rank of the point
 *
 * @details: The parameters 'igparam' are sorted as follows:
 * @details: - 0:(ndim-1)   : ranges in each Space direction
 * @details: - ndim:ngparam : rotation angles (=ndim or 1 in 2-D)
 */

void ShiftOpCs::_loadHHGrad(const AMesh *amesh,
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

    const CovAniso *covini = _getCova();
    CovAniso *cova = covini->clone();

    // Locally update the covariance for non-stationarity (if necessary)

    _updateCova(cova, ipref);
    const MatrixSquareGeneral &rotmat = cova->getAnisoRotMat();
    VectorDouble diag = VH::power(cova->getScales(), 2.);

    MatrixSquareSymmetric temp(ndim);
    if (igparam < ndim)
    {
      // Derivation with respect to the Range 'igparam'
      temp.fill(0.);
      temp.setValue(igparam, igparam, 2. * cova->getScale(igparam));
      hh.normMatrix(temp, rotmat);
    }
    else
    {
      // Derivation with respect to the Angle 'igparam'-ndim
      int ir = igparam - ndim;
      CovAniso *covaderiv = covini->clone();
      _updateCova(covaderiv, ipref);
      cova->setAnisoAngle(ir, covaderiv->getAnisoAngles(ir) + 90.);
      const MatrixSquareGeneral &drotmat = covaderiv->getAnisoRotMat();

      VH::divideConstant(diag, 180. / GV_PI); // Necessary as angles are provided in degrees. Factor 2 is for derivative
      temp.setDiagonal(diag);
      hh.innerMatrix(temp, drotmat, rotmat);
    }
    delete cova;
  }

  int number = amesh->getNApexPerMesh();
  hh.prodScalar(1. / number);
}

double ShiftOpCs::_computeGradLogDetHH(const AMesh *amesh,
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
    const CovAniso *covaini = _getCova();
    CovAniso *cova = covaini->clone();
    _updateCova(cova, ipref);
    const MatrixSquareGeneral &rotmat = cova->getAnisoRotMat();
    MatrixSquareSymmetric temp(ndim);
    temp.setDiagonal(cova->getScales());
    for (int idim = 0; idim < ndim; idim++)
    {
      if (idim != igparam)
      {
        temp.setValue(idim, idim, 0.);
      }
      else
      {
        temp.setValue(idim, idim, 2. * cova->getScale(idim) / number);
      }
    }

    work.normMatrix(temp, rotmat);
    work2.prodMatrix(work, invHH);
    double result = work2.trace();
    delete cova;
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
void ShiftOpCs::_loadHH(const AMesh *amesh,
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

void ShiftOpCs::_loadAux(VectorDouble &tab, const EConsElem &type, int imesh)
{
  if (! _isNoStat()) return;
  if (tab.empty()) return;
  int igrf = _getIgrf();
  int icov = _getIcov();
  int size = static_cast<int> (tab.size());

  VectorDouble tabloc(size, 0.);
  for (int i = 0; i < size; i++) tab[i] = 0;

  const ANoStat* nostat = _getModel()->getNoStat();
  for (int i = 0; i < (int) tab.size(); i++)
    if (nostat->isDefined(type, icov, i, -1, igrf))
      tab[i] = nostat->getValue(type, 0, imesh, icov, i, -1, igrf);
}

int ShiftOpCs::_preparMatrices(const AMesh *amesh,
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

cs* ShiftOpCs::_BuildTildeCGradfromMap(std::map< int, double> &tab) const
{
  std::map<int, double>::iterator it;

  cs* Striplet = cs_spalloc(0, 0, 1, 1, 1);
  int ip1_max = -1;

  it = tab.begin();
  while (it != tab.end())
  {
    int ip1 = it->first;
    if (!cs_entry(Striplet, ip1, ip1, it->second)) return nullptr;
    if (ip1 > ip1_max) ip1_max = ip1;
    it++;
  }


    // Add the fictitious value at maximum sparse matrix dimension

  cs_force_dimension(Striplet, getSize(), getSize());

    /* Optional printout */

  cs* S = cs_triplet(Striplet);
  if (S == nullptr) return nullptr;
  Striplet = cs_spfree(Striplet);

  return S;
}

int ShiftOpCs::_prepareMatricesSVariety(const AMesh* amesh,
                                         int imesh,
                                         VectorVectorDouble& coords,
                                         MatrixSquareSymmetric& matMtM,
                                         AMatrix& matP,
                                         double *deter)
{
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();

  amesh->getEmbeddedCoordinatesPerMeshInPlace(imesh, coords);

  MatrixRectangular matM(ndim, ncorner - 1);
  for (int icorn = 0; icorn < ncorner - 1; icorn++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      double val = coords[icorn][idim] - coords[ncorner - 1][idim];
      matM.setValue(idim, icorn, val);
    }
  }

  // Calculate M^t %*% M

  matMtM.normSingleMatrix(matM);
  *deter = matMtM.determinant();

  // Calculate (M^t %*% M)^{-1}

  if (matMtM.invert())
  {
    messerr("Problem for Mesh #%d", imesh + 1);
    amesh->printMesh(imesh);
    return 1;
  }

  // Calculate P = (M^t %*% M)^{-1} %*% M^t
  matP.prodMatrix(matMtM, matM, false, true);
  return 0;
}

int ShiftOpCs::_prepareMatricesSphere(const AMesh *amesh,
                                      int imesh,
                                      VectorVectorDouble &coords,
                                      AMatrixSquare &matMs,
                                      double *deter)
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
int ShiftOpCs::_buildS(const AMesh *amesh, double tol)
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
  MatrixSquareGeneral matMs(ndim);
  MatrixSquareSymmetric matPinvHPt(ncorner-1);
  VectorVectorDouble coords = amesh->getEmbeddedCoordinatesPerMesh();
  double detMtM = 0.;
  double dethh = 0.;

  // Define the global matrices

  int igrf = _getIgrf();
  int icov = _getIcov();
  if (_isGlobalHH(igrf, icov))
  {
    _loadHH(amesh, hh, 0);
    dethh = 1. / hh.determinant();
  }
  if (! _isNoStat())
    _loadAux(srot, EConsElem::SPHEROT, 0);

  _S = _preparSSparse(amesh);
  if (_S == nullptr) goto label_end;

  /* Loop on the active meshes */

  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    OptDbg::setCurrentIndex(imesh + 1);

    // Non stationary case

    if (_isNoStat())
    {
      const ANoStat* nostat = _getModel()->getNoStat();
      if (nostat->isDefinedforAnisotropy(icov, igrf))
      {
        _loadHH(amesh, hh, imesh);
        dethh = 1. / hh.determinant();
      }
      if (nostat->isDefined(EConsElem::SPHEROT, icov, -1, -1, igrf))
        _loadAux(srot, EConsElem::SPHEROT, imesh);
    }

    // Prepare M matrix

    bool flagSphere = (amesh->getVariety() != 0);
    flagSphere = false; // Modif DR a valider
    if (! flagSphere)
    {
      if (_prepareMatricesSVariety(amesh, imesh, coords, matMtM, matP, &detMtM))
        my_throw("Matrix inversion");
      matPinvHPt.normTMatrix(hh, matP);
    }
    else
    {
      if (_prepareMatricesSphere(amesh, imesh, coords, matMs, &detMtM))
        my_throw("Matrix inversion");
      matPinvHPt.normTMatrix(hh, matMs);
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
        cs_add_value(_S, ip0, ip1, vald);
      }
      int ip1 = amesh->getApex(imesh, ncorner - 1);
      cs_add_value(_S, ip0, ip1, -s);
      cs_add_value(_S, ip1, ip0, -s);
      S += s;
    }
    int ip0 = amesh->getApex(imesh, ncorner - 1);
    _TildeC[ip0] += ratio / 6.;
    cs_add_value(_S, ip0, ip0, S);
  }

  // Ending S construction

  cs_matvecnorm_inplace(_S, _TildeC.data(), 2);

  /* Set the error return code */

  error = 0;

  label_end: if (error) _S = cs_spfree(_S);
  return error;
}

cs* ShiftOpCs::_preparSSparse(const AMesh *amesh) const
{
  cs *Sret = nullptr;
  cs* Striplet = nullptr;
  cs* Sl = nullptr;
  int nmeshes = amesh->getNMeshes();
  int ncorner = amesh->getNApexPerMesh();

  // Define Sl as the sparse matrix giving the clutter of apices among vertices
  Striplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int imesh = 0; imesh < nmeshes; imesh++)
  {
    for (int ic = 0; ic < ncorner; ic++)
    {
      int iapex = amesh->getApex(imesh, ic);
      if (! cs_entry(Striplet, iapex, imesh, 1.)) goto label_end;
    }
  }
  Sl = cs_triplet(Striplet);
  if (Sl == nullptr) goto label_end;

  // Operate the product Sl * t(Sl) to get the final matrix Sret
  Sret = cs_prod_norm_single(2, Sl);

  // Blank out the contents of the sparse matrix
  cs_set_cste(Sret, 0.);

label_end:
  Striplet = cs_spfree(Striplet);
  Sl = cs_spfree(Sl);
  return Sret;
}

bool ShiftOpCs::_cond(int indref, int igparam, int ipref)
{
  if (ipref == indref && igparam == 0) return true;
  return false;
}

/**
 * Calculate the private member "_SGrad" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSGrad(const AMesh *amesh, double tol)
{
  const CovAniso* cova = _getCova();
  _nModelGradParam = cova->getGradParamNumber();
  int number = _nModelGradParam * getSize();
  VectorT<std::map<int, double> > tab(number);
  std::vector<std::map<std::pair<int, int>, double> > Mtab(number);

  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int ngparam = _nModelGradParam;

  // Initialize the arrays

  MatrixSquareSymmetric hh(ndim);
  MatrixSquareSymmetric work(ndim);
  MatrixSquareSymmetric work2(ndim);
  MatrixSquareSymmetric hhGrad(ndim);
  MatrixSquareSymmetric matMtM(ncorner-1);
  MatrixRectangular matP(ncorner-1,ndim);
  MatrixSquareGeneral matMs(ndim);
  MatrixSquareSymmetric matPGradHPt(ncorner-1);
  MatrixSquareSymmetric matPHHPt(ncorner-1);

  double detMtM = 0.;
  double dethh = 0.;
  double gradLogDetHH = 0.;
  // Define the global matrices

  int igrf = _getIgrf();
  int icov = _getIcov();
  if (_isGlobalHH(igrf, icov))
    _loadHH(amesh, hhGrad, 0);

  /* Loop on the meshes */

  VectorVectorDouble coords = amesh->getEmbeddedCoordinatesPerMesh();
  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setCurrentIndex(imesh + 1);

    // Prepare M matrix
    _loadHH(amesh, hh, imesh);
    bool flagSphere = (amesh->getVariety() == 1);
    flagSphere = false; // Modif DR a valider
    if (! flagSphere)
    {
      if (_prepareMatricesSVariety(amesh, imesh, coords, matMtM, matP, &detMtM))
        my_throw("Matrix inversion");
      matPHHPt.normTMatrix(hh, matP);
    }

    else
    {
      if (_prepareMatricesSphere(amesh, imesh, coords, matMs, &detMtM))
        my_throw("Matrix inversion");
      matPHHPt.normTMatrix(hh, matMs);
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
          matPGradHPt.normTMatrix(hhGrad, matP);
        else
          matPGradHPt.normTMatrix(hhGrad, matMs);

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
  cs* tildeCGradMat = nullptr;
  cs* A = nullptr;
  cs* At = nullptr;
  for (int ipar = 0; ipar < _nModelGradParam; ipar++)
  {
    for (int iap = 0; iap < getSize(); iap++)
    {
      VectorDouble tildeCGrad = csd_extract_diag_VD(_TildeCGrad[ind], 1);

      VH::multiplyInPlace(tildeCGrad, tempVec);
      cs_matvecnorm_inplace(_SGrad[ind], invSqrtTildeC.data(), 0);

      tildeCGradMat = cs_diag(tildeCGrad);
      A = cs_multiply(_S, tildeCGradMat);
      cs_spfree(tildeCGradMat);
      At = cs_transpose(A, 1);
      A = cs_add_and_release(A, At, 1., 1., 1);
      cs_spfree(At);
      _SGrad[ind] = cs_add_and_release(_SGrad[ind], A, 1., 1., 1);
      cs_spfree(A);

      ind++;
    }
  }

  return 0;
}

void ShiftOpCs::_mapTildeCUpdate(std::map<int, double> &tab,
                                 int ip0,
                                 double value,
                                 double tol) const
{
  std::pair<std::map<int,double>::iterator, bool> ret;

  if (ABS(value) < tol) return;
  ret = tab.insert(std::pair<int, double>(ip0, value));
  if (!ret.second) ret.first->second += value;
}

VectorT<std::map<int,double>> ShiftOpCs::_mapTildeCCreate()const
{
  int number = _ndim * getSize();
  VectorT<std::map<int,double>> tab(number);
  for (int i = 0; i < number; i++)
  {
    tab.push_back(std::map<int,double>());
  }
  return tab;
}

VectorT<std::map<int, double>> ShiftOpCs::_mapCreate() const
{
  int size = getSize();
  if (size <= 0) my_throw("_mapCreate");
  VectorT<std::map<int, double>> tab(size);
  return tab;
}

VectorT<VectorT<std::map<int, double>>> ShiftOpCs::_mapVectorCreate() const
{
  int number = _nModelGradParam * getSize();
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
void ShiftOpCs::_buildLambda(const AMesh *amesh)
{
  int ndim = getNDim();
  int nvertex = amesh->getNApices();
  int nmeshes = amesh->getNMeshes();
  const CovAniso* cova = _getCova();

  /* Load global matrices */

  _Lambda.clear();
  _Lambda.resize(nvertex, 0.);

  int icov = _getIcov();
  int igrf = _getIgrf();
  MatrixSquareSymmetric hh(ndim);
  double param = cova->getParam();
  bool flagSphere = (amesh->getVariety() == 1);
  bool flagSill = false;
  if (_isNoStat())
  {
    const ANoStat *nostat = _getModel()->getNoStat();
    if (nostat->isDefined(EConsElem::SILL, icov)) flagSill = true;
  }

  double sill = cova->getSill(0, 0);
  double correc = cova->getCorrec();
  double sqdethh = 0.;
  double factor = 1.;

  if (flagSphere)
  {
    const ASpace *space = getDefaultSpace();
    const SpaceSN *spaceSn = dynamic_cast<const SpaceSN*>(space);
    double r = 1.;
    if (spaceSn != nullptr) r = spaceSn->getRadius();
    correc = cova->evalCovOnSphere(0, 50, false) * pow(r, 2. * param) / sill;

    if (_isGlobalHH(igrf, icov))
    {
      _loadHH(amesh, hh, 0);
      sqdethh = sqrt(hh.determinant());
      // The next term compensate all what is needed to run on the Sphere
        factor = pow(sqdethh, - (2. * param  - 1.)/3.);
    }
  }

  /* Fill the array */

  if (_isNoStat() && flagSill)
  {
    VectorDouble cum(nvertex, 0.);
    const ANoStat *nostat = _getModel()->getNoStat();
    for (int imesh = 0; imesh < nmeshes; imesh++)
    {
      if (flagSphere && nostat->isDefinedforAnisotropy(icov, igrf))
      {
        _loadHH(amesh, hh, imesh);
        sqdethh = sqrt(hh.determinant());
        factor = pow(sqdethh, - (2. * param  - 1.)/3.);
      }
 
      sill = nostat->getValue(EConsElem::SILL, 0, imesh, icov);

      for (int ic = 0, ncorner = amesh->getNApexPerMesh(); ic < ncorner; ic++)
      {
        int ip = amesh->getApex(imesh, ic);
        _Lambda[ip] += sill / factor;
        cum[ip]++;
      }
    }

    for (int ip = 0; ip < nvertex; ip++)
    {
      if (cum[ip] > 0.) _Lambda[ip] /= cum[ip];
      _Lambda[ip] = sqrt(_TildeC[ip] * correc / _Lambda[ip]);
    }
  }
  else
  {
    for (int ip = 0; ip < nvertex; ip++)
     {
      _Lambda[ip] = sqrt(_TildeC[ip] * correc * factor / sill);
     }
  }
}

/**
 * Construct the _Lambda vector (Dimension: _napices)
 * @param amesh Description of the Mesh (New class)
 * @return
 */
bool ShiftOpCs::_buildLambdaGrad(const AMesh *amesh)
{
  int nvertex = amesh->getNApices();
  const CovAniso* covini = _getCova();
  CovAniso* cova = covini->clone();

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

  delete cova;
  return false;
}

/**
 * Project the coordinates of the mesh vertices on the sphere
 * @param amesh Mesh structure
 * @param srot  Rotation parameters
 * @param imesh Rank of the mesh of interest
 * @param coeff Coordinates of the projected vertices
 */
void ShiftOpCs::_projectMesh(const AMesh *amesh,
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
int ShiftOpCs::getSGradAddress(int iapex, int igparam) const
{
  int ngparam = _nModelGradParam;
  int napices = getSize();
  if (iapex < 0 || iapex >= napices)
  {
    mesArg("Mesh Apex index", iapex, napices);
    return -1;
  }
  if (igparam < 0 || igparam >= ngparam)
  {
    mesArg("Rank of the Model parameter", igparam, ngparam);
    return -1;
  }
  return napices * igparam + iapex;
}

double ShiftOpCs::getMaxEigenValue() const
{
  return cs_norm(getS());
}

bool ShiftOpCs::_isNoStat()
{
  const Model* model = _getModel();
  return model->isNoStat();
}

bool ShiftOpCs::_isGlobalHH(int igrf, int icov)
{
  if (! _isNoStat())
    return true;
  else
  {
    const ANoStat* nostat = _getModel()->getNoStat();
    if (! nostat->isDefinedforAnisotropy(icov, igrf)) return true;
  }
  return false;
}

int ShiftOpCs::getLambdaGradSize() const
{
  return _ndim;
}

const CovAniso* ShiftOpCs::_getCova()
{
  return _getModel()->getCova(_getIcov());
}

void ShiftOpCs::_mapGradUpdate(std::map<std::pair<int, int>, double> &tab,
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
cs* ShiftOpCs::_BuildSGradfromMap(std::map<std::pair<int, int>, double> &tab)
{
  std::map<std::pair<int, int>, double>::iterator it;

  cs* Striplet = cs_spalloc(0, 0, 1, 1, 1);
  int ip0_max = -1;
  int ip1_max = -1;

  it = tab.begin();
  while (it != tab.end())
  {
    int ip0 = it->first.first;
    int ip1 = it->first.second;
    if (!cs_entry(Striplet, ip0, ip1, it->second)) return nullptr;
    if (ip0 > ip0_max) ip0_max = ip0;
    if (ip1 > ip1_max) ip1_max = ip1;
    it++;
  }

  // Add the fictitious value at maximum sparse matrix dimension (if 'nmax' provided)

  cs_force_dimension(Striplet, getSize(), getSize());

  /* Optional printout */

  cs* S = cs_triplet(Striplet);
  if (S == nullptr) return nullptr;

  Striplet = cs_spfree(Striplet);

  return S;
}

void ShiftOpCs::_determineFlagNoStatByHH()
{
  _flagNoStatByHH = false;
  if (! _isNoStat()) return;
  _flagNoStatByHH = _getModel()->getNoStat()->isDefinedByType(EConsElem::TENSOR, _getIgrf());
}
