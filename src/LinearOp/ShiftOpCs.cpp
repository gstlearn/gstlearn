/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "csparse_f.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovAniso.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Model/ANoStat.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/Model.hpp"
#include "Space/SpaceSN.hpp"
#include "Space/ASpaceObject.hpp"

#include <math.h>

ShiftOpCs::ShiftOpCs()
    : ALinearOp(),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _variety(0),
      _model(nullptr),
      _igrf(0),
      _icov(0),
      _ndim(0)
{
}

ShiftOpCs::ShiftOpCs(const AMesh* amesh,
                     Model* model,
                     const Db* dbout,
                     int igrf,
                     int icov,
                     bool verbose)
    : ALinearOp(),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _variety(amesh->getVariety()),
      _model(model),
      _igrf(0),
      _icov(0),
      _ndim(amesh->getEmbeddedNDim())
{
  (void) initFromMesh(amesh, model, dbout, igrf, icov, verbose);
}

ShiftOpCs::ShiftOpCs(const cs* S,
                     const VectorDouble& TildeC,
                     const VectorDouble& Lambda,
                     Model* model,
                     bool verbose)
    : ALinearOp(),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _LambdaGrad(),
      _flagNoStatByHH(false),
      _model(model),
      _igrf(0),
      _icov(0),
      _ndim(0)
{
  _variety = 0;
  (void) initFromCS(S, TildeC, Lambda, model, verbose);
}

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
      _ndim(shift.getNDim())
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

/**
 *
 * @param s_mesh Mesh description (Old format)
 * @param model Pointer to the Model structure
 * @param dbout Pointer to Db which
 * @param flagAdvection When TRUE, S is replaced by G
 * @param verbose Verbose flag
 * @return Error return code
 */
int ShiftOpCs::initFromOldMesh(SPDE_Mesh* s_mesh,
                               Model* model,
                               Db* dbout,
                               bool flagAdvection,
                               bool verbose)
{
  double* units;

  int error = 0;
  units = (double *) nullptr;
  _setModel(model);

  try
  {
    // Not stationarity is not coded in this deprecated code

    if (model->isNoStat())
      my_throw("Non-stationarity is not coded in this deprecated Shiftop code");

    // Attach the Model

    if (spde_check(NULL, dbout, model, NULL, verbose, VectorDouble(), 0, 0, 0,
                   1, 1, 0, 0))
    my_throw("Problem with spde_check() method");

    // Converting Meshing

    MeshEStandard amesh;
    amesh.convertFromOldMesh(s_mesh, 0);

    // Calculate the meshes of the vertices

    units = _spde_get_mesh_dimension(&amesh);
    if (units == (double *) nullptr)
    my_throw("Problem with spde_get_mesh_dimension() method");

    // Construct G sparse Matrix (locally stored in _S)

    _S = _spde_fill_S(&amesh, model, units);
    if (_S == nullptr) my_throw("Problem with spde_fill_S() method");

    // Construct the TildeC vector

    _TildeC = _spde_fill_TildeC(&amesh, units);
    if (_TildeC.empty())
    my_throw("Problem with spde_fill_TildeC() method");

    // Construct the Lambda vector

    _Lambda = _spde_fill_Lambda(model, &amesh, _TildeC);
    if (_Lambda.empty())
    my_throw("Problem with _spde_fill_Lambda() method");

    // Construct the final Sparse matrix S

    if (!flagAdvection) cs_matvecnorm_inplace(_S, _TildeC.data(), 2);
  }

  catch (const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }

  units = (double *) mem_free((char * ) units);
  if (error) _reset();
  return error;
}

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
  // Initializations

  int error = 0;
  _setModel(model);
  _setIgrf(igrf);
  _setIcov(icov);
  _variety = amesh->getVariety();
  try
  {
    if (verbose) message(">>> Using the new calculation module <<<\n");
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

    // Attach the Model
//    if (spde_check(NULL, dbout, model, NULL, verbose, VectorDouble(), 0, 0, 0,
//                   1, 1, 0, 0))
//    my_throw("Problem with spde_check() method");

    // Identify the covariance
    CovAniso cova = *model->getCova(icov);

    // Construct S sparse Matrix
    if (amesh->getVariety() == 0)
    {
      if (!_isVelocity())
      {
        if (_buildSVariety(amesh))
        my_throw("Problem when buildS");
      }
      else
      {
        if (_buildSVel(amesh))
        my_throw("Problem when buildSVel");
      }
    }
    else
    {
      if (_buildSVariety(amesh)) // TODO or _buildSSphere
      my_throw("Problem when buildSVariety");
    }

    // Construct the TildeC vector

    if (amesh->getVariety() == 2)
    {
      if (_buildTildeC(amesh, units))
        my_throw("Problem with buildTildeC");
    }

    // Construct the Lambda vector

    _buildLambda(amesh);

    // Construct S sparse Matrix (locally stored in _S)

    if (! flagAdvection)
      cs_matvecnorm_inplace(_S, _TildeC.data(), 2);
  }

  catch (const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }

  if (error) _reset();
  return error;
}

/**
 * Initialize the environment for calculation of derivatives of S
 * @param amesh Meshing description (New format)
 * @param model Pointer to the Model structure
 * @param dbout Pointer to the Db structure
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance within the Model
 * @param verbose Verbose flag
 * @param tol Smallest value below which the value is not stored in sparse matrix
 * @return Error return code
 */
int ShiftOpCs::initGradFromMesh(const AMesh* amesh,
                                Model* model,
                                Db* dbout,
                                int igrf,
                                int icov,
                                bool verbose,
                                double tol)
{
  // Initializations

  int error = 0;
  _setModel(model);
  _setIgrf(igrf);
  _setIcov(icov);

  try
  {
    if (verbose) message(">>> Using the new calculation module <<<\n");

    // Attach the Non-stationary to Mesh and Db (optional)

    if (model->isNoStat())
    {
      if (model->getNoStat()->attachToMesh(amesh, verbose))
      {
        messerr("Problem when attaching 'mesh' to Non_stationary Parameters");
        return 1;
      }
    }

    if (spde_check(NULL, dbout, model, NULL, verbose, VectorDouble(), 0, 0, 0,
                   1, 1, 0, 0))
    my_throw("Problem with spde_check() method");

    // Identify the covariance
    CovAniso cova = *model->getCova(icov);

    // Construct S sparse Matrix
    if (_buildSGrad(amesh, tol))
      my_throw("Problem when buildSGrad");

    if (_buildLambdaGrad(amesh))
      my_throw("Problem when buildLambdaGrad");
  }

  catch (const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }

  if (error) _reset();
  return error;
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
  // Initializations

  int error = 0;
  _setModel(model);

  try
  {
    // Attach the Model

    if (spde_check(NULL, NULL, model, NULL, verbose, VectorDouble(), 0, 0, 0, 1,
                   1, 0, 0))
    my_throw("Problem with spde_check() method");

    // Store the TildeC & Lambda vectors

    _TildeC = TildeC;
    _Lambda = Lambda;
    _ndim = model->getDimensionNumber();

    // Duplicate the Shift Operator sparse matrix

    _S = cs_duplicate(S);
    if (_S == nullptr) my_throw("Problem when duplicating S sparse matrix");
  }

  catch (const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }

  return error;
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

void ShiftOpCs::prodLambdaOnSqrtTildeC(const VectorDouble& in,
                                     VectorDouble& out,
                                     double puis) const
{
  for (int i = 0; i < getSize(); i++)
    out[i] = in[i] * pow(_Lambda[i] / sqrt(_TildeC[i]), puis);
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
  int n = (int) x.size();
  cs_vecmult(_S, n, x.data(), y.data());
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

/**
 * Locally update the covariance (does nothing in Stationary Case)
 * @param cova Local CovAniso structure (updated here)
 * @param ip   Rank of the point
 */
void ShiftOpCs::_updateCova(CovAniso* cova, int ip)
{
  // Initializations
  if (! _isNoStat()) return;
  int igrf = _getIgrf();
  int icov = _getIcov();
  int ndim = getNDim();
  const ANoStat* nostat = _getModel()->getNoStat();

  // Third parameter
  if (nostat->isDefined(igrf, icov, EConsElem::PARAM, -1, -1))
  {
    double param = nostat->getValue(igrf, icov, EConsElem::PARAM, -1, -1, 0, ip);
    cova->setParam(param);
  }

  // Anisotropy coefficients
  if (nostat->isDefinedforRotation(igrf, icov))
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      if (nostat->isDefined(igrf, icov, EConsElem::RANGE, idim, -1))
      {
        double range = nostat->getValue(igrf, icov, EConsElem::RANGE, idim, -1, 0, ip);
        cova->setRange(idim, range);
      }
      if (nostat->isDefined(igrf, icov, EConsElem::SCALE, idim, -1))
      {
        double scale = nostat->getValue(igrf, icov, EConsElem::SCALE, idim, -1, 0, ip);
        cova->setScale(idim, scale);
      }
    }

    // Anisotropy Rotation
    for (int idim = 0; idim < ndim; idim++)
    {
      if (nostat->isDefined(igrf, icov, EConsElem::ANGLE, idim, -1))
      {
        double anisoAngle = nostat->getValue(igrf, icov, EConsElem::ANGLE, idim,-1, 0, ip);
        cova->setAnisoAngle(idim, anisoAngle);
      }
    }
  }
}

void ShiftOpCs::_updateHH(MatrixSquareSymmetric& hh, int ip)
{
  // Initializations
  if (! _isNoStat()) return;
  if (! _flagNoStatByHH) return;
  int igrf = _getIgrf();
  int icov = _getIcov();
  int ndim = getNDim();
  const ANoStat* nostat = _getModel()->getNoStat();

  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = idim; jdim < ndim; jdim++)
    {
      if (nostat->isDefined(igrf, icov, EConsElem::TENSOR, idim, jdim))
      {
       double value = nostat->getValue(igrf, icov, EConsElem::TENSOR, idim, jdim, 0, ip);
       hh.setValue(idim, jdim, value);
     }
   }
}

/**
 * Calculate HH matrix from parameters
 * Note that this function is also called in Stationary case...
 * So it must be workable without any updates
 * @param amesh AMesh structure
 * @param hh Output Array
 * @param ip Rank of the point
 */
void ShiftOpCs::_loadHHByApex(const AMesh *amesh,
                              MatrixSquareSymmetric& hh,
                              int ip)
{
  if (amesh->getVariety() == 0)
    _loadHHRegularByApex(hh, ip);
  else
    _loadHHVarietyByApex(hh, ip);
}

/**
 * Calculate HH matrix from parameters
 * Note that this function is also called in Stationary case...
 * So it must be workable without any updates
 * @param hh Output Array
 * @param ip Rank of the point
 */
void ShiftOpCs::_loadHHRegularByApex(MatrixSquareSymmetric& hh,
                                     int ip)
{
  int ndim = getNDim();
  const CovAniso* covini = _getCova();
  CovAniso* cova = covini->clone();

  if (_flagNoStatByHH)
  {
    _updateHH(hh, ip);
  }
  else
  {
    // Locally update the covariance for non-stationarity (if necessary)
    _updateCova(cova, ip);

    // Calculate the current HH matrix (using local covariance parameters)
    const MatrixSquareGeneral& rotmat = cova->getAnisoInvMat();
    VectorDouble diag = ut_vector_power(cova->getScales(), 2.);
    MatrixSquareSymmetric temp(ndim);
    temp.setDiagonal(diag);
    hh.normMatrix(temp, rotmat);
  }
  delete cova;
}

void ShiftOpCs::_loadHHVarietyByApex(MatrixSquareSymmetric& hh, int /*ip*/)
{
  int ndim = getNDim();
  const CovAniso* covini = _getCova();
  CovAniso* cova = covini->clone();

  if (_flagNoStatByHH)
  {
    messerr("To be implemented");
  }
  else
  {
//     Locally update the covariance for non-stationarity (if necessary)
   // _updateCova(cova, ip);

    // Calculate the current HH matrix (using local covariance parameters)
    VectorDouble diag = ut_vector_power(cova->getScales(), 2.);

    hh.fill(0.);
    for (int idim = 0; idim < ndim; idim++)
      hh.setValue(idim,idim, diag[0]);
  }
  delete cova;
}

/**
 * Calculate HH Gradient matrix from one of the Model parameters
 * for the given Apex.
 * @param hh      Output Array (updated here)
 * @param igparam Rank of the parameter for derivation
 * @param ip      Rank of the point
 *
 * @details: The parameters 'igparam' are sorted as follows:
 * @details: - 0:(ndim-1)   : ranges in each Space direction
 * @details: - ndim:ngparam : rotation angles (=ndim or 1 in 2-D)
 */
void ShiftOpCs::_loadHHGradByApex(MatrixSquareSymmetric& hh,
                                  int igparam,
                                  int ip)
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

    const CovAniso* covini = _getCova();
    CovAniso* cova = covini->clone();

    const MatrixSquareGeneral& rotmat = cova->getAnisoRotMat();
    VectorDouble diag = ut_vector_power(cova->getScales(), 2.);

    // Locally update the covariance for non-stationarity (if necessary)
    _updateCova(cova, ip);

    MatrixSquareSymmetric temp(ndim);
    if (igparam < ndim)
    {
      // Derivation with respect to the Range 'igparam'
      temp.fill(0);
      temp.setValue(igparam, igparam, 2. * cova->getScale(igparam));
      hh.normMatrix(temp, rotmat);
    }
    else
    {
      // Derivation with respect to the Angle 'igparam'-ndim
      int ir = igparam - ndim;
      CovAniso* dcova = cova;
      dcova->setAnisoAngle(ir, dcova->getAnisoAngles(ir) + 90.);
      const MatrixSquareGeneral& drotmat = dcova->getAnisoRotMat();
      ut_vector_divide_inplace(diag, 180. / GV_PI); // Necessary as angles are provided in degrees
      temp.setDiagonal(diag);
      hh.innerMatrix(temp, rotmat, drotmat);
    }
    delete cova;
  }
}

void ShiftOpCs::_loadAux(VectorDouble& tab,
                         const EConsElem& type,
                         int ip)
{
  if (tab.empty()) return;
  for (int i = 0; i < (int) tab.size(); i++) tab[i] = 0.;
  if (! _isNoStat()) return;
  int igrf = _getIgrf();
  int icov = _getIcov();

  const ANoStat* nostat = _getModel()->getNoStat();
  for (int i = 0; i < (int) tab.size(); i++)
    if (nostat->isDefined(igrf, icov, type, i, -1))
    {
      tab[i] = nostat->getValue(igrf, icov, type, i, -1, 0, ip);
    }
}

/**
 * Constitue HH (only in non-stationary case)
 * @param hh   Returned HH symmetric matrix
 * @param amesh Pointer to the meshing
 * @param imesh Rank of the mesh
 */
void ShiftOpCs::_loadHHPerMesh(const AMesh* amesh,
                               MatrixSquareSymmetric& hh,
                               int imesh)
{
  int number = amesh->getNApexPerMesh();
  int ndim = _ndim;
  MatrixSquareSymmetric hhloc(ndim);
  hh.fill(0.);

  // HH per mesh is obtained as the average of the HH per apex of the mesh
  for (int rank = 0; rank < number; rank++)
  {
    int ip = amesh->getApex(imesh, rank);
    _loadHHByApex(amesh, hhloc, ip);
    hh.add(hhloc);
  }
  hh.prodScalar(1. / number);
}

/**
 * Calculate the derivative of HH matrix with respect to
 * - the Model parameter 'igparam' and the Apex 'igp0'
 * @param hh      Resulting HH derivative matrix
 * @param amesh   Meshing structure
 * @param igp0    Rank of the Apex for derivation (between 0 to nApices-1)
 * @param igparam Rank of the Model parameter (from 0 to ngparam-1)
 */
void ShiftOpCs::_loadHHGradPerMesh(MatrixSquareSymmetric& hh,
                                   const AMesh* amesh,
                                   int igp0,
                                   int igparam)
{
  int number = amesh->getNApexPerMesh();
  hh.fill(0.);
  _loadHHGradByApex(hh, igparam, igp0);
  hh.prodScalar(1. / number);
}

void ShiftOpCs::_loadAuxPerMesh(const AMesh* amesh,
                                VectorDouble& tab,
                                const EConsElem& type,
                                int imesh)
{
  if (tab.empty()) return;
  int number = amesh->getNApexPerMesh();
  int size = static_cast<int> (tab.size());

  VectorDouble tabloc(size, 0.);
  for (int i = 0; i < size; i++) tab[i] = 0;

  for (int rank = 0; rank < number; rank++)
  {
    int ip = amesh->getApex(imesh, rank);
    _loadAux(tabloc, type, ip);
    ut_vector_add_inplace(tab, tabloc);
  }
  ut_vector_divide_inplace(tab, (double) number);
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
    amesh->printMeshes(imesh);
    return 1;
  }

  for (int icorn = 0; icorn < ncorner; icorn++)
    for (int idim = 0; idim < ndim; idim++)
      matw.setValue(idim, icorn, matu.getValue(icorn, idim));

  return 0;
}

/**
 * Transform the Map into a square cparse matrix
 * @param tab   Input Map
 * @param nmax  Dimension of the matrix (if provided)
 * @return
 *
 * @remark When 'nmax' is provided, make sure that the resulting sparse matrix
 * @remark has the full dimension (by adding a fictitious value equal to eps)
 */
cs* ShiftOpCs::_BuildSfromMap(std::map<std::pair<int, int>, double> &tab, int nmax)
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

  if (nmax > 0)
    cs_force_dimension(Striplet, nmax, nmax);

  /* Optional printout */

  cs* S = cs_triplet(Striplet);
  if (S == nullptr) return nullptr;

  Striplet = cs_spfree(Striplet);

  return S;
}

/**
 * Calculate the private member "_SGrad" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSGrad(const AMesh *amesh,
                           double tol)
{
  // Store the number of derivation parameters for the model as member
  const CovAniso* cova = _getCova();
  _nModelGradParam = cova->getGradParamNumber();
  int number = _nModelGradParam * getSize();
  std::vector<std::map<std::pair<int, int>, double> > Mtab(number);

  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int ngparam = _nModelGradParam;

  // Initialize the arrays

  VectorDouble matv(ncorner);
  MatrixSquareSymmetric hh(ndim);
  MatrixSquareGeneral matu(ncorner);
  MatrixSquareGeneral mat(ncorner);
  MatrixRectangular matw(ndim, ncorner);

  // Define the global matrices

  int igrf = _getIgrf();
  int icov = _getIcov();
  if (_isGlobalHH(igrf, icov))
    _loadHHByApex(amesh, hh, 0);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setIndex(imesh + 1);
    double meshSize = amesh->getMeshSize(imesh);

    // Calculate geometry

    if (_preparMatrices(amesh, imesh, matu, matw))
    my_throw("Problem in matrix inversion");

    // Loop on the derivative terms

    for (int j2 = 0; j2 < ncorner; j2++)
    {
      int igp0 = amesh->getApex(imesh, j2);
      for (int igparam = 0; igparam < ngparam; igparam++)
      {
        int iad = getSGradAddress(igp0, igparam);

        // Loop on apices of the current mesh

        for (int j0 = 0; j0 < ncorner; j0++)
          for (int j1 = 0; j1 < ncorner; j1++)
          {
            int ip0 = amesh->getApex(imesh, j0);
            int ip1 = amesh->getApex(imesh, j1);
            _loadHHGradPerMesh(hh, amesh, igp0, igparam);
            mat.normMatrix(hh, matw);

            double vald = mat.getValue(j0, j1) * meshSize;
            _mapUpdate(Mtab[iad], ip0, ip1, vald, tol);
          }
        }
      }
  }

  // Construct the SGrad member
  _resetGrad();
  _SGrad.resize(number);
  for (int i = 0; i < (int) Mtab.size(); i++)
  {
    _SGrad[i] = cs_spfree(_SGrad[i]);
    _SGrad[i] = _BuildSfromMap(Mtab[i]);
    if (_SGrad[i] == nullptr) goto label_end;
    cs_matvecnorm_inplace(_SGrad[i], _TildeC.data(), 2);
  }

  /* Set the error return code */

  error = 0;

  label_end:
  if (error)
    for (int i = 0; i < (int) _SGrad.size(); i++)
    _SGrad[i] = cs_spfree(_SGrad[i]);
  return error;
}

void ShiftOpCs::_mapUpdate(std::map<std::pair<int, int>, double>& tab,
                           int ip0,
                           int ip1,
                           double value,
                           double tol)
{
  std::pair<std::map<std::pair<int, int>, double>::iterator, bool> ret;

  if (ABS(value) < tol) return;
  std::pair<int, int> key(ip0, ip1);
  ret = tab.insert(std::pair<std::pair<int, int>, double>(key, value));
  if (!ret.second) ret.first->second += value;
}

/**
 * Calculate the private member "_S" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildS(const AMesh *amesh, double tol)
{
  std::map<std::pair<int, int>, double> tab;

  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();

  // Initialize the arrays

  MatrixSquareSymmetric hh(ndim);
  MatrixSquareGeneral matu(ncorner);
  MatrixSquareGeneral mat(ncorner);
  MatrixRectangular matw(ndim, ncorner);

  // Define the global HH matrix

  int igrf = _getIgrf();
  int icov = _getIcov();
  if (_isGlobalHH(igrf, icov))
    _loadHHByApex(amesh, hh, 0);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setIndex(imesh + 1);
    double meshSize = amesh->getMeshSize(imesh);

    // Case of Euclidean geometry

    if (_preparMatrices(amesh, imesh, matu, matw))
    my_throw("Problem in matrix inversion");

    // Non stationary case

    if (_isNoStat())
    {
      const ANoStat* nostat = _getModel()->getNoStat();
      if (nostat->isDefinedforAnisotropy(igrf, icov))
        _loadHHPerMesh(amesh, hh, imesh);
    }
    mat.normMatrix(hh, matw);

    for (int j0 = 0; j0 < ncorner; j0++)
      for (int j1 = 0; j1 < ncorner; j1++)
      {
        int ip0 = amesh->getApex(imesh, j0);
        int ip1 = amesh->getApex(imesh, j1);
        double vald = mat.getValue(j0, j1) * meshSize;
        _mapUpdate(tab, ip0, ip1, vald, tol);
      }
  }

  _S = cs_spfree(_S);
  _S = _BuildSfromMap(tab);
  if (_S == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end:
  if (error) _S = cs_spfree(_S);
  return error;
}

/**
 * Calculate the private member "_S" directly from the Mesh on a Variety in 3-D space
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 *
 * @remark TildeC is calculated at the same time
 */
int ShiftOpCs::_buildSVariety(const AMesh *amesh, double tol)
{
  std::map<std::pair<int, int>, double> tab;

  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int napices = amesh->getNApices();

  _TildeC.clear();
  _TildeC.resize(napices, 0.);

  // Initialize the arrays

  VectorDouble srot(2);
  MatrixSquareSymmetric hh(ndim);
  MatrixSquareSymmetric matMtM(ncorner-1);
  MatrixRectangular matP(ncorner-1,ndim);
  MatrixSquareGeneral matMs(ndim);
  MatrixSquareSymmetric matPinvHPt(ncorner-1);
  double detMtM = 0.;
  // Define the global matrices

  int igrf = _getIgrf();
  int icov = _getIcov();
  double dethh = 0.;
  if (_isGlobalHH(igrf, icov))
  {
    _loadHHByApex(amesh, hh, 0);
    dethh = 1./hh.determinant();
  }
  if (! _isNoStat())
    _loadAux(srot, EConsElem::SPHEROT, 0);

  /* Loop on the meshes */

  VectorVectorDouble coords = amesh->getEmbeddedCoordinatesPerMesh();
  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setIndex(imesh + 1);

    // Non stationary case

    if (_isNoStat())
    {
      const ANoStat* nostat = _getModel()->getNoStat();
      if (nostat->isDefinedforAnisotropy(igrf, icov))
      {
        _loadHHPerMesh(amesh, hh, imesh);
        dethh = 1./hh.determinant();

      }
      if (nostat->isDefined(igrf, icov, EConsElem::SPHEROT, -1, -1))
        _loadAuxPerMesh(amesh, srot, EConsElem::SPHEROT, imesh);
    }

    // Load M matrix [ndim, ncorner-1]
    amesh->getEmbeddedCoordinatesPerMesh(imesh, coords);

    MatrixRectangular matM(ndim, ncorner-1);

    for (int icorn = 0; icorn < ncorner-1; icorn++)
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        double val = coords[icorn][idim] - coords[ncorner-1][idim];
        if(amesh->getVariety()==1)
        {
          matM.setValue(idim, icorn, val);
        }
        else
        {
          matMs.setValue(idim, icorn, val);
        }
      }
    }

    if(amesh->getVariety() == 1)
    {
      // Calculate M^t %*% M
      matMtM.normSingleMatrix(matM);
      detMtM = matMtM.determinant();

      // Calculate (M^t %*% M)^{-1}

      if (matMtM.invert())
      {
        messerr("Problem for Mesh #%d", imesh + 1);
        amesh->printMeshes(imesh);
        my_throw("Matrix inversion");
      }
      // Calculate P = (M^t %*% M)^{-1} %*% M^t
      matM.transposeInPlace();
      matP.prodMatrix(matMtM, matM);

      // Calculate P %*% HH^{-1} %*% P^t

      matPinvHPt.normTMatrix(hh, matP);
    }
    else
    {
      double detM = matMs.determinant();
      detMtM = detM * detM;
      if(matMs.invert())
      {
        messerr("Problem for Mesh #%d", imesh + 1);
        amesh->printMeshes(imesh);
        my_throw("Matrix inversion");

      }
      matPinvHPt.normTMatrix(hh, matMs);
    }
    // Storing in the Map

    double ratio = sqrt(dethh * detMtM);
    double S = 0.;
    for (int j0 = 0; j0 < ncorner-1; j0++)
    {
      // Update TildeC

      int ip0 = amesh->getApex(imesh, j0);
      _TildeC[ip0] += ratio / 6.;

      double s = 0.;
      for (int j1 = 0; j1 < ncorner-1; j1++)
      {
        int ip1 = amesh->getApex(imesh, j1);
        double vald = matPinvHPt.getValue(j0, j1) * ratio / 2.;
        s += vald;
        _mapUpdate(tab, ip0, ip1, vald, tol);
      }
      int j1 = ncorner - 1;
      int ip1 = amesh->getApex(imesh, j1);
      _mapUpdate(tab, ip0, ip1, -s, tol);
      _mapUpdate(tab, ip1, ip0, -s, tol);
      S += s;
    }
    int j0 = ncorner-1;

    int ip0 = amesh->getApex(imesh, j0);
    _TildeC[ip0] += ratio / 6.;
    _mapUpdate(tab, ip0, ip0, S, tol);
  }

  _S = cs_spfree(_S);
  _S = _BuildSfromMap(tab, napices);
  if (_S == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: if (error) _S = cs_spfree(_S);
  return error;
}

/**
 * Calculate the private member "_S" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSSphere(const AMesh *amesh,
                             double tol)
{
  std::map<std::pair<int, int>, double> tab;
  double coeff[3][2];

  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();

  // Initialize the arrays

  VectorDouble srot(2), axe1(3), axe2(3), vel(3), matv(ncorner);
  MatrixSquareSymmetric hh(ndim);
  MatrixSquareGeneral matu(ncorner);
  MatrixSquareGeneral mat(ncorner);
  MatrixRectangular matw(ndim, ncorner);

  // Define the global matrices

  int igrf = _getIgrf();
  int icov = _getIcov();
  if (_isGlobalHH(igrf, icov))
    _loadHHByApex(amesh, hh, 0);
  if (! _isNoStat())
  {
    _loadAux(srot, EConsElem::SPHEROT, 0);
    _loadAux(vel, EConsElem::VELOCITY, 0);
  }

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setIndex(imesh + 1);

    // Non stationary case

    if (_isNoStat())
    {
      const ANoStat* nostat = _getModel()->getNoStat();
      if (nostat->isDefinedforAnisotropy(igrf, icov))
        _loadHHPerMesh(amesh, hh, imesh);
      if (nostat->isDefined(igrf, icov, EConsElem::SPHEROT, -1, -1))
        _loadAuxPerMesh(amesh, srot, EConsElem::SPHEROT, imesh);
      if (nostat->isDefined(igrf, icov, EConsElem::VELOCITY, -1, -1))
        _loadAuxPerMesh(amesh, vel, EConsElem::VELOCITY, imesh);
    }

    // Case of Spherical geometry

    _projectMesh(amesh, srot, imesh, coeff);

    for (int icorn = 0; icorn < ncorner; icorn++)
    {
      for (int idim = 0; idim < ndim; idim++)
        matu.setValue(idim, icorn, coeff[icorn][idim]);
      matu.setValue(ncorner - 1, icorn, 1.);
    }

    if (matu.invert())
    {
      messerr("Problem for Mesh #%d", imesh + 1);
      amesh->printMeshes(imesh);
      my_throw("Matrix inversion");
    }

    for (int icorn = 0; icorn < ncorner; icorn++)
      for (int idim = 0; idim < ndim; idim++)
        matw.setValue(idim, icorn, matu.getValue(icorn, idim));

    // Update for Advection (non-stationary)

    if (_isVelocity())
    {
      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        matv[icorn] = 0.;
        for (int idim = 0; idim < ndim; idim++)
          matv[icorn] += vel[idim] * matw.getValue(idim, icorn);
      }
    }
    else
    {
      mat.normMatrix(hh, matw);
    }

    for (int j0 = 0; j0 < ncorner; j0++)
      for (int j1 = 0; j1 < ncorner; j1++)
      {
        int ip0 = amesh->getApex(imesh, j0);
        int ip1 = amesh->getApex(imesh, j1);
        double vald = (_isVelocity()) ? matv[j1] : mat.getValue(j0, j1);
        _mapUpdate(tab, ip0, ip1, vald * amesh->getMeshSize(imesh), tol);
      }
  }

  _S = cs_spfree(_S);
  _S = _BuildSfromMap(tab);
  if (_S == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: if (error) _S = cs_spfree(_S);
  return error;
}

/**
 * Calculate the private member "_S" directly from the Mesh for Velocity
 * @param amesh Description of the Mesh (New class)
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSVel(const AMesh *amesh,
                          double tol)
{
  std::map<std::pair<int, int>, double> tab;

  int error = 1;
  int ndim = getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int igrf = _getIgrf();
  int icov = _getIcov();

  // Initialize the arrays

  VectorDouble vel(2);
  VectorDouble matv(ncorner);
  MatrixSquareGeneral matu(ncorner);
  MatrixRectangular matw(ndim, ncorner);

  // Define the global HH matrix

  _loadAux(vel, EConsElem::VELOCITY, 0);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    OptDbg::setIndex(imesh + 1);
    double meshSize = amesh->getMeshSize(imesh);

    // Non stationary case

    const ANoStat* nostat = _getModel()->getNoStat();
    if (nostat->isDefined(igrf, icov, EConsElem::VELOCITY, -1, -1))
      _loadAuxPerMesh(amesh, vel, EConsElem::VELOCITY, imesh);

    // Case of Euclidean geometry

    if (_preparMatrices(amesh, imesh, matu, matw))
      my_throw("Problem in matrix inversion");

    // Update for Advection (non-stationary)

    for (int icorn = 0; icorn < ncorner; icorn++)
    {
      matv[icorn] = 0.;
      for (int idim = 0; idim < ndim; idim++)
        matv[icorn] += vel[idim] * matw.getValue(idim, icorn);
    }

    for (int j0 = 0; j0 < ncorner; j0++)
      for (int j1 = 0; j1 < ncorner; j1++)
      {
        int ip0 = amesh->getApex(imesh, j0);
        int ip1 = amesh->getApex(imesh, j1);
        double vald = matv[j1] * meshSize;
        _mapUpdate(tab, ip0, ip1, vald, tol);
      }
  }

  _S = cs_spfree(_S);
  _S = _BuildSfromMap(tab);
  if (_S == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end:
  if (error) _S = cs_spfree(_S);
  return error;
}

/**
 * Calculate _TildeC directly from the Mesh (Dimension: _NApices)
 * @param amesh Description of the Mesh (New class)
 * @param units Array of sizes for all meshes
 * @return Error return code
 */
int ShiftOpCs::_buildTildeC(const AMesh *amesh, const VectorDouble& units)
{
  int nvertex = amesh->getNApices();
  int ncorner = amesh->getNApexPerMesh();
  double factor = (double) ncorner;

  /* Core allocation */

  VectorDouble cumunit(nvertex, 0.);
  _TildeC.clear();

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {

    /* Loop on the vertices */

    for (int icorn = 0; icorn < ncorner; icorn++)
    {
      int jp = amesh->getApex(imesh, icorn);
      cumunit[jp] += units[imesh];
    }
  }

  /* Scale */

  for (int ip = 0; ip < nvertex; ip++)
  {
    double value = cumunit[ip] / factor;
    if (ABS(value) <= 0.)
    {
      messerr("Meshing unit (%d) has a zero volume", ip + 1);
      _TildeC.clear();
      return 1;
    }
    _TildeC.push_back(value);
  }
  return 0;
}

/**
 * Construct the _Lambda vector (Dimension: _NApices)
 * @param amesh Description of the Mesh (New class)
 * @return
 */
void ShiftOpCs::_buildLambda(const AMesh *amesh)
{
  double sqdeth = 0.;

  int ndim = getNDim();
  int nvertex = amesh->getNApices();
  const CovAniso* cova = _getCova();
  double param = cova->getParam();
  double r = 1.;
  if( amesh->getVariety() == 1)
  {
    const ASpace* space = ASpaceObject::getDefaultSpace();
    const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
    r = spaceSn->getRadius();
  }

  /* Load global matrices */

  _Lambda.clear();
  int igrf = _getIgrf();
  int icov = _getIcov();
  MatrixSquareSymmetric hh(ndim);
  double correc  = cova->getCorrec();
  if (_isGlobalHH(igrf, icov))
  {
    _loadHHByApex(amesh, hh, 0);
    sqdeth = sqrt(hh.determinant());
   if(amesh->getVariety() == 1)
   {
     correc = cova->evalCovOnSphere(0,50,false);
   }
  }

  /* Fill the array */

  double sill = cova->getSill(0, 0);
  for (int ip = 0; ip < nvertex; ip++)
  {
    if (_isNoStat())
    {
      const ANoStat* nostat = _getModel()->getNoStat();
      if (nostat->isDefinedforAnisotropy(igrf, icov))
      {
        _loadHHByApex(amesh, hh, ip);
        sqdeth = sqrt(hh.determinant());
      }
      if (nostat->isDefined(igrf, icov, EConsElem::SILL, -1, -1))
        sill = nostat->getValue(igrf, icov, EConsElem::SILL, -1, -1, 0, ip);
    }
    if (amesh->getVariety() != 1)
    {
      _Lambda.push_back(sqrt((_TildeC[ip]) * correc / sill));
    }
    else
    {
        _Lambda.push_back(sqrt((_TildeC[ip]) *  correc * pow(r, 2. * param)  * pow(sqdeth, - (2. * param  - 1.)/3.)/ (sill*sill)));
    }
  }
}

/**
 * Construct the _Lambda vector (Dimension: _NApices)
 * @param amesh Description of the Mesh (New class)
 * @return
 */
bool ShiftOpCs::_buildLambdaGrad(const AMesh *amesh)
{
  int ndim = getNDim();
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

  /* Fill the array */

  if (_flagNoStatByHH)
  {

    // Filling by HH

    MatrixSquareSymmetric hh(ndim);
    for (int ip = 0; ip < nvertex; ip++)
    {
      // Update HH locally
      _updateHH(hh, ip);
      hh.invert();

      int ecr = 0;
      for(int idim = 0; idim < ndim; idim++)
      {
        for(int jdim = 0; jdim <= idim; jdim++, ecr++)
        {
          double ratio = (idim == jdim)? 1./4. : 1./2.;
          _LambdaGrad[ecr][ip] = - _Lambda[ip] * hh.getValue(idim,jdim) * ratio ;
        }
      }
    }
  }
  else
  {

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
    util_convert_sph2cart(amesh->getCoor(imesh, icorn, 0),
                          amesh->getCoor(imesh, icorn, 1), &xyz[icorn][0],
                          &xyz[icorn][1], &xyz[icorn][2]);
    for (int i = 0; i < 3; i++)
      center[i] += xyz[icorn][i];
  }
  double ratio = ut_vector_norm(center);
  ut_vector_divide_inplace(center, sqrt(ratio));

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
  VectorDouble v1 = ut_vector_cross_product(center, w);
  ut_vector_norm(v1);

  // V2 = Center ^ V1: second axis
  VectorDouble v2 = ut_vector_cross_product(center, v1);
  ut_vector_norm(v2);

  // Get the end points from Unit vectors
  VectorDouble axe1 = ut_vector_add(center, v1);
  VectorDouble axe2 = ut_vector_add(center, v2);

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
  int napex = getSize();
  if (iapex < 0 || iapex >= napex)
  {
    mesArg("Mesh Apex index", iapex, napex);
    return -1;
  }
  if (igparam < 0 || igparam >= ngparam)
  {
    mesArg("Rank of the Model parameter", igparam, ngparam);
    return -1;
  }
  return napex * igparam + iapex;
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
    if (! nostat->isDefinedforAnisotropy(igrf, icov)) return true;
  }
  return false;
}

bool ShiftOpCs::_isVelocity()
{
  if (! _isNoStat()) return false;
  const ANoStat* nostat = _getModel()->getNoStat();
  int igrf = _getIgrf();
  int icov = _getIcov();
  return nostat->isDefined(igrf, icov, EConsElem::VELOCITY, -1, -1);
}

int ShiftOpCs::getLambdaGradSize() const
{
  if (_flagNoStatByHH)
    return _nModelGradParam;
  else
    return _ndim;
}

void ShiftOpCs::_determineFlagNoStatByHH()
{
  _flagNoStatByHH = false;
  if (! _isNoStat()) return;
  _flagNoStatByHH = _getModel()->getNoStat()->isDefinedByType(_getIgrf(), EConsElem::TENSOR);
}

const CovAniso* ShiftOpCs::_getCova()
{
  return _getModel()->getCova(_getIcov());
}
