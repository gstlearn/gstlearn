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
#include "MatrixC/MatrixCSGeneral.hpp"
#include "MatrixC/MatrixCRectangular.hpp"
#include "MatrixC/MatrixCSSym.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_e.h"

#include "Basic/AException.hpp"
#include "Covariances/CovAniso.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Model/ANoStat.hpp"

ShiftOpCs::ShiftOpCs()
    : ALinearOp(),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _LambdaGrad(),
      _dim(0)
{
}

ShiftOpCs::ShiftOpCs(AMesh* amesh,
                     Model* model,
                     const Db* dbout,
                     ANoStat* nostat,
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
      _dim(amesh->getNDim())
{
  (void) initFromMesh(amesh, model, dbout, nostat, igrf, icov, verbose);
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
      _dim(0)
{
  (void) initFromCS(S, TildeC, Lambda, model, verbose);
}

ShiftOpCs::ShiftOpCs(const ShiftOpCs &shift)
    : ALinearOp(shift),
      _TildeC(),
      _Lambda(),
      _S(nullptr),
      _nModelGradParam(0),
      _SGrad(),
      _dim(shift.getDim())
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
 * @param dbout Pointer to Db which contains the NoStat information
 * @param nostat Non-stationary structure (optional)
 * @param flagAdvection When TRUE, S is replaced by G
 * @param verbose Verbose flag
 * @return Error return code
 */
int ShiftOpCs::initFromOldMesh(SPDE_Mesh* s_mesh,
                             Model* model,
                             Db * dbout,
                             ANoStat* nostat,
                             bool flagAdvection,
                             bool verbose)
{
  double* units;

  // Initializations

  int error = 0;
  units = (double *) nullptr;

  try
  {
    // Not stationarity is not coded in this deprecated code
    if (nostat == nullptr)
    my_throw("Non-stationarity is not coded in this deprecated Shiftop code");

    // Attach the Model

    if (spde_check(NULL, dbout, model, NULL, verbose, VectorDouble(), 0, 0, 0,
                   1, 1, 0, 0))
    my_throw("Problem with spde_check() method");

    // Converting Meshing

    MeshEStandard amesh;
    amesh.convertFromOldMesh(s_mesh, 0);

    // Calculate the meshes of the vertices

    units = spde_get_mesh_dimension(&amesh);
    if (units == (double *) nullptr)
    my_throw("Problem with spde_get_mesh_dimension() method");

    // Construct G sparse Matrix (locally stored in _S)

    _S = spde_fill_S(&amesh, model, units);
    if (_S == (cs *) NULL) my_throw("Problem with spde_fill_S() method");

    // Construct the TildeC vector

    _TildeC = spde_fill_TildeC(&amesh, units);
    if (_TildeC.empty())
    my_throw("Problem with spde_fill_TildeC() method");

    // Construct the Lambda vector

    _Lambda = spde_fill_Lambda(model, &amesh, _TildeC);
    if (_Lambda.empty())
    my_throw("Problem with spde_fill_Lambda() method");

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
 * @param dbout Pointer to the Db structure (used for nostat)
 * @param nostat Pointer to the ANoStat for non-stationary parameters
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance within the Model
 * @param flagAdvection When TRUE, S is replaced by G
 * @param verbose Verbose flag
 * @return Error return code
 */
int ShiftOpCs::initFromMesh(AMesh* amesh,
                            Model* model,
                            const Db* dbout,
                            ANoStat* nostat,
                            int igrf,
                            int icov,
                            bool flagAdvection,
                            bool verbose)
{
  // Initializations

  int error = 0;

  try
  {
    if (verbose) message(">>> Using the new calculation module <<<\n");

    // Attach the Non-stationary to Mesh and Db (optional)

    if (nostat != nullptr)
    {
      if (nostat->attachModel(model))
      {
        messerr("Model and Non Stationary parameters are incompatible");
        return 1;
      }
      NoStatArray* nostatarray = dynamic_cast<NoStatArray*>(nostat);
      if (nostatarray->attachDb(dbout, 2, verbose))
      {
        messerr("Problem when attaching 'dbout' to Non-Stationary Parameters");
        return 1;
      }
      if (nostatarray->attachMesh(dbout, amesh, verbose))
      {
        messerr("Problem when attaching 'mesh' to Non_stationary Parameters");
        return 1;
      }
    }

    // Calculating and storing the mesh sizes
    VectorDouble units = amesh->getMeshSizes();

    // Attach the Model
    if (spde_check(NULL, dbout, model, NULL, verbose, VectorDouble(), 0, 0, 0,
                   1, 1, 0, 0))
    my_throw("Problem with spde_check() method");
    int flag_vel = 0;
    if (nostat != NULL)
      flag_vel = nostat->isDefined(igrf, icov, CONS_VELOCITY, -1, -1);

    // Identify the covariance
    CovAniso cova = *model->getCova(icov);

    // Construct S sparse Matrix
    if (amesh->getVariety() == 0)
    {
      if (!flag_vel)
      {
        if (_buildS(amesh, cova, igrf, icov, nostat))
        my_throw("Problem when buildS");
      }
      else
      {
        if (_buildSVel(amesh, cova, igrf, icov, nostat))
        my_throw("Problem when buildSVel");
      }
    }
    else
    {
      if (_buildSSphere(amesh, cova, igrf, icov, nostat))
      my_throw("Problem when buildSSphere");
    }

    // Construct the TildeC vector

    if (_buildTildeC(amesh, units))
    my_throw("Problem with buildTildeC");

    // Construct the Lambda vector

    _buildLambda(amesh, cova, igrf, icov, nostat);

    // Construct S sparse Matrix (locally stored in _S)

    if (! flagAdvection)
      cs_matvecnorm_inplace(_S, _TildeC.data(), 2);

    _dim = amesh->getNDim();
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
 * @param dbout Pointer to the Db structure (used for nostat)
 * @param nostat Pointer to the ANoStat for non-stationary parameters
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance within the Model
 * @param verbose Verbose flag
 * @param tol Smallest value below which the value is not stored in sparse matrix
 * @return Error return code
 */
int ShiftOpCs::initGradFromMesh(AMesh* amesh,
                                Model* model,
                                Db* dbout,
                                ANoStat* nostat,
                                int igrf,
                                int icov,
                                bool verbose,
                                double tol)
{
  // Initializations

  int error = 0;

  try
  {
    if (verbose) message(">>> Using the new calculation module <<<\n");

    // Attach the Non-stationary to Mesh and Db (optional)

    if (nostat != nullptr)
    {
      if (nostat->attachModel(model))
      {
        messerr("Model and Non Stationary parameters are incompatible");
        return 1;
      }
      NoStatArray* nostatarray = dynamic_cast<NoStatArray*>(nostat);
      if (nostatarray->attachMesh(dbout, amesh, verbose))
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
    if (_buildSGrad(amesh, cova, igrf, icov, nostat, tol))
      my_throw("Problem when buildSGrad");

    if (_buildLambdaGrad(amesh, cova, igrf, icov, nostat))
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

  try
  {
    // Attach the Model

    if (spde_check(NULL, NULL, model, NULL, verbose, VectorDouble(), 0, 0, 0, 1,
                   1, 0, 0))
    my_throw("Problem with spde_check() method");

    // Store the TildeC & Lambda vectors

    _TildeC = TildeC;
    _Lambda = Lambda;
    _dim = model->getDimensionNumber();
    // Duplicate the Shift Operator sparse matrix

    _S = cs_duplicate(S);
    if (_S == (cs *) NULL) my_throw("Problem when duplicating S sparse matrix");
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
                           ENUM_POPTS power) const
{
  if (power == POPT_ONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * _TildeC[i];
  }
  else if (power == POPT_MINUSONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / _TildeC[i];
  }
  else if (power == POPT_HALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * sqrt(_TildeC[i]);
  }
  else if (power == POPT_MINUSHALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / sqrt(_TildeC[i]);
  }
  else if (power == POPT_LOG)
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
                         ENUM_POPTS power) const
{
  if (power == POPT_ONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * _Lambda[i];
  }
  else if (power == POPT_MINUSONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] / _Lambda[i];
  }
  else if (power == POPT_HALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] = x[i] * sqrt(_Lambda[i]);
  }
  else if (power == POPT_MINUSHALF)
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
  cs_vecmult(_S, x.data(), y.data());
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
  for (int i = 0; i < (int) _SGrad.size(); i++)
    _SGrad[i] = cs_duplicate(shift._SGrad[i]);
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
 * Locally update the covariance
 * @param cova Local CovAniso structure
 * @param igrf Rank of the underlying GRF
 * @param icov Rank for the target Covariance (within Model)
 * @param ip   Rank of the point
 * @param ndim Space dimension
 * @param nostat Non-stationary parameters
 */
void ShiftOpCs::_updateCova(CovAniso& cova,
                            int igrf,
                            int icov,
                            int ip,
                            int ndim,
                            ANoStat *nostat)
{
 if (nostat->isDefined(igrf, icov, CONS_PARAM, -1, -1))
 {
   double param = nostat->getValue(igrf, icov, CONS_PARAM, -1, -1, 0, ip);
   cova.setParam(param);
 }

 // Anisotropy coefficients
   for (int idim = 0; idim < ndim; idim++)
   {
     if (nostat->isDefined(igrf, icov, CONS_RANGE, idim, -1))
     {
       double range = nostat->getValue(igrf, icov, CONS_RANGE, idim, -1, 0, ip);
       cova.setRange(idim, range);
     }
     if (nostat->isDefined(igrf, icov, CONS_SCALE, idim, -1))
     {
       double scale = nostat->getValue(igrf, icov, CONS_SCALE, idim, -1, 0, ip);
       cova.setScale(idim, scale);
     }
   }

 // Anisotropy Rotation
   for (int idim = 0; idim < ndim; idim++)
   {
     if (nostat->isDefined(igrf, icov, CONS_ANGLE, idim, -1))
     {
       double anisoAngle = nostat->getValue(igrf, icov, CONS_ANGLE, idim, -1, 0, ip);
       cova.setAnisoAngle(idim, anisoAngle);
     }
 }
}

/**
 * Calculate HH matrix from parameters
 * @param hh Output Array
 * @param covini CovAniso structure
 * @param igrf Rank of the Grf
 * @param icov Rank of the covariance
 * @param ip   Rank of the point
 * @param nostat Non-stationary parameters
 */
void ShiftOpCs::_loadHHByApex(MatrixCSGeneral& hh,
                              const CovAniso& covini,
                              int igrf,
                              int icov,
                              int ip,
                              ANoStat* nostat)
{
  CovAniso cova = covini;
  bool flagNostat = (nostat != nullptr);
  int ndim = hh.getNSize();

  // Locally update the covariance for non-stationarity

  if (flagNostat)
    _updateCova(cova, igrf, icov, ip, ndim, nostat);

  const MatrixCSGeneral& rotmat = cova.getAnisoRotMat();
  VectorDouble diag = ut_vector_power(cova.getScales(), 2.);
  MatrixCSSym temp(ndim);
  temp.setDiagonal(diag);
  hh.normMatrix(temp, rotmat);
}

/**
 * Calculate HH Gradient matrix from one of the Model parameters
 * for the given Apex.
 * @param hh Output Array
 * @param covini Covariance structure
 * @param igparam Rank of the parameter for derivation
 * @param igrf Rank of the Grf
 * @param icov Rank of the covariance
 * @param ip   Rank of the point
 * @param nostat Non-stationary parameters
 * @param flagFormal True for Formal calculations; False for Numerical Derivation
 *
 * @details: The parameters 'igparam' are sorted as follows:
 * @details: - 0:(ndim-1)   : ranges in each Space direction
 * @details: - ndim:ngparam : rotation angles (=ndim or 1 in 2-D)
 */
void ShiftOpCs::_loadHHGradByApex(MatrixCSGeneral& hh,
                                  const CovAniso& covini,
                                  int igparam,
                                  int igrf,
                                  int icov,
                                  int ip,
                                  ANoStat* nostat,
                                  bool flagFormal)
{
  CovAniso cova = covini;
  bool flagNostat = (nostat != nullptr);
  int ndim = hh.getNSize();

  // Locally update the covariance for non-stationarity

  if (flagNostat)
    _updateCova(cova, igrf, icov, ip, ndim, nostat);
  const MatrixCSGeneral& rotmat = cova.getAnisoRotMat();
  VectorDouble diag = ut_vector_power(cova.getScales(), 2.);

  // Main dispatch

  MatrixCSSym temp(ndim);

  if (!flagFormal)
  {

    // Numerical integration

    double eps = 1.e-05;
    MatrixCSGeneral hhp = MatrixCSGeneral(ndim);
    MatrixCSGeneral hhm = MatrixCSGeneral(ndim);

    if (igparam >= ndim)
    {
     // Angle
      int ir = igparam - ndim;

      CovAniso covap = cova;
      covap.setAnisoAngle(ir,covap.getAnisoAngles(ir) + eps);
      const MatrixCSGeneral& rotmatp = covap.getAnisoRotMat();
      temp.setDiagonal(diag);
      hhp.normMatrix(temp, rotmatp);

      CovAniso covam = cova;
      covam.setAnisoAngle(ir,covam.getAnisoAngles(ir) - eps);
      const MatrixCSGeneral& rotmatm = covam.getAnisoRotMat();
      temp.setDiagonal(diag);
      hhm.normMatrix(temp, rotmatm);
    }
    else
    {
      // Scale
      CovAniso covap = cova;
      covap.setScale(igparam, covap.getScale(igparam) + eps);
      VectorDouble diagp = ut_vector_power(covap.getScales(), 2.);
      temp.setDiagonal(diagp);
      hhp.normMatrix(temp, rotmat);

      CovAniso covam = cova;
      covam.setScale(igparam, covam.getScale(igparam) - eps);
      VectorDouble diagm = ut_vector_power(covam.getScales(), 2.);
      temp.setDiagonal(diagm);
      hhm.normMatrix(temp, rotmat);
    }

    // Calculate the finite difference

    for (int i = 0; i < ndim; i++)
      for (int j = 0; j < ndim; j++)
        hh.setValue(i, j, (hhp.getValue(i, j) - hhm.getValue(i, j)) / (2. * eps));
  }
  else
  {

    // Formal derivation

    if (igparam < ndim)
    {
      // Derivation with respect to the Range 'igparam'
      temp.fill(0);
      temp.setValue(igparam,igparam, 2. * cova.getScale(igparam));
      hh.normMatrix(temp, rotmat);
    }
    else
    {
      // Derivation with respect to the Angle 'igparam'-ndim
      int ir = igparam - ndim;
      CovAniso dcova = cova;
      dcova.setAnisoAngle(ir, dcova.getAnisoAngles(ir) + 90.);
      const MatrixCSGeneral& drotmat = dcova.getAnisoRotMat();
      ut_vector_divide_inplace(diag, 180. / GV_PI); // Necessary as angles are provided in degrees
      temp.setDiagonal(diag);
      hh.innerMatrix(temp, rotmat, drotmat);
    }
  }
}

void ShiftOpCs::_loadAux(VectorDouble& tab,
                         int igrf,
                         int icov,
                         ENUM_CONS type,
                         int ip,
                         ANoStat *nostat)
{
  if (tab.empty()) return;
  for (int i = 0; i < (int) tab.size(); i++)
    tab[i] = 0.;

  if (nostat == nullptr) return;

  for (int i = 0; i < (int) tab.size(); i++)
    if (nostat->isDefined(igrf, icov, type, i, -1))
    {
      tab[i] = nostat->getValue(igrf, icov, type, i, -1, 0, ip);
    }
}

void ShiftOpCs::_loadHHPerMesh(MatrixCSGeneral& hh,
                               AMesh* amesh,
                               const CovAniso& cova,
                               int igrf,
                               int icov,
                               int imesh,
                               ANoStat* nostat)
{
  int number = amesh->getNApexPerMesh();
  int ndim = amesh->getNDim();
  MatrixCSGeneral hhloc(ndim);
  hh.fill(0.);

  // HH per mesh is obtained as the average of the HH per apex of the mesh
  for (int rank = 0; rank < number; rank++)
  {
    int ip = amesh->getApex(imesh, rank);
    _loadHHByApex(hhloc, cova, igrf, icov, ip, nostat);
    hh.add(hhloc);
  }
  hh.prodScalar(1. / number);

}

/**
 * Calculate the derivative of HH matrix with respect to
 * - the Model parameter 'igparam' and the Apex 'igp0'
 * @param hh      Resulting HH derivative matrix
 * @param amesh   Meshing structure
 * @param cova    CovAniso structure
 * @param igp0    Rank of the Apex for derivation (between 0 to nApices-1)
 * @param igparam Rank of the Model parameter (from 0 to ngparam-1)
 * @param igrf    Rank of the current GRF
 * @param icov    Rank of the current Covariance (within Model)
 * @param imesh   Rank of the current mesh
 * @param nostat  Nostat structure
 * @param flagFormal True for Formal calculations; False for Numerical Derivation
 */
void ShiftOpCs::_loadHHGradPerMesh(MatrixCSGeneral& hh,
                                   AMesh* amesh,
                                   const CovAniso& cova,
                                   int igp0,
                                   int igparam,
                                   int igrf,
                                   int icov,
                                   int imesh,
                                   ANoStat* nostat,
                                   bool flagFormal)
{
  int number = amesh->getNApexPerMesh();
  hh.fill(0.);
  _loadHHGradByApex(hh, cova, igparam, igrf, icov, igp0, nostat, flagFormal);
  hh.prodScalar(1. / number);
}

void ShiftOpCs::_loadAuxPerMesh(VectorDouble& tab,
                              AMesh* amesh,
                              int igrf,
                              int icov,
                              ENUM_CONS type,
                              int imesh,
                              ANoStat* nostat)
{
  if (tab.empty()) return;
  int number = amesh->getNApexPerMesh();
  int size = tab.size();

  VectorDouble tabloc(size, 0.);
  for (int i = 0; i < size; i++)
    tab[i] = 0;

  for (int rank = 0; rank < number; rank++)
  {
    int ip = amesh->getApex(imesh, rank);
    _loadAux(tabloc, igrf, icov, type, ip, nostat);
    ut_vector_add_inplace(tab, tabloc);
  }
  ut_vector_divide_inplace(tab, (double) number);
}

int ShiftOpCs::_preparMatrices(AMesh *amesh,
                             int imesh,
                             MatrixCSGeneral& matu,
                             MatrixCRectangular& matw) const
{
  int ndim = amesh->getNDim();
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

cs* ShiftOpCs::_BuildSfromMap(std::map<std::pair<int, int>, double> &tab)
{
  std::map<std::pair<int, int>, double>::iterator it;

  cs* Striplet = cs_spalloc(0, 0, 1, 1, 1);

  it = tab.begin();
  while (it != tab.end())
  {
    int ip0 = it->first.first;
    int ip1 = it->first.second;
    if (!cs_entry(Striplet, ip0, ip1, it->second)) return nullptr;
    it++;
  }

  /* Optional printout */

  cs* S = cs_triplet(Striplet);
  if (S == (cs *) NULL) return nullptr;

  Striplet = cs_spfree(Striplet);

  return S;
}

/**
 * Calculate the private member "_SGrad" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param cova Description of the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance within the Model
 * @param nostat Structure for Non-Stationary parameters
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSGrad(AMesh *amesh,
                           const CovAniso& cova,
                           int igrf,
                           int icov,
                           ANoStat* nostat,
                           double tol)
{
  // Store the number of derivation parameters for the model as member
  _nModelGradParam = cova.getGradParamNumber();
  int number = _nModelGradParam * getSize();
  std::vector<std::map<std::pair<int, int>, double> > Mtab(number);

  int error = 1;
  int ndim = amesh->getNDim();
  int ncorner = amesh->getNApexPerMesh();
  int ngparam = _nModelGradParam;
  bool flag_nostat = (nostat != NULL);

  // Initialize the arrays

  VectorDouble matv(ncorner);
  MatrixCSGeneral hh = MatrixCSGeneral(ndim);
  MatrixCSGeneral matu = MatrixCSGeneral(ncorner);
  MatrixCSGeneral mat = MatrixCSGeneral(ncorner);
  MatrixCRectangular matw = MatrixCRectangular(ndim, ncorner);

  // Define the global HH matrix

  if (! flag_nostat)
  {
    _loadHHByApex(hh, cova, igrf, icov);
  }

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    debug_index(imesh + 1);
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
            _loadHHGradPerMesh(hh, amesh, cova, igp0, igparam, igrf, icov, imesh, nostat);
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
  if (error) for (int i = 0; i < (int) _SGrad.size(); i++)
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
 * @param cova  Description of the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance within the Model
 * @param nostat Structure for Non-Stationary parameters
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildS(AMesh *amesh,
                       const CovAniso& cova,
                       int igrf,
                       int icov,
                       ANoStat* nostat,
                       double tol)
{
  std::map<std::pair<int, int>, double> tab;

  int error = 1;
  int ndim = amesh->getNDim();
  int ncorner = amesh->getNApexPerMesh();
  bool flag_nostat = (nostat != NULL);

  // Initialize the arrays

  MatrixCSGeneral hh = MatrixCSGeneral(ndim);
  MatrixCSGeneral matu = MatrixCSGeneral(ncorner);
  MatrixCSGeneral mat = MatrixCSGeneral(ncorner);
  MatrixCRectangular matw = MatrixCRectangular(ndim, ncorner);

  // Define the global HH matrix

  if (!flag_nostat)
    _loadHHByApex(hh, cova, igrf, icov);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    debug_index(imesh + 1);
    double meshSize = amesh->getMeshSize(imesh);

    // Case of Euclidean geometry

    if (_preparMatrices(amesh, imesh, matu, matw))
    my_throw("Problem in matrix inversion");

    // Non stationary case

    if (flag_nostat)
    {
      if (nostat->isDefinedforAnisotropy(-1, icov))
        _loadHHPerMesh(hh, amesh, cova, igrf, icov, imesh, nostat);
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
 * Calculate the private member "_S" directly from the Mesh for Velocity
 * @param amesh Description of the Mesh (New class)
 * @param cova Description of the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance within the Model
 * @param nostat Structure for Non-Stationary parameters
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSVel(AMesh *amesh,
                          const CovAniso& cova,
                          int igrf,
                          int icov,
                          ANoStat* nostat,
                          double tol)
{
  std::map<std::pair<int, int>, double> tab;

  int error = 1;
  int ndim = amesh->getNDim();
  int ncorner = amesh->getNApexPerMesh();

  // Initialize the arrays

  VectorDouble vel(2);
  VectorDouble matv(ncorner);
  MatrixCSGeneral matu = MatrixCSGeneral(ncorner);
  MatrixCRectangular matw = MatrixCRectangular(ndim, ncorner);

  // Define the global HH matrix

  _loadAux(vel, igrf, icov, CONS_VELOCITY);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    debug_index(imesh + 1);
    double meshSize = amesh->getMeshSize(imesh);

    // Non stationary case

    if (nostat->isDefined(igrf, icov, CONS_VELOCITY, -1, -1))
      _loadAuxPerMesh(vel, amesh, igrf, icov, CONS_VELOCITY, imesh, nostat);

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
 * Calculate the private member "_S" directly from the Mesh
 * @param amesh Description of the Mesh (New class)
 * @param cova Description of the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance within the Model
 * @param nostat Structure for Non-Stationary parameters
 * @param tol Tolerance beyond which elements are not stored in S matrix
 * @return Error return code
 */
int ShiftOpCs::_buildSSphere(AMesh *amesh,
                             const CovAniso& cova,
                             int igrf,
                             int icov,
                             ANoStat* nostat,
                             double tol)
{
  std::map<std::pair<int, int>, double> tab;
  double coeff[3][2];

  int error = 1;
  int ndim = amesh->getNDim();
  int ncorner = amesh->getNApexPerMesh();
  bool flag_nostat = (nostat != NULL);
  bool flag_vel = false;
  if (flag_nostat)
    flag_vel = nostat->isDefined(igrf, icov, CONS_VELOCITY, -1, -1);

  // Initialize the arrays

  VectorDouble srot(2), center[3], axe1(3), axe2(3), vel(3), matv(ncorner);
  MatrixCSGeneral hh = MatrixCSGeneral(ndim);
  MatrixCSGeneral matu = MatrixCSGeneral(ncorner);
  MatrixCSGeneral mat = MatrixCSGeneral(ncorner);
  MatrixCRectangular matw = MatrixCRectangular(ndim, ncorner);

  // Define the global HH matrix

  if (!flag_nostat)
  {
    _loadHHByApex(hh, cova, igrf, icov);
    _loadAux(srot, igrf, icov, CONS_SPHEROT);
    _loadAux(vel, igrf, icov, CONS_VELOCITY);
  }

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {
    debug_index(imesh + 1);

    // Non stationary case

    if (flag_nostat)
    {
      if (nostat->isDefinedforAnisotropy(igrf, icov))
        _loadHHPerMesh(hh, amesh, cova, igrf, icov, imesh, nostat);
      if (nostat->isDefined(igrf, icov, CONS_SPHEROT, -1, -1))
        _loadAuxPerMesh(srot, amesh, igrf, icov, CONS_SPHEROT, imesh, nostat);
      if (nostat->isDefined(igrf, icov, CONS_VELOCITY, -1, -1))
        _loadAuxPerMesh(vel, amesh, igrf, icov, CONS_VELOCITY, imesh, nostat);
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

    if (flag_vel)
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
        double vald = (flag_vel) ? matv[j1] : mat.getValue(j0, j1);
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
 * Calculate _TildeC directly from the Mesh (Dimension: _NApices)
 * @param amesh Description of the Mesh (New class)
 * @param units Array of sizes for all meshes
 * @return Error return code
 */
int ShiftOpCs::_buildTildeC(AMesh *amesh, const VectorDouble& units)
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
 * @param cova Structure containing the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance
 * @param nostat Structure for Non-Stationary parameters
 * @return
 */
void ShiftOpCs::_buildLambda(AMesh *amesh,
                             const CovAniso& cova,
                             int igrf,
                             int icov,
                             ANoStat* nostat)
{
  double sqdeth = 0.;

  int ndim = amesh->getNDim();
  int nvertex = amesh->getNApices();
  bool flag_nostat = (nostat != NULL);

  /* Core allocation */

  _Lambda.clear();

  MatrixCSGeneral hh = MatrixCSGeneral(ndim);
  if (!flag_nostat)
  {
    _loadHHByApex(hh, cova, icov);
    sqdeth = sqrt(matrix_determinant(ndim, hh.getValues().data()));
  }

  /* Fill the array */

  double sill = cova.getSill(0, 0);
  for (int ip = 0; ip < nvertex; ip++)
  {
    if (flag_nostat)
    {
      if (nostat->isDefinedforAnisotropy(igrf, icov))
      {
        _loadHHByApex(hh, cova, igrf, icov, ip, nostat);
        sqdeth = sqrt(matrix_determinant(ndim, hh.getValues().data()));
      }
      if (nostat->isDefined(igrf, icov, CONS_SILL, -1, -1))
        sill = nostat->getValue(igrf, icov, CONS_SILL, -1, -1, 0, ip);
    }
    _Lambda.push_back(sqrt((_TildeC[ip]) / (sqdeth * sill)));
  }
}

/**
 * Construct the _Lambda vector (Dimension: _NApices)
 * @param amesh Description of the Mesh (New class)
 * @param cova Structure containing the Covariance
 * @param igrf Rank of the GRF
 * @param icov Rank of the covariance
 * @param nostat Structure for Non-Stationary parameters
 * @return
 */
bool ShiftOpCs::_buildLambdaGrad(AMesh *amesh,
                             CovAniso& cova,
                             int igrf,
                             int icov,
                             ANoStat* nostat)
{


  int ndim = amesh->getNDim();
  int nvertex = amesh->getNApices();


  /* Core allocation */

  _LambdaGrad.clear();
  VectorDouble temp(nvertex);
  for(int i = 0; i< ndim; i++)
  {
    _LambdaGrad.push_back(temp);
  }


  /* Fill the array */


  for (int ip = 0; ip < nvertex; ip++)
  {

     _updateCova(cova, igrf, icov, ip, ndim, nostat);


    for(int idim = 0;idim < ndim ;idim++)
    {
      _LambdaGrad[idim][ip] = - _Lambda[ip] / (2. * cova.getScale(idim));
    }

   // dsqdeth = sqrt(matrix_determinant(ndim, hh.getValues().data()));
   // _Lambda.push_back(sqrt((_TildeC[ip]) / (dsqdeth * sill)));
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
void ShiftOpCs::_projectMesh(AMesh *amesh,
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
  double phi = srot[1] * GV_PI / 180.;
  double theta = srot[0] * GV_PI / 180.;
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
