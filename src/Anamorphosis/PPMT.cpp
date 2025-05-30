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
#include "geoslib_define.h"

#include "Anamorphosis/PPMT.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Db/Db.hpp"
#include "Stats/Classical.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/MathFunc.hpp"

#include <math.h>

PPMT::PPMT(int ndir,
           bool flagPreprocessing,
           const EDirGen& methodDir,
           const EGaussInv& methodTrans,
           int nbpoly,
           double alpha)
    : AStringable(),
      _niter(0),
      _ndir(ndir),
      _nbpoly(nbpoly),
      _alpha(alpha),
      _methodDir(methodDir),
      _methodTrans(methodTrans),
      _flagPreprocessing(flagPreprocessing),
      _isFitted(false),
      _ndim(0),
      _serieAngle(),
      _serieScore(),
      _dirmat(nullptr),
      _anams(),
      _initAnams(),
      _initSphering(nullptr)
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _niter(m._niter),
      _ndir(m._ndir),
      _nbpoly(m._nbpoly),
      _alpha(m._alpha),
      _methodDir(m._methodDir),
      _methodTrans(m._methodTrans),
      _flagPreprocessing(m._flagPreprocessing),
      _isFitted(m._isFitted),
      _ndim(m._ndim),
      _serieAngle(m._serieAngle),
      _serieScore(m._serieScore),
      _dirmat(m._dirmat),
      _anams(m._anams),
      _initAnams(m._initAnams),
      _initSphering(m._initSphering)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _niter = m._niter;
    _ndir = m._ndir;
    _nbpoly = m._nbpoly;
    _ndim = m._ndim;
    _alpha = m._alpha;
    _methodDir = m._methodDir;
    _methodTrans = m._methodTrans;
    _flagPreprocessing = m._flagPreprocessing;
    _isFitted = m._isFitted;
    _ndim = m._ndim;
    _serieAngle = m._serieAngle;
    _serieScore = m._serieScore;
    _dirmat = m._dirmat;
    _anams = m._anams;
    _initAnams = m._initAnams;
    _initSphering = m._initSphering;
  }
  return *this;
}

PPMT::~PPMT()
{
  delete _dirmat;
  if (! _anams.empty())
  {
    for (int i = 0; i < (int) _anams.size(); i++)
      delete _anams[i];
  }
  if (! _initAnams.empty())
  {
    for (int i = 0; i < (int) _initAnams.size(); i++)
      delete _initAnams[i];
  }
  delete _initSphering;
}

/**
 * Create the Multivariate Gaussian anamorphosis
 * @param ndir Number of Directions to be tested at each iteration
 * @param flagPreprocessing True for pre-processing (Normal Score and Sphering)
 * @param methodDir Method for Direction Generation
 * @param methodTrans Method for Gaussian Transformation
 * @param nbpoly Number of Polynomial (only used for "hermite" transformation)
 * @param alpha Distance exponent
 * @return An instance of PPMT class
 */
PPMT* PPMT::create(int ndir,
                   bool flagPreprocessing,
                   const EDirGen& methodDir,
                   const EGaussInv& methodTrans,
                   int nbpoly,
                   double alpha)
{
  return new PPMT(ndir, flagPreprocessing, methodDir, methodTrans, nbpoly, alpha);
}

String PPMT::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);

  std::stringstream sstr;

  sstr << toTitle(1, "PPMT");
  if (_flagPreprocessing)
    sstr << "- Initial Anamorphosis per component and Sphering" << std::endl;
  if (getMethodDir() == EDirGen::VDC)
    sstr << "- Using Van der Corput method for generating Directions" << std::endl;
  else
    sstr << "- Using Uniform method for generating Directions" << std::endl;
  sstr << "- Number of iterations = " << getNiter() << std::endl;
  sstr << "- Number of Directions = " << getNdir() << std::endl;
  if (getMethodTrans() == EGaussInv::HMT)
    sstr << "- Number of Hermite Polynomials = " << getNbpoly() << std::endl;
  sstr << "- Exponent value       = " << getAlpha() << std::endl;

  sstr << std::endl;
  if (isFitted())
    sstr << "Fitting has been performed" << std::endl;
  else
    sstr << "Fitting has not been performed yet" << std::endl;

  return sstr.str();
}

double PPMT::_gaussianizeForward(double Yi,
                                 int rank,
                                 const AnamHermite *anam,
                                 const VectorDouble &N0)
{
  double theo = 0.;
  if (anam != nullptr)
    theo = anam->rawToTransformValue(Yi);
  else
    theo = N0[rank];
  return (theo - Yi);
}

/**
 * Calculate the Inverse normal score transform. This is only available for Hermite
 * @param Yi   Input value
 * @param anam Anamorphosis
 * @return The back-anamorphosed value
 */
double PPMT::_gaussianizeBackward(double Yi, const AnamHermite *anam)
{
  double theo = anam->transformToRawValue(Yi);
  return (theo - Yi);
}

void PPMT::_initGaussianizeForward(AMatrix* Y)
{
  int ncol = getNdim();

  for (int icol = 0; icol < ncol; icol++)
  {
    VectorDouble Zvec = Y->getColumn(icol);
    VectorDouble Yvec = _initAnams[icol]->rawToGaussianVector(Zvec);
    Y->setColumn(icol, Yvec);
  }
}

void PPMT::_initGaussianizeBackward(AMatrix* Y)
{
  int ncol = getNdim();

  for (int icol = 0; icol < ncol; icol++)
  {
    VectorDouble Zvec = Y->getColumn(icol);
    VectorDouble Yvec = _initAnams[icol]->gaussianToRawVector(Zvec);
    Y->setColumn(icol, Yvec);
  }
}

double PPMT::_getGaussianDistance(const VectorDouble &Yi,
                                  const VectorInt &Ri,
                                  const VectorDouble &N0) const
{
  int np = (int) Yi.size();

  double di = 0.;
  for (int ip = 0; ip < np; ip++)
  {
    double value = _gaussianizeForward(Yi[ip], Ri[ip], nullptr, N0);
    di += pow(ABS(value), getAlpha());
  }
  di /= (double) np;
  return di;
}

void PPMT::_iterationFit(AMatrix *Y, const VectorDouble& N0)
{
  int np = Y->getNRows();

  // Initialization

  VectorDouble Y0(np, TEST);
  VectorDouble Yi(np, TEST);
  VectorInt    R0(np, ITEST);
  int    idmax = -1;
  double ddmax = -1.e30;

  // Loop on directions

  AnamHermite* anam = nullptr;
  const bool flagHermite = (getMethodTrans() == EGaussInv::HMT);
  if (flagHermite) anam = new AnamHermite(getNbpoly());

  for (int id = 0; id < getNdir(); id++)
  {
    // Projection of the data set on the target direction
    _projectOnDirection(Y, id, Yi);

    // Preparing the Normal scoring for the target distance
    VectorInt Ri = VH::sortRanks(Yi, true, np);

    // Calculate the distance on the projected axis
    double di = _getGaussianDistance(Yi, Ri, N0);

    // Keep the Direction corresponding to the largest Gaussian scPreore
    if (ddmax < di)
    {
      idmax = id;
      ddmax = di;
      Y0 = Yi;
      R0 = Ri;
    }
  }

  if (flagHermite) anam->fitFromArray(Y0);
  _shiftForward(Y, idmax, anam, Y0, R0, N0);

  // Returning arguments
  _serieAngle.push_back(idmax);
  _serieScore.push_back(ddmax);
  if (flagHermite) _anams.push_back(anam);
}

void PPMT::_shiftForward(AMatrix *Y,
                         int id,
                         const AnamHermite *anam,
                         const VectorDouble &Y0,
                         const VectorInt    &R0,
                         const VectorDouble &N0) const
{
  int np   = Y->getNRows();
  int ndim = getNdim();

  for (int ip = 0; ip < np; ip++)
  {
    double scale = _gaussianizeForward(Y0[ip], R0[ip], anam, N0);
    for (int idim = 0; idim < ndim; idim++)
    {
      double value = Y->getValue(ip, idim) + scale * _dirmat->getValue(id, idim);
      Y->setValue(ip, idim, value);
    }
  }
}

void PPMT::_shiftBackward(AMatrix *Y,
                          int id,
                          const AnamHermite *anam,
                          const VectorDouble &Y0) const
{
  int np   = Y->getNRows();
  int ndim = getNdim();

  for (int ip = 0; ip < np; ip++)
  {
    double scale = _gaussianizeBackward(Y0[ip], anam);
    for (int idim = 0; idim < ndim; idim++)
    {
      double value = Y->getValue(ip, idim) - scale * _dirmat->getValue(id, idim);
      Y->setValue(ip, idim, value);
    }
  }
}

/**
 * Project the data set on the direction rank 'id'
 * @param Y   Matrix containing the data set
 * @param id  Rank of the target direction
 * @param Y0  Vector of projected coordinates
 */
void PPMT::_projectOnDirection(const AMatrix* Y, int id, VectorDouble& Y0)
{
  int np   = Y->getNRows();
  int ndim = getNdim();

  for (int ip = 0; ip < np; ip++)
  {
    double value = 0.;
    for (int idim = 0; idim < ndim; idim++)
      value += Y->getValue(ip, idim) * _dirmat->getValue(id, idim);
    Y0[ip] = value;
  }
}

void PPMT::_iterationForward(AMatrix *Y, const VectorDouble& N0, int iter)
{
  int np    = Y->getNRows();
  int idmax = _serieAngle[iter];
  VectorDouble Y0(np, TEST);

  // Projection of the data set on the optimal direction
  _projectOnDirection(Y, idmax, Y0);

  // Preparing the Normal scoring for the target distance
  VectorInt R0 = VH::sortRanks(Y0, true, np);

  // Forward Shift
  _shiftForward(Y, idmax, _anams[iter], Y0, R0, N0);
}

void PPMT::_iterationBackward(AMatrix *Y, const VectorDouble& N0, int iter)
{
  DECLARE_UNUSED(N0);
  int np    = Y->getNRows();
  int idmax = _serieAngle[iter];
  VectorDouble Y0(np, TEST);

  // Forward Shift
  _shiftBackward(Y, idmax, _anams[iter], Y0);

  // Projection of the data set on the optimal direction
  _projectOnDirection(Y, idmax, Y0);
}

void PPMT::_generateAllDirections()
{
  int ndir = getNdir();

  MatrixDense* Umat;
  if (getMethodDir() == EDirGen::VDC)
  {
    Umat = vanDerCorput(ndir, _ndim);
  }
  else
  {
    VectorDouble X = VH::simulateUniform(ndir * _ndim);
    Umat = MatrixDense::createFromVD(X, ndir, _ndim);
  }
  _dirmat = GeometryHelper::getDirectionsInRn(Umat);
  delete Umat;
}

void PPMT::_fitInitHermite(AMatrix* Y)
{
  int ncol = getNdim();

  for (int icol = 0; icol < ncol; icol++)
  {
    VectorDouble Yvec = Y->getColumn(icol);
    AnamHermite* anam = new AnamHermite(getNbpoly());
    anam->fitFromArray(Yvec);
    _initAnams.push_back(anam);
  }
}

int PPMT::fitFromMatrix(AMatrix *Y, int niter, bool verbose)
{
  if (Y == nullptr)
  {
    messerr("Input Argument 'Y' (matrix) should be provided. Nothing is done");
  }
  _ndim  = Y->getNCols();
  _niter = niter;

  // Clearing the storage
  _serieAngle.clear();
  _serieScore.clear();
  _initAnams.clear();
  _anams.clear();

  // Creating the directions
  _generateAllDirections();

  int np = Y->getNRows();
  VectorDouble sequence = VH::sequence(1., np, 1., 1. + np);
  VectorDouble N0 = VH::qnormVec(sequence);

  // Optional Pre-processing
  if (_flagPreprocessing)
  {
    if (verbose) message("Pre-processing:\n");

    // Anamorphosis transform of each component of the input vector
    if (verbose) message("- Normal scoring each component\n");
    _fitInitHermite(Y);
    _initGaussianizeForward(Y);

    // Sphering
    if (verbose) message("- Sphering\n");
    _initSphering = sphering(Y);
    Y->prodMatInPlace(_initSphering);
  }

  // Loop on the iterations

  if (verbose) message("\nLoop on iterations to find best direction:\n");
  for (int iter = 0; iter < niter; iter++)
  {
    _iterationFit(Y, N0);

    if (verbose)
    {
      message("Iteration %3d/%3d: Score = %lf\n",iter+1,niter, _serieScore[iter]);
    }
  }

  // Set the flag saying that the PPMT model has been fitted correctly
  _isFitted = true;

  return 0;
}

int PPMT::fit(Db *db,
              const VectorString &names,
              bool flagStoreInDb,
              int niter,
              bool verbose,
              const NamingConvention &namconv)
{
  VectorString exp_names = db->expandNameList(names);
  MatrixDense Y = db->getColumnsAsMatrix(exp_names, true);
  if (Y.empty())
  {
    messerr("This Multivariate Transform requires several variables to be defined");
    return 1;
  }

  // Fitting the PPMT model
  if (fitFromMatrix(&Y, niter, verbose)) return 1;

  // Add the newly created information in the Db (optional)

  if (flagStoreInDb)
  {
    int iptr = db->addColumns(Y.getValues(), String(), ELoc::UNKNOWN, 0, true);
    namconv.setNamesAndLocators(exp_names, db, iptr);
  }
  return 0;
}

int PPMT::rawToGaussian(Db *db,
                        const VectorString &names,
                        int niter,
                        const NamingConvention &namconv)
{
  // Extract the relevant information

  if (db == nullptr)
  {
    messerr("The argument 'db' must be provided");
    return 1;
  }
  VectorString exp_names = db->expandNameList(names);
  MatrixDense Y = db->getColumnsAsMatrix(exp_names, true);
  if (Y.empty())
  {
    messerr("This Multivariate Transform requires several variables to be defined");
    return 1;
  }
  if (! isFitted())
  {
    messerr("You must Fit PPMT beforehand");
    return 1;
  }
  int np = Y.getNRows();
  if (niter <= 0) niter = getNiter();
  niter = MIN(niter, getNiter());

  VectorDouble sequence = VH::sequence(1., np, 1., 1. + np);
  VectorDouble N0 = VH::qnormVec(sequence);

  // Pre-processing
  if (_flagPreprocessing)
  {
    _initGaussianizeForward(&Y);
    Y.prodMatInPlace(_initSphering);
  }

  // Loop on the iterations
  for (int iter = 0; iter < niter; iter++)
    _iterationForward(&Y, N0, iter);

  // Add the newly created information in the Db
  int iptr = db->addColumns(Y.getValues(), String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(exp_names, db, iptr);

  return 0;
}

int PPMT::gaussianToRaw(Db *db,
                        const VectorString &names,
                        int niter,
                        const NamingConvention &namconv)
{
  // Extract the relevant information

  if (db == nullptr)
  {
    messerr("The argument 'db' must be provided");
    return 1;
  }
  VectorString exp_names = db->expandNameList(names);
  MatrixDense Y = db->getColumnsAsMatrix(exp_names, true);
  if (Y.empty())
  {
    messerr("This Multivariate Back-Transform requires several variables to be defined");
    return 1;
  }
  if (getMethodTrans() != EGaussInv::HMT)
  {
    messerr("The PPMT back-trasform is only available when methodTrans = 'hermite'");
    return 1;
  }
  if (! isFitted())
  {
    messerr("You must Fit PPMT beforehand");
    return 1;
  }
  int np = Y.getNRows();
  if (niter <= 0) niter = getNiter();
  niter = MIN(niter, getNiter());

  VectorDouble sequence = VH::sequence(1., np, 1., 1. + np);
  VectorDouble N0 = VH::qnormVec(sequence);

  // Loop on the iterations (reverse order)
  for (int iter = niter-1; iter >= 0; iter--)
    _iterationBackward(&Y, N0, iter);

  // Post-processing
  if (_flagPreprocessing)
  {
    AMatrix* backSphering = _initSphering->transpose();
    Y.prodMatInPlace(backSphering);

    _initGaussianizeBackward(&Y);
  }

  // Add the newly created information in the Db
  int iptr = db->addColumns(Y.getValues(), String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(exp_names, db, iptr);

  return 0;
}

VectorDouble PPMT::getSerieScore(bool flagLog) const
{
  VectorDouble vec;
  for (int iter = 0; iter < getNiter(); iter++)
  {
    double value = _serieScore[iter];
    if (flagLog) value = log(value);
    vec.push_back(value);
  }
  return vec;
}
