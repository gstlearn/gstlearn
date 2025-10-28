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
#include "MultiLayers/MGibbs.hpp"

#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Core/Keypair.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovFactory.hpp"
#include "Db/Db.hpp"
#include "Enum/ELoadBy.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Space/SpaceSN.hpp"
#include "geoslib_define.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"
#include <cmath>
#include <cstring>

/* Global symbols for SPDE */

#define NBLIN_TERMS 10

/*! \cond */
#define VT_FREE  1
#define VT_GIBBS 2
#define VT_HARD  4

#define VT_INPUT  8
#define VT_OUTPUT 16
#define VT_OTHER  32

#define CASE_KRIGING  1
#define CASE_SIMULATE 2

#define GWORK(ilayer, iech)      (gwork[(ilayer) * _nechout + (iech)])
#define YDAT(nech, ilayer, iech) (ydat[(ilayer) * nech + (iech)])
#define YMEAN(ilayer, iech)      (ymean[(ilayer) * _nechc + (iech)])

namespace gstlrn
{
typedef struct
{
  std::vector<SPDE_Matelem> Matelems;
  Id ndata;                            /* Number of active data */
  VectorInt ndata1;                    /* Number of data per variable (icov=0) */
  VectorInt ntarget1;                  /* Number of target per variable (icov=0) */
  VectorDouble Csill;                  /* Array of LU of sill matrices */
  std::vector<MatrixSparse*> Bnugget;  /* Sparse matrices for nugget effect (nvs2) */
  std::vector<MatrixSparse*> BheteroD; /* Sparse matrices for heterotopy (_nvar)*/
  std::vector<MatrixSparse*> BheteroT; /* Sparse matrices for heterotopy (_nvar)*/
} SPDE_Environ;

typedef struct
{
  bool flag_est;     /* Perform Estimation */
  bool flag_std;     /* Perform Standard deviation */
  Id flag_case;      /* Perform: matrices(0), est(1) or simu(2) */
  bool flag_filnug;  /* Filtering the Nugget Effect */
  bool flag_several; /* Perform Kriging in iterative mode */
} SPDE_Decision;

typedef struct
{
  Id flag_sphere;
  double sqdeth;
  double correc;
  double R;
  VectorDouble blin;
  MatrixSquare hh;
  VectorDouble vv;
  VectorDouble srot;
} SPDE_Calcul;

/*! \endcond */
static bool DEBUG           = true;
static double FACDIM[]      = {0., 1., 2., 6.};
static Id SPDE_CURRENT_ICOV = 0;
static SPDE_Environ S_ENV;
static SPDE_Decision S_DECIDE;
static String string_encode;
static SPDE_Calcul Calcul;

MGibbs::MGibbs(Db* dbin, Db* dbout, Model* model)
  : AStringable()
  , _hasValidEnvironment(false)
  , _dbin(dbin)
  , _dbout(dbout)
  , _model(model)
  , _ndim(0)
  , _nvar(0)
  , _nechin(0)
  , _nechout(0)
  , _nechc(0)
  , _nlayer(0)
  , _niter(0)
  , _seed(0)
  , _nbsimu(0)
  , _icolPinch(-1)
  , _flagED(false)
  , _flagDrift(false)
  , _flagCE(false)
  , _flagCStd(false)
  , _verbose(false)
{
  _hasValidEnvironment = _checkArguments();
}

MGibbs::MGibbs(const MGibbs& m)
  : AStringable(m)
  , _hasValidEnvironment(m._hasValidEnvironment)
  , _dbin(m._dbin)
  , _dbout(m._dbout)
  , _model(m._model)
  , _ndim(m._ndim)
  , _nvar(m._nvar)
  , _nechin(m._nechin)
  , _nechout(m._nechout)
  , _nechc(m._nechc)
  , _nlayer(m._nlayer)
  , _niter(m._niter)
  , _seed(m._seed)
  , _nbsimu(m._nbsimu)
  , _icolPinch(m._icolPinch)
  , _flagED(m._flagED)
  , _flagDrift(m._flagDrift)
  , _flagCE(m._flagCE)
  , _flagCStd(m._flagCStd)
  , _verbose(m._verbose)
{
}

MGibbs& MGibbs::operator=(const MGibbs& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _hasValidEnvironment = m._hasValidEnvironment;
    _dbin                = m._dbin;
    _dbout               = m._dbout;
    _model               = m._model;
    _ndim                = m._ndim;
    _nvar                = m._nvar;
    _nechin              = m._nechin;
    _nechout             = m._nechout;
    _nechc               = m._nechc;
    _nlayer              = m._nlayer;
    _niter               = m._niter;
    _seed                = m._seed;
    _nbsimu              = m._nbsimu;
    _icolPinch           = m._icolPinch;
    _flagED              = m._flagED;
    _flagDrift           = m._flagDrift;
    _flagCE              = m._flagCE;
    _flagCStd            = m._flagCStd;
    _verbose             = m._verbose;
  }
  return *this;
}

MGibbs::~MGibbs()
{
}

/// Interface to AStringable
String MGibbs::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  return String();
}

/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the Cholesky decomposition
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  QC   QChol structure to be finalized
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
Id MGibbs::_qcholCholesky(bool verbose, QChol* QC)
{
  /* Check that the Q matrix has already been defined */

  if (QC->Q == nullptr) return (1);

  /* Cholesky decomposition is only valid for square matric */

  if (QC->Q->getNRows() != QC->Q->getNCols())
  {
    messerr("You wish to perform a Cholesky Decomposition of a Matrix");
    messerr("which is not square: %d x %d", QC->Q->getNRows(), QC->Q->getNCols());
    messerr("This must be an error");
    return (1);
  }

  if (verbose) message("  Cholesky Decomposition... ");

  // Perform the Cholesky decomposition (new style)

  if (QC->chol == nullptr)
  {
    QC->chol = new CholeskySparse(*QC->Q);
    if (QC->chol == nullptr)
    {
      messerr("Error in Cholesky decompostion (new version)");
      messerr("Error in Cholesky decompostion (new version)");
      return 1;
    }
  }

  if (verbose) message("Finished\n");

  return (0);
}

/****************************************************************************/
/*!
 **  Inversion using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 ** \param[in,out]  xcr  Current vector
 ** \param[in]  rhs      Current R.H.S. vector
 **
 *****************************************************************************/
void MGibbs::_cholInvert(QChol* qctt,
                         VectorDouble& xcr,
                         const VectorDouble& rhs)
{
  auto n = qctt->Q->getNCols();

  VectorDouble rhsVD(n);
  for (Id i = 0; i < n; i++) rhsVD[i] = rhs[i];
  VectorDouble xcrVD = qctt->chol->solveX(rhsVD);
  for (Id i = 0; i < n; i++) xcr[i] = xcrVD[i];
}

/****************************************************************************/
/*!
 **  Simulate using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 **
 ** \param[out] simu     Simulated array
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void MGibbs::_cholSimulate(QChol* qctt,
                           VectorDouble& simu,
                           const VectorDouble& work)
{
  auto n = qctt->Q->getNCols();

  VectorDouble simuVD(n);
  VectorDouble workVD(n);
  for (Id i = 0; i < n; i++) workVD[i] = work[i];
  (void)qctt->chol->addSimulateToDest(workVD, simuVD);
  for (Id i = 0; i < n; i++) simu[i] = simuVD[i];
}

/****************************************************************************/
/*!
 **  Returns the index of a pair of variable ranks within the triangle
 **
 ** \return Absolute index
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks The calling function does not have to bother of the relative
 ** \remarks order between 'ivar' and 'jvar'
 **
 *****************************************************************************/
Id MGibbs::_getRank(Id ivar, Id jvar)
{
  if (jvar > ivar) return (jvar * (jvar + 1) / 2 + ivar);
  return (ivar * (ivar + 1) / 2 + jvar);
}

/****************************************************************************/
/*!
 **  Print the contents of one SP_Mat structure
 **
 ** \param[in] icov    Rank of the covariance
 **
 *****************************************************************************/
void MGibbs::_printMatelem(Id icov)
{
  static const char* NOK[] = {"OFF", "ON"};

  const SPDE_Matelem& Matelem = _getCurrentMatelem(icov);

  mestitle(1, "Contents of Matelem structure #%d", icov + 1);
  message("S is defined:      %s\n", NOK[Matelem.S != NULL]);
  message("Aproj is defined:  %s\n", NOK[Matelem.Aproj != NULL]);
  message("QC is defined:     %s\n", NOK[Matelem.QC != NULL]);
  message("QCov are defined:  %s\n", NOK[!Matelem.QCov.empty()]);
  message("Lambda is defined: %s\n", NOK[!Matelem.Lambda.empty()]);
  message("QCtt is defined:   %s\n", NOK[Matelem.QCtt != NULL]);
  message("QCtd is defined:   %s\n", NOK[Matelem.QCtd != NULL]);
  message("s_cheb is defined: %s\n", NOK[Matelem.s_cheb != NULL]);
  message("s_mesh is defined: %s\n", NOK[Matelem.amesh != NULL]);
}

/****************************************************************************/
/*!
 **  Get the pointer to the current SPDE_Matelem structure
 **
 ** \param[in] icov    Rank of the target Covariance (or -1)typedef struct
 **
 *****************************************************************************/
SPDE_Matelem& MGibbs::_getCurrentMatelem(Id icov)
{
  if (icov < 0)
    return (S_ENV.Matelems[SPDE_CURRENT_ICOV]);
  return (S_ENV.Matelems[icov]);
}

/****************************************************************************/
/*!
 **  Update a string to include the rank of the current GRF and Covariance
 **
 ** \param[in]  flag_icov  To add current COV
 ** \param[in]  rank       Rank of the highlight (see mestitle or -1)
 ** \param[in]  title      Input title
 **
 *****************************************************************************/
void MGibbs::_setTitle(bool flag_icov, Id rank, const char* title)
{
  if (flag_icov)
  {
    (void)gslSPrintf(string_encode, "(COV:%d) %s", SPDE_CURRENT_ICOV + 1, title);
  }
  else
  {
    (void)gslSPrintf(string_encode, "%s", title);
  }

  if (rank >= 0)
    mestitle(rank, string_encode.data());
  else
  {
    (void)gslSPrintf(string_encode, "%s\n",
                     string_encode.data());
    message(string_encode.data());
  }
}

/****************************************************************************/
/*!
 **  Returns the pointer to structure containing the Nugget Effect (or NULL)
 **
 *****************************************************************************/
CovAniso* MGibbs::_getNugget(void) const
{
  for (Id is = 0; is < _model->getNCov(); is++)
  {
    CovAniso* cova = _model->getCovAniso(is);
    if (cova->getType() == ECov::NUGGET) return (cova);
  }
  return (nullptr);
}

/****************************************************************************/
/*!
 **  Returns the pointer to the structure
 **
 *****************************************************************************/
CovAniso* MGibbs::_getCurrentCova(void) const

{
  Id is0 = SPDE_CURRENT_ICOV;

  for (Id icov = 0, jcov = 0; icov < _model->getNCov(); icov++)
  {
    CovAniso* cova = _model->getCovAniso(icov);
    if (cova->getType() == ECov::NUGGET) continue;
    if (is0 == jcov) return (cova);
    jcov++;
  }
  return (nullptr);
}

/****************************************************************************/
/*!
 **  Return if a nugget effect component must be filtered
 **
 *****************************************************************************/
bool MGibbs::_isFilterNugget(void)
{
  return (S_DECIDE.flag_filnug && S_DECIDE.flag_case == CASE_KRIGING);
}

/****************************************************************************/
/*!
 **  Defines if a nugget effect component must be filtered
 **
 ** \param[in]  flag_filnug  Flag to define if a nugget effect must be filtered
 **
 *****************************************************************************/
void MGibbs::_setFilterNugget(Id flag_filnug)
{
  if (DEBUG) _setTitle(false, -1, "(DEBUG) Set 'filnug'");
  S_DECIDE.flag_filnug = flag_filnug;
}

/****************************************************************************/
/*!
 **  Get the value of the Inverse of the sill for a given covariance and
 **  a pair of variables
 **
 *****************************************************************************/
double MGibbs::_getIsill(Id icov, Id ivar, Id jvar)
{
  const SPDE_Matelem& Maticov = _getCurrentMatelem(icov);
  return Maticov.Isill.getValue(ivar, jvar);
}

/****************************************************************************/
/*!
 **  Encode the status of the variable
 **
 ** \param[in]  auth    Status option
 **
 *****************************************************************************/
void MGibbs::_printStatus(Id auth)

{
  if (auth & VT_OTHER) message("OTHER ");
  if (auth & VT_FREE) message("FREE ");
  if (auth & VT_GIBBS) message("GIBBS ");
  if (auth & VT_HARD) message("DATA ");
  if (auth & VT_INPUT) message("INPUT ");
  if (auth & VT_OUTPUT) message("OUTPUT ");
}

/****************************************************************************/
/*!
 **  Print information on the Filter
 **
 ** \param[in]  title   Optional title
 ** \param[in]  auth    Filter option
 **
 *****************************************************************************/
void MGibbs::_printQcholFilter(const char* title, Id auth)
{
  message("%s = ", title);
  _printStatus(auth);
  message("\n");
}

/****************************************************************************/
/*!
 **  Construct the sparse matrix Q from another sparse matrix
 **
 ** \return The Q structure or NULL
 **
 ** \param[in]  Q_in      Input sparse matrix
 ** \param[in]  row_auth  Specification for rows extraction
 ** \param[in]  col_auth  Specification for columns extraction
 **
 ** \remarks The Cholesky decomposition is performed (if possible)
 **
 *****************************************************************************/
MatrixSparse* MGibbs::_extractQfromQ(MatrixSparse* Q_in, Id row_auth, Id col_auth)
{
  DECLARE_UNUSED(row_auth, col_auth);

  /* Initializations */

  Id n_in = Q_in->getNCols();
  VectorInt rank_rows(n_in);
  VectorInt rank_cols(n_in);

  /* Fill the index vectors */

  for (Id i = 0; i < n_in; i++)
  {
    rank_rows[i] = i;
    rank_cols[i] = i;
  }

  /* Extract the submatrix */

  return Q_in->extractSubmatrixByRanks(rank_rows, rank_cols);
}

/****************************************************************************/
/*!
 **  Manage the Sparse matrix and its cholesky decomposition
 **
 ** \return  Pointer to the QChol structure
 **
 ** \param[in]  mode       Management mode
 **                         1 : Allocation
 **                        -1 : Total deallocation (NULL is returned)
 ** \param[in]  QC         Pointer to the QChol structure (when mode < 0)
 **
 *****************************************************************************/
QChol* MGibbs::_manageQchol(Id mode, QChol* QC)
{

  /* Dispatch */

  switch (mode)
  {
    case 1: /* Allocation */
      QC       = new QChol;
      QC->Q    = nullptr;
      QC->chol = nullptr;
      break;

    case -1: /* Total deallocation */
      if (QC == nullptr) return (QC);
      delete QC->Q;
      delete QC->chol;
      delete QC;
      QC = nullptr;
      break;
  }
  return (QC);
}

/****************************************************************************/
/*!
 **  Construct the QChol sub-structure from another QChol structure
 **
 ** \return The QChol structure or NULL
 **
 ** \param[in]  title     Name of the QChol item
 ** \param[in]  QC_in     Input QChol structure
 ** \param[in]  row_auth  Specification for rows extraction
 ** \param[in]  col_auth  Specification for columns extraction
 **
 *****************************************************************************/
QChol* MGibbs::_extractQCfromQ(const char* title,
                               QChol* QC_in,
                               Id row_auth,
                               Id col_auth) const
{
  Id error  = 1;
  QChol* QC = _manageQchol(1, nullptr);

  /* Extract the submatrix */

  QC->Q = _extractQfromQ(QC_in->Q, row_auth, col_auth);
  if (QC->Q == nullptr) goto label_end;

  /* Optional printout */

  if (_verbose)
  {
    message("Extracting a part of Q for '%s'\n", title);
    _printQcholFilter("- Row authorization code   ", row_auth);
    _printQcholFilter("- Column authorization code", col_auth);
  }

  /* Set the error return code */

  error = 0;

label_end:
  if (error) QC = _manageQchol(-1, QC);
  return (QC);
}

/****************************************************************************/
/*!
 **  Return the sill of the Nugget Effect (or TEST)
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the variables
 **
 *****************************************************************************/
double MGibbs::_getNuggetSill(Id ivar, Id jvar)
{
  CovAniso* cova = _getNugget();
  if (cova == nullptr) return (TEST);
  return (cova->getSill(ivar, jvar));
}

/****************************************************************************/
/*!
 **  Return the sill of the model (or TEST)
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the structure or of the variables
 **
 *****************************************************************************/
double MGibbs::_getCovaSill(Id ivar, Id jvar) const
{
  CovAniso* cova = _getCurrentCova();
  return (cova->getSill(ivar, jvar));
}

/****************************************************************************/
/*!
 **  Return the param of the model
 **
 *****************************************************************************/
double MGibbs::_getCovaParam(void) const
{
  return (_getCurrentCova()->getParam());
}

/****************************************************************************/
/*!
 **  Returns the number of structures in the Model (nugget excluded)
 **
 *****************************************************************************/
Id MGibbs::_getNcovaWithoutNugget(void) const

{
  Id ncova = 0;
  for (Id is = 0; is < _model->getNCov(); is++)
  {
    CovAniso* cova = _model->getCovAniso(is);
    if (cova->getType() != ECov::NUGGET) ncova++;
  }
  return (ncova);
}

/****************************************************************************/
/*!
 **  Return the number of vertices for the current Matelem
 **
 ** \param[in] icov    Rank of the target Covariance (or -1)
 **
 *****************************************************************************/
Id MGibbs::_getNvertex(Id icov)
{
  return (_getCurrentMatelem(icov).amesh->getNApices());
}

/****************************************************************************/
/*!
 **  Get the normalized range
 **
 *****************************************************************************/
double MGibbs::_getCovaRange(void) const
{
  return (_getCurrentCova()->getRangeIso());
}

/****************************************************************************/
/*!
 **  Return the total sill
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the variables
 **
 *****************************************************************************/
double MGibbs::_getSillTotal(Id ivar, Id jvar) const
{
  double total = 0.;

  CovAniso* cova = _getNugget();
  if (cova != nullptr) total += cova->getSill(ivar, jvar);

  for (Id icov = 0; icov < _getNcovaWithoutNugget(); icov++)
  {
    SPDE_CURRENT_ICOV = icov;
    cova              = _getCurrentCova();
    total += cova->getSill(ivar, jvar);
  }
  return (total);
}

/****************************************************************************/
/*!
 **  Print the Matelem characteristics (for given GRF and COV)
 **
 ** \param[in]  title   Title to be printed
 **
 *****************************************************************************/
void MGibbs::_printAll(const char* title) const
{

  /* Initializations */

  CovAniso* cova = _getCurrentCova();

  /* Print the title */

  _setTitle(true, 1, title);

  /* Global parameters */

  message("Rank of the structure = %d\n", SPDE_CURRENT_ICOV + 1);
  message("Param                 = %lf\n", _getCovaParam());
  message("Alpha                 = %lf\n", _getCovaParam() + _ndim / 2.);
  message("Total Sill            = %lf\n", _getSillTotal(0, 0));
  message("Ranges                = ");
  for (Id idim = 0; idim < _ndim; idim++)
    message("%lf ", _getCovaRange() * cova->getAnisoCoeff(idim));
  message("\n");

  /* 'H' Rotation */

  print_matrix("Anisotropy H matrix", 0, 1, _ndim, _ndim, NULL,
               Calcul.hh.getValues().data());
  message("Square root of Determinant                    = %lf\n",
          Calcul.sqdeth);
  message("Correction factor                             = %lf\n",
          Calcul.correc);

  /* Linear combination */

  Id nblin = static_cast<Id>(Calcul.blin.size());
  message("Number of terms in Linear Combination         = %d\n", nblin);
  print_matrix("Coefficients of the Linear Combination", 0, 1, 1, nblin, NULL,
               Calcul.blin.data());
}

/****************************************************************************/
/*!
 **  Compute the variance correction term
 **  Store in the SPDE_Calcul structure
 **
 *****************************************************************************/
void MGibbs::_computeCorrec(void) const
{
  double param  = _getCovaParam();
  double ndims2 = (static_cast<double>(_ndim)) / 2.;
  double gammap = exp(loggamma(param));
  double gammaa = exp(loggamma(param + ndims2));
  double g0     = pow(4. * GV_PI, ndims2);
  Calcul.correc = gammap / (g0 * gammaa);
}

/****************************************************************************/
/*!
 **  Compute the coefficients of the linear combination
 **
 ** \remarks This function stores the coefficients 'blin' in SPDE_Calcul
 **
 *****************************************************************************/
void MGibbs::_computeBLin(void) const
{
  double param  = _getCovaParam();
  double ndims2 = (static_cast<double>(_ndim)) / 2.;
  double alpha  = param + ndims2;
  Id p          = static_cast<Id>(ceil(alpha));
  Id ndimp      = p + 1;
  double lambda = alpha - floor(alpha);
  double delta  = lambda - alpha;
  double correc = Calcul.correc;

  Calcul.blin.resize(NBLIN_TERMS, 0);

  if (lambda > 0.)
  {
    /* Core allocation */

    VectorDouble v1(ndimp, 0);
    VectorDouble v2(ndimp, 0);
    MatrixSquare m(ndimp);
    MatrixSquare tp = ut_pascal(ndimp);

    for (Id idim = 0; idim < ndimp; idim++)
    {
      v1[idim] = 1. / (2. * p - idim + delta);
      for (Id jdim = 0; jdim < ndimp; jdim++)
        m.setValue(idim, jdim, 1. / (2. * p - idim - jdim + lambda));
    }
    (void)m.invert();
    m.prodMatVecInPlace(v1, v2);
    tp.prodMatVecInPlace(v2, Calcul.blin);
  }
  else
  {
    for (Id i = 0; i <= p; i++)
      Calcul.blin[i] = ut_cnp(p, i) * correc;
  }

  Calcul.blin.resize(ndimp);
}

/****************************************************************************/
/*!
 **  Compute H matrix for anisotropic case and the square root of determinant
 **  Requires the knowledge of the actual parameters of the current Covariance
 **  Fills the SPDE_Calcul structure
 **
 *****************************************************************************/
void MGibbs::_computeHH() const
{

  /* Initializations */

  CovAniso* cova = _getCurrentCova();
  MatrixSquare temp(_ndim);

  /* Processing */

  for (Id i = 0; i < _ndim; i++)
  {
    double scale = cova->getScale(i);
    if (Calcul.flag_sphere) scale /= Calcul.R;
    temp.setValue(i, i, scale * scale);
  }
  Calcul.hh.prodNormMatMatInPlace(&cova->getAnisoRotMat(), &temp, false);
}

/****************************************************************************/
/*!
 **  Initialize the contents of SPDE_Calcul structure
 **
 *****************************************************************************/
void MGibbs::_calculInitialize() const
{
  Calcul.flag_sphere = isDefaultSpaceSphere();
  Calcul.sqdeth      = 0.;
  Calcul.correc      = 0.;
  Calcul.R           = 0.;
  Calcul.hh.resize(_ndim, _ndim);
  if (Calcul.flag_sphere)
  {
    const ASpace* space = getDefaultSpaceSh().get();
    const auto* spaceSn = dynamic_cast<const SpaceSN*>(space);
    Calcul.R            = spaceSn->getRadius();
    Calcul.srot.resize(2, 0.);
  }
  Calcul.vv.resize(_ndim, 0.);
}

/****************************************************************************/
/*!
 **  Update the contents of SPDE_Calcul structure
 **
 *****************************************************************************/
void MGibbs::_calculUpdate(void)
{
  // Check that the structure has already been initiated

  if (Calcul.hh.size() <= 0)
    my_throw("You should run '_calcul_init' beforehand");

  // Calculate the 'correc' term (from 'param')
  _computeCorrec();

  // Calculate the set of 'blin' coefficients (from 'param' and 'correc')
  _computeBLin();

  // Calculate the 'HH' matrix
  _computeHH();

  // Calculate the determinant of HH
  Calcul.sqdeth = sqrt(Calcul.hh.determinant());
}

/****************************************************************************/
/*!
 **  Modify the Exponential into a Matern
 **
 ** \param[in]  cova         Covariance sructure
 **
 *****************************************************************************/
void MGibbs::_convertExponentialToMatern(CovAniso* cova) const
{
  double scale_exp, range_exp, scale_bes, range_bes;

  if (cova->getType() != ECov::EXPONENTIAL) return;

  range_exp = cova->getRangeIso();
  scale_exp = range2scale(ECov::EXPONENTIAL, range_exp, 0.);

  scale_bes = scale_exp;
  range_bes = scale2range(ECov::MATERN, scale_bes, 0.5);

  cova->setType(ECov::MATERN);
  cova->setParam(0.5);
  cova->setRangeIsotropic(range_bes);

  /* Optional printout */

  if (_verbose)
  {
    message("Convert from Exponential to Matern\n");
    message("- Exponential: Range=%lf Scale=%lf\n", range_exp, scale_exp);
    message("- Matern     : Range=%lf Scale=%lf\n", range_bes, scale_bes);
  }
}

/****************************************************************************/
/*!
 **  Check that the Model is authorized for SPDE
 **
 ** \return Error returned code
 **
 *****************************************************************************/
bool MGibbs::_checkValidAuxiliary() const
{
  if (_dbc->getNLoc(ELoc::Z) != _nvar)
  {
    messerr(
      "Model (%d) and Input Db (%d) must refer to the same number of variables",
      _nvar, _dbc->getNLoc(ELoc::Z));
    return false;
  }

  /* Check incompatibility between non-stationary and multivariate */

  if (_getNcovaWithoutNugget() > 1 || _nvar > 1 || _model->hasNugget())
    S_DECIDE.flag_several = true;

  return true;
}

/****************************************************************************/
/*!
 **  Identify a parameter among the non-stationary ones
 **
 ** \return The rank of the parameter of -1 (not found)
 **
 ** \param[in]  type0     Type of parameter (EConsElem)
 ** \param[in]  icov0     Rank of the target covariance
 ** \param[in]  ivar0     Rank of the target variable (only when type=EConsElem::SILL)
 ** \param[in]  jvar0     Rank of the target variable (only when type=EConsElem::SILL)
 **
 ** \remark The covariance are ranked from 0 for non-nugget ones
 ** \remark The nugget effect corresponds to rank (-1)
 **
 *****************************************************************************/
Id MGibbs::_identifyNostatParam(const EConsElem& type0,
                                Id icov0,
                                Id ivar0,
                                Id jvar0)
{
  DECLARE_UNUSED(type0);
  return icov0 + ivar0 + jvar0;
}

/****************************************************************************/
/*!
 **  Calculate the array of dimensions of the meshes
 **
 ** \return Pointer to the newly allocated array containing mesh dimensions
 **
 ** \param[in]  amesh    MeshEStandard structure
 **
 ** \remark The array returned by this function must be deallocated
 **
 *****************************************************************************/
VectorDouble MGibbs::_getMeshDimension(AMesh* amesh) const

{
  Id nmesh = amesh->getNMeshes();
  VectorDouble units(nmesh);

  /* Dispatch */

  if (Calcul.flag_sphere)
  {
    for (Id imesh = 0; imesh < nmesh; imesh++)
    {
      units[imesh] = GH::geodeticTriangleSurface(amesh->getCoor(imesh, 0, 0),
                                                 amesh->getCoor(imesh, 0, 1),
                                                 amesh->getCoor(imesh, 1, 0),
                                                 amesh->getCoor(imesh, 1, 1),
                                                 amesh->getCoor(imesh, 2, 0),
                                                 amesh->getCoor(imesh, 2, 1));
    }
  }
  else
  {
    MatrixSquare mat(_ndim);
    Id ncorner = amesh->getNApexPerMesh();
    for (Id imesh = 0; imesh < nmesh; imesh++)
    {
      for (Id icorn = 1; icorn < ncorner; icorn++)
        for (Id idim = 0; idim < _ndim; idim++)
          mat.setValue(idim, icorn - 1,
                       amesh->getCoor(imesh, icorn, idim) - amesh->getCoor(imesh, 0, idim));
      units[imesh] = ABS(mat.determinant()) / FACDIM[_ndim];
    }
  }
  return (units);
}

/****************************************************************************/
/*!
 **  Fill the Isill matrix linked to the covariance of the Model
 **
 ** \return Error returned code
 **
 ** \param[in] Matelem Matelem structure
 **
 ** \remark The matrix 'Isill' is dimensioned to _nvar * _nvar where
 **
 *****************************************************************************/
Id MGibbs::_fillIsill(SPDE_Matelem& Matelem) const
{
  MatrixSquare mcova(_nvar);

  /* Load the sill of the covariance */

  for (Id ivar = 0; ivar < _nvar; ivar++)
    for (Id jvar = 0; jvar < _nvar; jvar++)
      mcova.setValue(ivar, jvar, _getCovaSill(ivar, jvar));

  /* Loop on the structures to invert the sill matrices */

  mcova.invert();

  /* Optional printout */

  if (_verbose) message("Calculation of Isill\n");

  Matelem.Isill = mcova;
  return 0;
}

/****************************************************************************/
/*!
 **  Fill the Aproj matrix
 **
 ** \return Error returned code
 **
 ** \param[in] Matelem SPDE_Matelem structure
 **
 *****************************************************************************/
Id MGibbs::_fillAproj(SPDE_Matelem& Matelem) const
{

  Matelem.Aproj = dynamic_cast<MatrixSparse*>(Matelem.amesh->createProjMatrix(_dbc, -1, false));
  return (Matelem.Aproj == nullptr);
}

/****************************************************************************/
/*!
 **  Fill the Csill matrix linked to the continuous parts of the Model
 **
 ** \return Error returned code
 **
 ** \param[in] Matelem SPDE_Matelem structure
 ** \param[in] icov    Rank of the covariance
 **
 ** \remark The matrix 'Csill' is dimensioned to ncova * _nvar * (_nvar+1)/2 where
 ** \remark - ncova designates the number of continuous structures of the Model
 **
 *****************************************************************************/
Id MGibbs::_fillCsill(SPDE_Matelem& Matelem, Id icov) const
{

  /* Load the sills of continuous covariance elements */

  CholeskyDense chol(_model->getSills(icov));
  Matelem.Csill = chol.getLowerTriangle();

  /* Optional printout */

  if (_verbose) message("Calculation of Csill\n");

  return 0;
}

/****************************************************************************/
/*!
 **  Fill the Bnugget sparse matrix linked to nugget effect
 **
 ** \return Error returned code
 **
 ** \remark This function allocates 'nvs2' sparse matrices of dimension 'ndata'.
 ** \remark where nvs2 is the product _nvar * (_nvar+1) / 2
 **
 *****************************************************************************/
Id MGibbs::_fill_Bnugget() const
{
  VectorDouble local;
  VectorDouble local0;
  VectorDouble mat;
  VectorInt ind;
  Id ndata, nvs2, nvar2, size, ecr, nvr, ivar, jvar, iad;
  Id flag_nostat_sillnug;
  DECLARE_UNUSED(nvr, iad, flag_nostat_sillnug)
  std::vector<MatrixSparse*> Bnugget;

  /* Initializations */

  ndata = _dbc->getNSample(true);
  nvar2 = _nvar * _nvar;
  nvs2  = _nvar * (_nvar + 1) / 2;

  /* In the non-stationary case, identify the rank of the parameter */
  /* which corresponds to the sill of the nugget effect */

  flag_nostat_sillnug = _identifyNostatParam(EConsElem::SILL) >= 0;
  /*  if (flag_nostat_sillnug)
    {
      messerr("Non-stationarity on nugget sill values not programmed yet");
      goto label_end;
    } */

  /* Core allocation */

  size = ndata * nvs2;
  local.resize(nvar2);
  local0.resize(nvar2);
  ind.resize(ndata);
  mat.resize(size, 0);

  /* Loop on the active samples */

  ecr = 0;
  /* for (Id iech = 0; iech < dbin->getNSample(); iech++)
  {
    if (!dbin->isActive(iech)) continue;
 */
  /* Check the heterotopy for the nugget effect */

  /*   nvr = 0;
    for (ivar = 0; ivar < _nvar; ivar++)
    {
      if (FFFF(dbin->getZVariable(iech, ivar))) continue;
      ind[nvr] = ivar;
      nvr++;
    }
    if (nvr <= 0)
    {
      messerr("For sample %#d, no variable is defined", iech + 1);
      goto label_end;
    } */

  /* Dispatch */

  /*   if (nvr == _nvar && !flag_nostat_sillnug)
    { */

  /* Isotopic case: Store the sill partial matrix */

  /*    for (ivar = 0; ivar < _nvar; ivar++)
       for (jvar = 0; jvar <= ivar; jvar++)
       {
         iad = _get_rank(ivar, jvar);
         mat[iad * ndata + ecr] = LOCAL0(ivar, jvar);
       }
   }
   else
   { */

  /* Constitute the sill matrix for the nugget effect */

  /*    for (Id ivr = 0; ivr < nvr; ivr++)
       for (Id jvr = 0; jvr < nvr; jvr++)
         LOCAL(ivr,jvr) = _get_nugget_sill(ind[ivr], ind[jvr]);
*/
  /* Invert the sill partial matrix */

  /*    if (matrix_invert(local, nvr, -1))
     {
       messerr("Problem when inverting Nugget matrix of sill at sample #%d",
               iech + 1);
       goto label_end;
     } */

  /* Store the sill partial matrix */

  /*   for (Id ivr = 0; ivr < nvr; ivr++)
      for (Id jvr = 0; jvr <= ivr; jvr++)
      {
        ivar = ind[ivr];
        jvar = ind[jvr];
        iad = _get_rank(ivar, jvar);
        mat[iad * ndata + ecr] = LOCAL(ivr, jvr);
      }
  }
  ecr++;
}
*/
  /* Define the sparse matrices */

  ecr = 0;
  for (ivar = 0; ivar < _nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ecr++)
    {
      VectorDouble diag = VH::initVDouble(&mat[ecr * ndata], ndata);
      Bnugget.push_back(MatrixSparse::diagVec(diag));
    }

  /* Optional printout */

  if (_verbose) message("Calculation of Bnugget (%d sparse matrices)\n", nvs2);

  S_ENV.Bnugget = Bnugget;
  S_ENV.ndata   = ndata;
  return 0;
}

/****************************************************************************/
/*!
 **  Return the list of (target+data) indices for a given mesh
 **
 ** \return An array of vertex identification (Dimension: nvertex) or NULL
 **
 ** \param[in]  amesh       AMesh structure
 **
 ** \remarks The array ranks is filled as follows:
 ** \remarks - Its contents follows the mesh numbering
 ** \remarks - If positive, its value provides the rank of the data
 ** \remarks - If negative, its absolute value provides the rank of the target
 ** \remarks - If zero, theses are Steiner points
 ** \remarks Warning: Ranks are counted from 1
 **
 ** \remarks The returned array must be freed by the calling function
 ** \remarks Dimension: nvertex
 **
 *****************************************************************************/
VectorInt MGibbs::_getVertexRanks(AMesh* amesh) const
{
  Id nvertex = amesh->getNApices();
  Id n_in    = (_dbc != nullptr) ? _dbc->getNSample(true) : 0;
  Id n_out   = _dbout->getNSample(true);
  if (nvertex < (n_in + n_out))
    messageAbort("Nvertex(%d) must be larger than n_in(%d) + n_out(%d)",
                 nvertex, n_in, n_out);

  /* Core allocation */

  VectorInt ranks(nvertex, 0);

  /* Identify the vertices */

  Id ecr = 0;
  if (_dbc != nullptr)
    for (Id i = 0; i < _dbc->getNSample(); i++)
    {
      if (!_dbc->isActive(i)) continue;
      ranks[ecr++] = (i + 1);
    }

  for (Id i = 0; i < _dbout->getNSampleActive(); i++)
  {
    if (!_dbout->isActive(i)) continue;
    ranks[ecr++] = -(i + 1);
  }
  return (ranks);
}

/****************************************************************************/
/*!
 **  Fill some matrices for Kriging in the case of a model without nugget effect
 **  Constitute the Bhetero sparse matrices
 **
 ** \return Error returned code
 **
 *****************************************************************************/
Id MGibbs::_fillBhetero() const

{
  VectorInt ranks;
  VectorInt ndata1;
  VectorInt ntarget1;
  Id ndata, ecrT, nvertex, flag_add, iech;
  double value;
  std::vector<MatrixSparse*> BheteroD;
  std::vector<MatrixSparse*> BheteroT;
  AMesh* amesh;

  /* Initializations */

  ndata              = _dbc->getNSample(true);
  SPDE_Matelem& Mat1 = _getCurrentMatelem(0);
  amesh              = Mat1.amesh;
  nvertex            = amesh->getNApices();

  /* Core allocation */

  ranks = _getVertexRanks(amesh);

  /* Define the sparse matrices */

  ndata1.resize(_nvar, 0.);
  ntarget1.resize(_nvar, 0);
  BheteroD.resize(_nvar);
  for (Id ivar = 0; ivar < _nvar; ivar++)
    BheteroD[ivar] = nullptr;
  BheteroT.resize(_nvar);
  for (Id ivar = 0; ivar < _nvar; ivar++)
    BheteroT[ivar] = nullptr;

  /**************************************************************/
  /* Creating the sparse matrix for handling heterotopy on data */
  /**************************************************************/
  /* This matrix dimension is [NDmax , Nvertex]                 */
  /* where NDmax is the number of active samples of Dbin        */
  /* regardless of their contents (heterotopy)                  */
  /* A line (dbin sample) contains the barycenter weights       */
  /* assigned to the vertices of the mesh to which it belongs   */
  /* (if the variable is defined for this sample); 0 otherwise  */

  for (Id ivar = 0; ivar < _nvar; ivar++)
  {
    NF_Triplet Btriplet;
    for (Id i = 0; i < nvertex; i++)
    {
      if (ranks[i] <= 0) continue; // Target or Steiner
      ndata1[ivar]++;
      iech  = ranks[i] - 1;
      value = (FFFF(_dbc->getZVariable(iech, ivar))) ? 0. : 1.;
      Btriplet.add(iech, i, value);
    }
    // Add a fictitious sample (zero value) as a dimension constraint
    Btriplet.force(ndata1[ivar], nvertex);
    BheteroD[ivar] = MatrixSparse::createFromTriplet(Btriplet);
  }

  /* Optional printout */

  if (_verbose)
    message("Calculation of Bhetero for Data (%d sparse matrices)\n", _nvar);

  /**********************************************************************/
  /* Creating the sparse matrix for handling heterotopy on target       */
  /**********************************************************************/
  /* Per variable, this matrix should mask off all mesh vertices which  */
  /* correspond to a data sample which is defined for this variable but */
  /* not defined for all variables                                      */

  for (Id ivar = 0; ivar < _nvar; ivar++)
  {
    NF_Triplet Btriplet;
    ecrT = 0;
    for (Id i = 0; i < nvertex; i++)
    {
      flag_add = 0;
      if (ranks[i] <= 0)
      {
        // This is a target (from dbout) or a Steiner
        flag_add = 1;
      }
      else
      {
        // This is a data
        iech = ranks[i] - 1;

        // Could the data be considered as a target (heterotopic case)
        if (FFFF(_dbc->getZVariable(iech, ivar)))
        {
          // The sample is not defined for the current variable: it is a target
          flag_add = 1;
        }
      }
      if (!flag_add) continue;
      ntarget1[ivar]++;
      Btriplet.add(ecrT, i, 1.);
      ecrT++;
    }

    // Add a fictitious sample (zero value) as a dimension constraint
    Btriplet.force(ntarget1[ivar], nvertex);
    BheteroT[ivar] = MatrixSparse::createFromTriplet(Btriplet);
  }

  /* Optional printout */

  if (_verbose)
    message("Calculation of Bhetero for Target (%d sparse matrices)\n", _nvar);

  /* Set the error return code */

  S_ENV.BheteroD = BheteroD;
  S_ENV.BheteroT = BheteroT;
  S_ENV.ndata    = ndata;
  S_ENV.ndata1   = ndata1;
  S_ENV.ntarget1 = ntarget1;
  return 0;
}

/****************************************************************************/
/*!
 **  Get the 3-D coordinates of the center of a triangle on the sphere
 **
 ** \param[in]  amesh    MeshEStandard structure
 ** \param[in]  ncorner  Number of vertices per element
 ** \param[in]  imesh    Rank of the current mesh
 **
 ** \param[out] center   Coordinates of the center point (Dimension: 3)
 ** \param[out] xyz      Coordinate of the point (Dimension: 3x3)
 **
 *****************************************************************************/
void MGibbs::_calculateTriangleCenter(AMesh* amesh,
                                      Id ncorner,
                                      Id imesh,
                                      double center[3],
                                      double xyz[3][3])
{
  double ratio;

  for (Id i = 0; i < 3; i++)
    center[i] = 0.;
  for (Id icorn = 0; icorn < ncorner; icorn++)
  {
    GH::convertSph2Cart(amesh->getCoor(imesh, icorn, 0),
                        amesh->getCoor(imesh, icorn, 1),
                        &xyz[icorn][0], &xyz[icorn][1], &xyz[icorn][2]);
    for (Id i = 0; i < 3; i++)
      center[i] += xyz[icorn][i];
  }

  ratio = 0.;
  for (Id i = 0; i < 3; i++)
  {
    center[i] /= 3.;
    ratio += center[i] * center[i];
  }
  ratio = 1. / sqrt(ratio);
  for (Id i = 0; i < 3; i++)
    center[i] *= ratio;
}

/****************************************************************************/
/*!
 **  Project a point on the tangent plane
 **
 ** \param[in]  center    Coordinates of the reference point (Dimension: 3)
 ** \param[in]  axes      Coordinates of the endpoints (Dimension: 2 * 3)
 ** \param[in]  xyz       Coordinates of the target point (Dimension: 3)
 **
 ** \param[out] coeff     Coordinate of point in the local system (Dimension: 2)
 **
 *****************************************************************************/
void MGibbs::_calculateProjectPlane(double center[3],
                                    double axes[2][3],
                                    double xyz[3],
                                    double coeff[2])
{
  double v[3];

  /* Projection */

  for (Id j = 0; j < 2; j++)
  {
    coeff[j] = 0.;
    for (Id i = 0; i < 3; i++)
      coeff[j] += (axes[j][i] - center[i]) * (xyz[i] - center[i]);
  }

  /* Projected vector */

  for (Id i = 0; i < 3; i++)
  {
    v[i] = 0.;
    for (Id j = 0; j < 2; j++)
      v[i] += coeff[j] * (axes[j][i] - center[i]);
  }

  /* Returned coordinates */

  VH::addInPlace(center, v, xyz, 3);
}

/****************************************************************************/
/*!
 **  Get the coordinates of the axis endpoints in the tangent plane
 **
 ** \param[in]  center  Coordinates of the reference point (Dimension: 3)
 ** \param[in]  srot    Rotation angles on sphere (Dimension: 2)
 **
 ** \param[out] axes    Coordinates of the endpoints (Dimension: 2 * 3)
 **
 *****************************************************************************/
void MGibbs::_calculateTangent(double center[3],
                               const double srot[2],
                               double axes[2][3])
{
  double sinphi, cosphi, sintet, costet, theta, phi, v[3], w[3];

  // Center gives the vector joining the origin to the center of triangle
  phi    = srot[1] * GV_PI / 180.;
  theta  = srot[0] * GV_PI / 180.;
  sinphi = sin(phi);
  cosphi = cos(phi);
  sintet = sin(theta);
  costet = cos(theta);
  // W is the Pole vector
  w[0] = sinphi * costet;
  w[1] = sinphi * sintet;
  w[2] = cosphi;
  // V = Center ^ w: first axis
  VH::crossProduct3DInPlace(center, w, v);
  VH::normalize(v, 3);
  // W = Center ^ V: second axis
  VH::crossProduct3DInPlace(center, v, w);
  VH::normalize(w, 3);
  // Get the end points from Unit vectors
  VH::addInPlace(center, v, axes[0], 3);
  VH::addInPlace(center, w, axes[1], 3);
}

/****************************************************************************/
/*!
 **  Fill the sparse matrix S linked to mesh vertices
 **
 ** \return G sparse matrix
 **
 ** \param[in]  amesh     MeshEStandard_Mesh structure
 ** \param[in]  units     Array containing the mesh dimensions
 **
 *****************************************************************************/
MatrixSparse* MGibbs::_fillS(AMesh* amesh, const VectorDouble& units)
{
  double xyz[3][3], center[3], axes[2][3], coeff[3][2], vald;
  Id errcod;
  long ip1, ip2;
  MatrixSparse* G = nullptr;
  std::map<std::pair<Id, Id>, double> tab;
  std::pair<std::map<std::pair<Id, Id>, double>::iterator, bool> ret;
  std::map<std::pair<Id, Id>, double>::iterator it;

  /* Initializations */

  Id error   = 1;
  Id ncorner = amesh->getNApexPerMesh();
  NF_Triplet Gtriplet;
  bool flag_sphere = isDefaultSpaceSphere();
  bool flag_nostat = false;
  if (!flag_nostat) _calculUpdate();
  MatrixSquare mat(ncorner);
  MatrixSquare matu(ncorner);
  MatrixDense matw(_ndim, ncorner);
  MatrixDense mat1(ncorner, _ndim);
  VectorDouble matv(3);

  /* Loop on the meshes */

  for (Id imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {

    /* Get parameters in the non-stationary case */
    // Not coded in this version
    // if (flag_nostat) _calcul_update_nostat(amesh, imesh);

    // Processing on the Sphere

    if (flag_sphere)
    {

      // Case of the calculations on the Sphere

      _calculateTriangleCenter(amesh, ncorner, imesh, center, xyz);
      if (ncorner < 0 || ncorner > 3)
      {
        messerr("Error in _triangle_center: wrong number or corners: %d", ncorner);
        goto label_end;
      }

      /* Look for the tangent plane and its axes */

      _calculateTangent(center, Calcul.srot.data(), axes);

      /* Project corner points on the Tangent plane */

      for (Id icorn = 0; icorn < ncorner; icorn++)
        _calculateProjectPlane(center, axes, xyz[icorn], coeff[icorn]);

      for (Id icorn = 0; icorn < ncorner; icorn++)
      {
        for (Id idim = 0; idim < _ndim; idim++)
          matu.setValue(icorn, idim, coeff[icorn][idim]);
        matu.setValue(icorn, ncorner - 1, 1.);
      }
    }
    else
    {

      // Case of Euclidean geometry

      for (Id icorn = 0; icorn < ncorner; icorn++)
      {
        for (Id idim = 0; idim < _ndim; idim++)
          matu.setValue(icorn, idim, amesh->getCoor(imesh, icorn, idim));
        matu.setValue(icorn, ncorner - 1, 1.);
      }
    }

    /* Invert the matrix 'matu'*/

    errcod = matu.invert();
    if (errcod)
    {
      messerr("Error in Mesh #%3d - Its volume is zero", imesh + 1);
      for (Id icorn = 0; icorn < ncorner; icorn++)
      {
        message("Sample #%4d - Coordinates (", amesh->getApex(imesh, icorn));
        for (Id idim = 0; idim < _ndim; idim++)
          message(" %lf", amesh->getCoor(imesh, icorn, idim));
        message(")\n");
      }
      print_matrix("MATU", 0, 1, ncorner, ncorner, NULL, matu.getValues().data());
    }
    else
    {
      for (Id icorn = 0; icorn < ncorner; icorn++)
        for (Id idim = 0; idim < _ndim; idim++)
          matw.setValue(idim, icorn, matu.getValue(idim, icorn));

      mat1.prodMatMatInPlace(&matw, &Calcul.hh, true, false);
      if (flag_nostat)
        matw.prodMatVecInPlace(Calcul.vv, matv, true);
      mat.prodMatMatInPlace(&mat1, &matw);

      for (Id j0 = 0; j0 < ncorner; j0++)
        for (Id j1 = 0; j1 < ncorner; j1++)
        {
          ip1 = amesh->getApex(imesh, j0);
          ip2 = amesh->getApex(imesh, j1);
          std::pair<Id, Id> key(ip1, ip2);
          vald = units[imesh] * mat.getValue(j0, j1);
          if (flag_nostat) vald += matv[j1] * units[imesh];
          ret = tab.insert(std::pair<std::pair<Id, Id>, double>(key, vald));
          if (!ret.second) ret.first->second += vald;
        }
    }
  }

  it = tab.begin();
  while (it != tab.end())
  {
    ip1 = it->first.first;
    ip2 = it->first.second;
    Gtriplet.add(ip1, ip2, it->second);
    it++;
  }

  /* Optional printout */

  G = MatrixSparse::createFromTriplet(Gtriplet);

  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    delete G;
    G = nullptr;
  }
  return (G);
}

/****************************************************************************/
/*!
 **  Fill the vector TildeC (Dimension: nvertex)
 **
 ** \return Error return code
 **
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  units     Array containing the element units
 **
 *****************************************************************************/
VectorDouble MGibbs::_fillTildeC(AMesh* amesh, const VectorDouble& units)
{
  VectorDouble tildec, cumunit;
  Id nvertex = amesh->getNApices();
  Id ncorner = amesh->getNApexPerMesh();
  cumunit.resize(nvertex, 0);

  /* Loop on the meshes */

  for (Id imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {

    /* Loop on the vertices */

    for (Id icorn = 0; icorn < ncorner; icorn++)
    {
      Id ip = amesh->getApex(imesh, icorn);
      cumunit[ip] += units[imesh];
    }
  }

  /* Scale */

  auto factor = static_cast<double>(ncorner);
  for (Id ip = 0; ip < nvertex; ip++)
  {
    double value = cumunit[ip] / factor;
    if (ABS(value) <= 0.)
    {
      messerr("Meshing unit (%d) has a zero volume", ip + 1);
      return VectorDouble();
    }
    tildec.push_back(value);
  }
  return tildec;
}

/****************************************************************************/
/*!
 **  Fill the vector for sill correction factors
 **  Works for both stationary and non-stationary cases
 **
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  TildeC    Vector TildeC
 **
 *****************************************************************************/
VectorDouble MGibbs::_fillLambda(AMesh* amesh, const VectorDouble& TildeC)
{
  VectorDouble Lambda;
  Id nvertex  = amesh->getNApices();
  double sill = _getCovaSill(0, 0);

  /* Fill the array */

  double sqdeth = Calcul.sqdeth;
  for (Id ip = 0; ip < nvertex; ip++)
    Lambda.push_back(sqrt((TildeC[ip]) / (sqdeth * sill)));

  return (Lambda);
}

/****************************************************************************/
/*!
 **  Extract the sparse matrix from the Q matrix (case of nugget effect)
 **
 ** \return The extracted sparse matrix or NULL
 **
 ** \param[in]  row_var    Rank of the variable for the row
 ** \param[in]  col_var    Rank of the variable for the column
 **
 ** \param[out] nrows      Number of rows
 ** \param[out] ncols      Number of columns
 **
 ** \remarks Extracts a part of Bnugget matrix for:
 ** \remarks - a given pair of variables
 ** \remarks - for Data-Data operators
 **
 *****************************************************************************/
MatrixSparse* MGibbs::_extractQ1Nugget(Id row_var,
                                       Id col_var,
                                       Id* nrows,
                                       Id* ncols)
{
  MatrixSparse* B0 = S_ENV.Bnugget[_getRank(row_var, col_var)]->clone();
  if (B0 != nullptr)
  {
    *nrows = B0->getNRows();
    *ncols = B0->getNCols();
  }
  return (B0);
}

/****************************************************************************/
/*!
 **  Extract the sparse matrix from the Q matrix (coninuous structure)
 **
 ** \return The extracted sparse matrix or NULL
 **
 ** \param[in]  row_var    Rank of the variable for the row
 ** \param[in]  col_var    Rank of the variable for the column
 ** \param[in]  row_oper   Operator type for row (1:Data or 2:Target)
 ** \param[in]  col_oper   Operator type for column (1:Data or 2:Target)
 **
 ** \param[out] nrows      Number of rows
 ** \param[out] ncols      Number of columns
 **
 ** \remarks Extracts a part of Q matrix (for the first structure) for:
 ** \remarks - a given pair of variables
 ** \remarks - a given pair of operators (Data or Target)
 ** \remarks The returned matrix is multipled by the inverse of the Sill
 **
 *****************************************************************************/
MatrixSparse* MGibbs::_extracQ1Hetero(Id row_var,
                                      Id col_var,
                                      Id row_oper,
                                      Id col_oper,
                                      Id* nrows,
                                      Id* ncols)
{
  Id error;
  MatrixSparse *Q, *Brow, *Bcol, *B1, *Bt, *Qn;

  /* Initializations */

  error = 1;
  Q = Brow = Bcol = B1 = Bt = Qn = nullptr;
  SPDE_Matelem& Matelem1         = _getCurrentMatelem(0);

  /* Identify the operating matrices */

  Brow = (row_oper == 1) ? S_ENV.BheteroD[row_var] : S_ENV.BheteroT[row_var];
  if (Brow == nullptr) goto label_end;
  Bcol = (col_oper == 1) ? S_ENV.BheteroD[col_var] : S_ENV.BheteroT[col_var];
  if (Bcol == nullptr) goto label_end;
  Bt = Bcol->transpose();
  if (Bt == nullptr) goto label_end;
  B1 = MatrixFactory::prodMatMat<MatrixSparse>(Brow, Matelem1.QC->Q);
  if (B1 == nullptr) goto label_end;
  Qn = MatrixFactory::prodMatMat<MatrixSparse>(B1, Bt);
  if (Qn == nullptr) goto label_end;

  /* Multiply by the corresponding sill */

  Q = MatrixSparse::addMatMat(Qn, Qn, _getIsill(0, row_var, col_var), 0.);
  if (Q == nullptr) goto label_end;

  /* Set the error return code */

  error  = 0;
  *nrows = (row_oper == 1) ? S_ENV.ndata1[row_var] : S_ENV.ntarget1[row_var];
  *ncols = (col_oper == 1) ? S_ENV.ndata1[col_var] : S_ENV.ntarget1[col_var];

label_end:
  delete B1;
  delete Bt;
  delete Qn;
  if (error) delete Q;
  return (Q);
}

/****************************************************************************/
/*!
 **  Construct the sparse matrix QCov (used in multistructure - multivariable)
 **
 ** \return Error return code
 **
 ** \param[in]  Matelem     SPDE_Matelem structure
 **
 ** \remarks This function requires the Q matrices to be established already,
 ** \remarks as well as the Aproj matrices.
 ** \remarks In case of presence of nugget effect, we also need 'Bnugget'
 ** \remarks Otherwise, we need 'BheteroD' and 'BheteroT'
 **
 *****************************************************************************/
Id MGibbs::_buildQCov(SPDE_Matelem& Matelem)

{
  Id error, icov0, nrows, ncols;
  MatrixSparse *B0, *Bi;
  std::vector<QChol*> QCov;

  /* Initializations */

  if (!S_DECIDE.flag_several) return (0);
  error = 1;
  Bi = B0 = nullptr;
  icov0   = SPDE_CURRENT_ICOV;

  /* Core allocation */

  QCov.resize(_nvar);
  for (Id ivar = 0; ivar < _nvar; ivar++)
    QCov[ivar] = _manageQchol(1, NULL);

  /* Dispatch */

  if (_model->hasNugget())
  {

    /****************************************/
    /* Case when a nugget effect is present */
    /****************************************/

    if (Matelem.Aproj == NULL || S_ENV.Bnugget.empty()) return (1);

    for (Id ivar = 0; ivar < _nvar; ivar++)
    {
      // Sill(icov)_ii * Q(icov) + A^t(icov) * E_ii * A(icov)
      B0 = _extractQ1Nugget(ivar, ivar, &nrows, &ncols);
      if (B0 == nullptr) goto label_end;
      Bi = prodNormMatMat(B0, Matelem.Aproj, true);
      if (Bi == nullptr) goto label_end;
      QCov[ivar]->Q = MatrixSparse::addMatMat(Matelem.QC->Q, Bi, _getIsill(icov0, ivar, ivar));
      if (QCov[ivar]->Q == nullptr) goto label_end;
      delete Bi;
      delete B0;
    }
  }
  else
  {

    /***************************************/
    /* Case when there is no nugget effect */
    /***************************************/

    if (Matelem.Aproj == NULL || S_ENV.BheteroD.empty() || S_ENV.BheteroT.empty())
      return (1);

    for (Id ivar = 0; ivar < _nvar; ivar++)
    {
      if (icov0 == 0)
      {
        // Q1_tt_ii
        QCov[ivar]->Q = _extracQ1Hetero(ivar, ivar, 2, 2, &nrows, &ncols);
        if (QCov[ivar]->Q == nullptr) goto label_end;
      }
      else
      {
        // Sill(icov)_ii * Q(icov) + A^t(icov) * Q1_dd_ii * A(icov)
        B0 = _extracQ1Hetero(ivar, ivar, 1, 1, &nrows, &ncols);
        if (B0 == nullptr) goto label_end;
        Bi = prodNormMatMat(B0, Matelem.Aproj, true);
        if (Bi == nullptr) goto label_end;
        QCov[ivar]->Q = MatrixSparse::addMatMat(Matelem.QC->Q, Bi, _getIsill(icov0, ivar, ivar));
        if (QCov[ivar]->Q == nullptr) goto label_end;
        delete Bi;
        delete B0;
      }
    }
  }
  Matelem.QCov = QCov;

  /* Optional printout */

  if (_verbose) message("Building QCov (%d sparse matrices)\n", _nvar);

  /* Set the error return code */

  error = 0;

label_end:
  delete B0;
  delete Bi;
  if (error)
  {
    if (!QCov.empty())
    {
      for (Id ivar = 0; ivar < _nvar; ivar++)
        QCov[ivar] = _manageQchol(-1, QCov[ivar]);
    }
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Construct the final sparse matrix Q from the Model
 **
 ** \return Error return code
 **
 ** \param[in]  Matelem    SPDE_Matelem structure
 **
 *****************************************************************************/
Id MGibbs::_buildQ(SPDE_Matelem& Matelem) const

{
  Id nblin   = static_cast<Id>(Calcul.blin.size());
  Id nvertex = Matelem.S->getNCols();
  if (nvertex <= 0)
  {
    messerr("You must define a valid Meshing beforehand");
    return 1;
  }
  if (nblin <= 0)
  {
    messerr("You must have a set of already available 'blin' coefficients");
    messerr("These coefficients come from the decomposition in series for Q");
    messerr("This decomposition is available only if 'alpha' is an integer");
    messerr("where: alpha = param + ndim/2");
    return 1;
  }

  //  Build Q within QC

  Matelem.QC       = _manageQchol(1, NULL);
  MatrixSparse* Q  = MatrixSparse::diagConstant(nvertex, Calcul.blin[0]);
  MatrixSparse* Bi = Matelem.S->clone();

  /* Loop on the different terms */

  for (Id iterm = 1; iterm < nblin; iterm++)
  {
    Q->addMat(*Bi, 1., Calcul.blin[iterm]);
    if (iterm < nblin - 1)
      Bi->prodMat(Matelem.S);
  }
  delete Bi;

  /* Final scaling */

  Q->prodNormDiagVecInPlace(Matelem.Lambda, 1);
  Matelem.QC->Q = Q;

  /* Optional printout */

  if (_verbose) message("Building Global Q matrix\n");

  return 0;
}

/****************************************************************************/
/*!
 **  Build all matrices needed for establishing the Q sparse matrix
 **
 ** \return Error return code
 **
 ** \param[in]  Matelem    Matelem structure
 **
 ** \remarks Contents of SP_MAT (sparse matrices or vectors) is allocated here
 ** \remarks It must be freed by the calling functions
 **
 *****************************************************************************/
Id MGibbs::_buildMatrices(SPDE_Matelem& Matelem)
{
  AMesh* amesh = Matelem.amesh;

  /* Calculate the units of the meshes */

  VectorDouble units = _getMeshDimension(amesh);

  /* Fill S sparse matrix */

  Matelem.S = _fillS(amesh, units);
  if (Matelem.S == nullptr) return 1;
  if (_verbose) message("Filling S Sparse Matrix performed successfully\n");

  /* Fill the TildeC vector */

  VectorDouble tildeC = _fillTildeC(amesh, units);
  if (_verbose) message("Filling TildeC Sparse Matrix performed successfully\n");

  /* Construct the matrix for the sill correction array */

  Matelem.Lambda = _fillLambda(amesh, tildeC);
  if (_verbose) message("Filling Lambda Sparse Matrix performed successfully\n");

  /* Build the sparse matrix B */

  Matelem.S->prodNormDiagVecInPlace(tildeC, 2);

  /* Build the sparse matrix Q */

  if (_buildQ(Matelem)) return 1;

  return 0;
}

/****************************************************************************/
/*!
 **  Internal function used for the Chebychev approximation
 **
 ** \return  Returned value
 **
 ** \param[in]  x           Input value
 ** \param[in]  power       Parameter used in the Chebychev approximation
 ** \param[in]  blin        Array of coefficients for Linear combination
 **
 *****************************************************************************/
double MGibbs::_chebychevFunction(double x,
                                  double power,
                                  const VectorDouble& blin)
{
  double value, total;

  value = 1.;
  total = blin[0];
  for (Id i = 1, nblin = static_cast<Id>(blin.size()); i < nblin; i++)
  {
    value *= x;
    total += blin[i] * value;
  }
  if (power == 0.) return (log(total));
  return (pow(total, power));
}

/****************************************************************************/
/*!
 **  Evaluate the number of coefficients necessary to evaluate a function
 **  (at a sample location) at a given approximation
 **
 ** \return Error return code
 **
 ** \param[in]  cheb_elem  Cheb_Elem structure to be filled
 ** \param[in]  blin       Array of coefficients for Linear combination
 **
 *****************************************************************************/
Id MGibbs::_chebychevCoefficients(Cheb_Elem* cheb_elem,
                                  const VectorDouble& blin) const

{
  double a = cheb_elem->a;
  double b = cheb_elem->b;
  Id ndisc = cheb_elem->ndisc;

  /* Calculate the polynomials */

  cheb_elem->coeffs.resize(cheb_elem->ncmax);

  /* Evaluate the coefficients of the Chebychev approximation */

  if (ut_chebychev_coeffs(_chebychevFunction, cheb_elem, blin)) return 1;

  /* Loop on some discretized samples of the interval */

  Id number = 0;
  for (Id idisc = 1; idisc < ndisc; idisc++)
  {
    double value = a + (b - a) * idisc / ndisc;
    Id numloc    = ut_chebychev_count(_chebychevFunction, cheb_elem, value, blin);
    if (numloc > number) number = numloc;
  }

  /* Optional printout */

  if (_verbose)
  {
    message("Chebychev Polynomial Approximation:\n");
    message("- Power = %lf\n", cheb_elem->power);
    message("- Performed using %d terms\n", number);
    message("- between %lf and %lf (Nb. discretization steps=%d)\n", a, b, ndisc);
    message("- with a tolerance of %lg\n", cheb_elem->tol);
  }

  /* Core Reallocation */

  cheb_elem->coeffs.resize(number);
  cheb_elem->ncoeffs = number;
  return 0;
}

/****************************************************************************/
/*!
 **  Manage Cheb_Elem structure
 **
 ** \return  Error return code
 **
 ** \param[in]  Matelem    Matelem structure
 ** \param[in]  power      Parameter passed to Chebychev function
 ** \param[in]  blin       Array of coefficients for Linear combinaison
 ** \param[in]  S          Shift operator
 **
 ** \remarks Arguments 'power', 'nblin', 'blin' and 'B' are used if mode=1
 ** \remarks Argument 'cheb_old' is used if mode=-1
 **
 *****************************************************************************/
Id MGibbs::_manageCheb(SPDE_Matelem& Matelem,
                       double power,
                       const VectorDouble& blin,
                       MatrixSparse* S)
{
  auto* cheb_elem = new Cheb_Elem();

  // Allocation

  Id ncmax   = static_cast<Id>(get_keypone("Number_Polynomials_Chebychev", 10001.));
  Id ndisc   = static_cast<Id>(get_keypone("Number_Discretization_Chebychev", 100.));
  double tol = get_keypone("Chebychev_Tolerance", 5.e-3);

  /* Calculate key values */

  double a  = 0.;
  double b  = S->L1Norm();
  double v1 = 2. / (b - a);
  double v2 = -(b + a) / (b - a);

  /* Store the values */

  cheb_elem->a       = a;
  cheb_elem->b       = b;
  cheb_elem->v1      = v1;
  cheb_elem->v2      = v2;
  cheb_elem->power   = power;
  cheb_elem->ncmax   = ncmax;
  cheb_elem->ndisc   = ndisc;
  cheb_elem->tol     = tol;
  cheb_elem->ncoeffs = 0;

  /* Get the optimal count of Chebychev coefficients */

  if (_chebychevCoefficients(cheb_elem, blin))
  {
    delete cheb_elem;
    return 1;
  }

  Matelem.s_cheb = cheb_elem;
  return 0;
}

/****************************************************************************/
/*!
 **  Initialize one SP_Mat structure
 **
 ** \param[in] mode    Type of the action
 **                    1 for allocation;
 **                    0 for partial deallocation (of current Matelem)
 **                   -1 for deallocation
 **
 ** \remarks This function is called when the current IGRF has been chosen
 **
 *****************************************************************************/
void MGibbs::_manageMatelem(Id mode)

{
  auto ncova = _getNcovaWithoutNugget();

  /* Dispatch */

  switch (mode)
  {
    case 1: // Allocation
      S_ENV.Matelems.resize(ncova);

      for (Id is = 0; is < ncova; is++)
      {
        SPDE_Matelem& Matelem = S_ENV.Matelems[is];
        Matelem.S             = nullptr;
        Matelem.Aproj         = nullptr;
        Matelem.QC            = nullptr;
        Matelem.QCtt          = nullptr;
        Matelem.QCtd          = nullptr;
        Matelem.s_cheb        = nullptr;
        Matelem.amesh         = nullptr;
      }
      break;

    case -1: // Deallocation
      for (Id icov = 0; icov < ncova; icov++)
      {
        SPDE_Matelem& Matelem = _getCurrentMatelem(icov);
        delete Matelem.S;
        delete Matelem.Aproj;
        Matelem.QC = _manageQchol(-1, Matelem.QC);
        if (!Matelem.QCov.empty())
        {
          for (Id ivar = 0; ivar < _nvar; ivar++)
            Matelem.QCov[ivar] = _manageQchol(-1, Matelem.QCov[ivar]);
        }
        Matelem.QCtt = _manageQchol(-1, Matelem.QCtt);
        Matelem.QCtd = _manageQchol(-1, Matelem.QCtd);
        Matelem.Isill.clear();
        Matelem.Csill.clear();
        delete Matelem.s_cheb;
        delete Matelem.amesh;
        Matelem.amesh = nullptr;
      }
      break;
  }
}

/****************************************************************************/
/*!
 **  Load the meshes
 **
 ** \return  Error return code
 **
 ** \param[in]  Matelem    Matelem structure
 **
 ** \remarks The option 'flag_force' forces to use the regular meshing rather
 ** \remarks than the Turbo one
 **
 *****************************************************************************/
Id MGibbs::_createMeshes(SPDE_Matelem& Matelem)
{
  // Processing

  if (isDefaultSpaceSphere())
  {

    /* Particular case of data on the sphere */

    messerr("This is not possible in Standard Meshing technique");
    messerr("Use MeshEStandardExt meshing technique instead");
    return 1;
  }

  /* Standard case */
  auto* dbgrid = dynamic_cast<DbGrid*>(_dbout);
  if (dbgrid == nullptr)
  {
    messerr("This type of Meshing is not available in Standard");
    messerr("Use MeshEStandardExt meshing technique instead");
    return 1;
  }

  /* Regular meshing */

  if (_verbose) message("Using Turbo Meshing\n");
  MeshETurbo* mesh = MeshETurbo::createFromGrid(dbgrid, false, _verbose);
  mesh->setPolarized(false);
  Matelem.amesh = mesh;
  return 0;
}

/****************************************************************************/
/*!
 **  Preparation using SPDE (for all GRF and COV)
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id MGibbs::_prepar()
{
  _calculInitialize();

  /* Title (optional) */

  if (_verbose) _setTitle(false, 1, "Preparing the Environment");

  /* Prepare the array of inverse of nugget sill matrices */

  if (S_DECIDE.flag_several && _model->hasNugget())
  {
    if (_fill_Bnugget()) return 1;
  }

  /* Loop on the covariances */

  for (Id icov = 0, ncova = _getNcovaWithoutNugget(); icov < ncova; icov++)
  {
    SPDE_CURRENT_ICOV     = icov;
    SPDE_Matelem& Matelem = _getCurrentMatelem(icov);

    /* Title (optional) */

    if (_verbose) _setTitle(true, 1, "Preparing the Process");

    /* Load the AMesh structure */

    if (_createMeshes(Matelem)) return 1;

    /* Prepare the array of sparse matrices (without nugget effect) */

    if (S_DECIDE.flag_several && !_model->hasNugget())
    {
      if (_fillBhetero()) return 1;
    }

    /* Prepare the projection matrix */

    if (_fillAproj(Matelem)) return 1;

    /* Prepare the kriging environment per structure */

    if (S_DECIDE.flag_several)
    {
      if (_fillIsill(Matelem)) return 1;
    }

    /* Prepare the simulation environment per structure */

    if (S_DECIDE.flag_case == CASE_SIMULATE)
    {
      if (_fillCsill(Matelem, icov)) return 1;
    }

    /* Build all relevant matrices */

    if (_buildMatrices(Matelem)) return 1;

    /* Build additional matrices */

    if (_buildQCov(Matelem)) return 1;

    /* Partially free the SP_Mat structure */

    _manageMatelem(0);

    /* Building simulation or Kriging environment */

    Matelem.QCtt = _extractQCfromQ("f_f", _getCurrentMatelem(-1).QC,
                                   VT_FREE, VT_FREE);
    if (Matelem.QCtt == nullptr) return 1;

    /* Prepare the Chebychev simulation environment */

    if (_manageCheb(Matelem, -0.5, Calcul.blin, Matelem.S)) return 1;

    /* Verbose output (optional) */

    if (DEBUG && _verbose) _printMatelem(icov);
  }

  SPDE_CURRENT_ICOV = 0;
  return 0;
}

/****************************************************************************/
/*!
 **  Define the main options
 **
 ** \return Error return code
 **
 ** \remarks This function initiates the Matelem structures
 **
 *****************************************************************************/
bool MGibbs::_checkSPDE()
{
  S_DECIDE.flag_est     = false;
  S_DECIDE.flag_std     = false;
  S_DECIDE.flag_several = false;
  S_DECIDE.flag_case    = CASE_SIMULATE;
  if (S_DECIDE.flag_est || S_DECIDE.flag_std) S_DECIDE.flag_case = CASE_KRIGING;

  /* Checks */

  // Check the consistency of the constaints Db newly created
  if (!_checkValidAuxiliary()) return false;

  // Initialize the calculation environment
  _calculInitialize();
  S_ENV.ndata = 0;
  S_ENV.ndata1.clear();
  S_ENV.ntarget1.clear();
  _manageMatelem(1);
  Id ncova = _getNcovaWithoutNugget();

  for (Id icov = 0; icov < ncova; icov++)
  {
    SPDE_CURRENT_ICOV = icov;
    _calculUpdate();
    if (_verbose) _printAll("Model (Stationary) Parameters");
  }
  if (S_DECIDE.flag_std && _nvar > 1)
  {
    messerr(
      "Calculation of Kriging Variance is incompatible with Multivariate");
    return false;
  }

  return true;
}

/****************************************************************************/
/*!
 **  Check the pinchout variable
 **
 ** \return Error returned code
 **
 *****************************************************************************/
bool MGibbs::_checkValidPinchout()
{
  VectorDouble tab = _dbout->getColumnByUID(_icolPinch);

  // Check that values are within [0,1] interval

  for (Id iech = 0; iech < _nechout; iech++)
  {
    if (!_dbout->isActive(iech)) continue;
    if (FFFF(tab[iech])) continue;
    if (tab[iech] < 0 || tab[iech] > 1)
    {
      messerr("Pinchout variable should lie in [0,1]");
      messerr("At grid node %d/%d, the value is %lf", iech + 1, _nechout, tab[iech]);
      return false;
    }
  }
  return true;
}

/****************************************************************************/
/*!
 **  Get the elevation within bounds
 **
 ** \return The value assigned to this inequality
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  lower       Lower bound
 ** \param[in]  upper       Upper bound
 **
 *****************************************************************************/
double MGibbs::_drawElevation(M2D_Environ& m2denv,
                              double lower,
                              double upper)
{
  double mean   = m2denv.zmean;
  double stdv   = m2denv.zstdv;
  double lowloc = lower;
  double upploc = upper;
  double value  = 0.;
  if (!FFFF(lower)) lowloc = (lower - mean) / stdv;
  if (!FFFF(upper)) upploc = (upper - mean) / stdv;
  if (!FFFF(lower) && !FFFF(upper))
    value = mean + stdv * law_gaussian_between_bounds(lowloc, upploc);
  else if (FFFF(lower) && FFFF(upper))
    value = mean;
  else if (FFFF(lower))
    value = mean + stdv * law_gaussian_between_bounds(TEST, upploc);
  else if (FFFF(upper))
    value = mean + stdv * law_gaussian_between_bounds(lowloc, TEST);

  return (value);
}

/****************************************************************************/
/*!
 **  Print (concatenate) the printout of an interval
 **
 ** \param[in]  title       Optional title
 ** \param[in]  lower       Lower bound or FFFF
 ** \param[in]  upper       Upper bound or FFFF
 ** \param[in]  tail        0: blank character; 1: "\n" character
 **
 ** \remarks The printed string starts with a blank character
 ** \remarks It ends with either a blank or a <SR/LF> character (see 'tail')
 **
 *****************************************************************************/
void MGibbs::_printConcatenateInterval(const char* title,
                                       double lower,
                                       double upper,
                                       Id tail)
{
  if (title != nullptr) message("%s", title);
  message(" [");
  if (FFFF(lower))
    message("    NA");
  else
    message("%8.4lf", lower);
  message(" ; ");
  if (FFFF(upper))
    message("    NA");
  else
    message("%8.4lf", upper);
  message("]");

  if (tail == 0)
    message(" ");
  else
    message("\n");
}

/****************************************************************************/
/*!
 **  Print the constraints information for a single point
 **
 ** \param[in]  ilayer      Rank of the layer
 ** \param[in]  iech        Rank of the sample
 ** \param[in]  value       Current value
 ** \param[in]  drift       Drift value (or TEST)
 ** \param[in]  vgaus       Current Gaussian value (or TEST)
 ** \param[in]  lower       Lower bound or FFFF
 ** \param[in]  upper       Upper bound or FFFF
 **
 *****************************************************************************/
void MGibbs::_printConstraintsPerPoint(Id ilayer,
                                       Id iech,
                                       double value,
                                       double drift,
                                       double vgaus,
                                       double lower,
                                       double upper)
{
  message("Sample (%d) - Layer (%3d) in", iech + 1, ilayer + 1);
  _printConcatenateInterval(NULL, lower, upper, 0);
  if (!FFFF(drift)) message("- Drift=%8.3lf ", drift);
  if (!(FFFF(value) && FFFF(vgaus)))
  {
    message("->");
    if (FFFF(value))
      message("       NA");
    else
      message(" %8.4lf", value);
    if (!FFFF(vgaus)) message(" (Gaus=%8.4lf)", vgaus);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Check the validity of the Mean and Variance values
 **
 ** \return Error return code
 **
 ** \param[in]  db            Db structure containing the constraints
 ** \param[in]  ilayer        Rank of the layer of interest
 ** \param[in]  iech          Rank of the sample of interest
 ** \param[in]  flag_positive Positivity check
 ** \param[in]  flag_verbose  Verbose output
 ** \param[in]  M             Value for the Mean
 ** \param[in]  S             Value for the Variance
 ** \param[in]  eps           Tolerance
 **
 *****************************************************************************/
bool MGibbs::_checkValidityMS(Db* db,
                              Id ilayer,
                              Id iech,
                              bool flag_positive,
                              bool flag_verbose,
                              double M,
                              double S,
                              double eps)
{
  Id error = 0;
  if (FFFF(M) || FFFF(S)) error = 1;
  if (flag_positive)
  {
    if (M < eps || S < eps) error = 1;
  }
  if (error == 0) return true;

  if (flag_verbose)
  {
    messerr("Error at Sample #%d/%d for Layer #%d", iech + 1,
            db->getNSample(), ilayer + 1);
    if (FFFF(M))
      messerr("- Mean is undefined");
    else
    {
      if (flag_positive && M < eps)
        messerr("- Mean has a too small value (%lf)", M);
    }
    if (FFFF(S))
      messerr("- Variance is undefined");
    else
    {
      if (flag_positive && S < eps)
        messerr("- Variance has a too small value (%lf)", S);
    }
  }
  return false;
}

/****************************************************************************/
/*!
 **  Returns the value of the drift increment at a sample (mean)
 **
 ** \return The mean value or TEST value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  ilayer      Rank of the layer of interest
 ** \param[in]  iech        Rank of the sample of interest
 **
 *****************************************************************************/
double MGibbs::_getM(M2D_Environ& m2denv,
                     Db* db,
                     Id type,
                     Id ilayer,
                     Id iech)
{
  Id iatt = (type == 1) ? m2denv.iatt_fd : m2denv.iatt_fg;
  return db->getArray(iech, iatt + ilayer);
}

/****************************************************************************/
/*!
 **  Returns the value of the gaussian standard deviation
 **
 ** \return The mean value or TEST value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  ilayer      Rank of the layer of interest
 ** \param[in]  iech        Rank of the sample of interest
 **
 *****************************************************************************/
double MGibbs::_getS(M2D_Environ& m2denv,
                     Db* db,
                     Id type,
                     Id ilayer,
                     Id iech)
{
  DECLARE_UNUSED(db);
  DECLARE_UNUSED(type);
  DECLARE_UNUSED(ilayer);
  DECLARE_UNUSED(iech);
  return m2denv.ystdv;
}

/****************************************************************************/
/*!
 **  At a point, returns the external drift increment from previous layer
 **
 ** \return The external drift increment
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  ilayer0     Rank of the layer of interest
 ** \param[in]  iech0       Rank of the sample of interest
 **
 *****************************************************************************/
double MGibbs::_incrementED(M2D_Environ& m2denv,
                            Db* db,
                            Id ilayer0,
                            Id iech0)
{
  double previous;
  double value = db->getLocVariable(ELoc::F, iech0, ilayer0);
  if (FFFF(value)) return (TEST);
  if (ilayer0 > 1)
    previous = db->getLocVariable(ELoc::F, iech0, ilayer0 - 1);
  else
    previous = m2denv.dmini;
  if (FFFF(previous)) return (TEST);
  value -= previous;
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the drift contribution at a sample
 **  This value is a weighted combinaison of constant and external drift term
 **
 ** \return The drift interval value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  ilayer0     Rank of the layer of interest
 ** \param[in]  iech0       Rank of the sample of interest
 **
 *****************************************************************************/
double MGibbs::_getDrift(M2D_Environ& m2denv,
                         Db* db,
                         Id ilayer0,
                         Id iech0) const
{
  double drift;
  double coeff = m2denv.dcoef[ilayer0];
  if (_flagED)
    drift = _incrementED(m2denv, db, ilayer0, iech0);
  else
    drift = 1.;
  if (FFFF(drift)) return (TEST);
  double value = coeff * drift;
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the drift increment in a Db
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  icol_pinch  Pointer to the pinchout variabe
 ** \param[in]  db          Db structure
 ** \param[in]  iatt        Pointer to the drift vector
 **
 *****************************************************************************/
void MGibbs::_setM(M2D_Environ& m2denv,
                   Id icol_pinch,
                   Db* db,
                   Id iatt) const
{
  double drift;
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (db->isActive(iech))
      {
        drift = _getDrift(m2denv, db, ilayer, iech);
        if (!FFFF(drift) && ilayer > 0 && icol_pinch >= 0)
          drift *= db->getArray(iech, icol_pinch);
      }
      else
      {
        drift = TEST;
      }
      db->setArray(iech, iatt + ilayer, drift);
    }
  }
}

/****************************************************************************/
/*!
 **  Locally migrate the pinchout distance from grid to point
 **
 ** \return  Address of the newly added vector in 'dbc'
 **
 *****************************************************************************/
Id MGibbs::_migratePinchToPoint() const
{
  if (_icolPinch < 0) return 0;

  // Add an attribute

  Id iptr = _dbc->addColumnsByConstant(1, TEST);
  if (iptr < 0) return 1;

  // Core allocation

  VectorDouble tab(_nechc);

  // Migrate information from grid to point

  if (migrateByAttribute(_dbout, _dbc, {_icolPinch}, 0, VectorDouble(), false, false))
  {
    _dbc->deleteColumnByUID(iptr);
    return 1;
  }

  // Store the resulting array in the file

  _dbc->setColumnByUID(tab, iptr);

  // Set the error returned code

  return (iptr);
}

/****************************************************************************/
/*!
 **  Calculate and store drift value per point in constraints and output Db
 **  Check the validity of the drift at points
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  mode        1 adding; -1 deleting
 **
 *****************************************************************************/
Id MGibbs::_manageDriftIncrement(M2D_Environ& m2denv, Id mode)
{
  Id iptr = -1;

  /* Dispatch */

  if (mode > 0)
  {

    /* Identify the drift at the constraining samples */

    m2denv.iatt_fd = _dbc->addColumnsByConstant(_nlayer, TEST);
    if (m2denv.iatt_fd < 0) return (1);

    /* If pinch-out is defined, interpolate it at well data */

    iptr = _migratePinchToPoint();
    _setM(m2denv, iptr, _dbc, m2denv.iatt_fd);
    if (iptr >= 0) _dbc->deleteColumnByUID(iptr);

    /* Check validity of drift at data points */

    for (Id iech = 0; iech < _dbc->getNSample(); iech++)
    {
      if (!_dbc->isActive(iech)) continue;
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {
        double M = _getM(m2denv, _dbc, 1, ilayer, iech);
        double S = _getS(m2denv, _dbc, 1, ilayer, iech);
        if (!_checkValidityMS(_dbc, ilayer, iech, true, true, M, S)) return (1);
      }
    }

    /* Identify the drift at the target grid nodes */

    m2denv.iatt_fg = _dbout->addColumnsByConstant(_nlayer, TEST);
    if (m2denv.iatt_fg < 0) return (1);
    _setM(m2denv, _icolPinch, _dbout, m2denv.iatt_fg);
  }
  else
  {

    /* Deleting the drift at the constraining samples */

    if (m2denv.iatt_fd >= 0)
      _dbc->deleteColumnsByUIDRange(m2denv.iatt_fd, _nlayer);

    /* Deleting the drift at the target grid */

    if (m2denv.iatt_fg >= 0)
      _dbout->deleteColumnsByUIDRange(m2denv.iatt_fg, _nlayer);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate global statistics on elevations
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  percent     Extension percentage for min and max values
 **
 *****************************************************************************/
void MGibbs::_statsInit(M2D_Environ& m2denv, double percent)
{
  double nb   = 0.;
  double mm   = 0.;
  double vv   = 0.;
  double mini = MAXIMUM_BIG;
  double maxi = MINIMUM_BIG;

  /* Loop on the layers */

  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {

    /* Loop on the samples */

    for (Id iech = 0; iech < _nechin; iech++)
    {
      if (!_dbin->isActive(iech)) continue;
      double lower = _dbin->getLocVariable(ELoc::L, iech, ilayer);
      double upper = _dbin->getLocVariable(ELoc::U, iech, ilayer);

      // Process the minimum bound

      if (!FFFF(lower))
      {
        nb += 1.;
        mm += lower;
        vv += lower * lower;
        if (lower < mini) mini = lower;
        if (lower > maxi) maxi = lower;
      }

      // Process the maximum bound

      if (!FFFF(upper))
      {
        nb += 1.;
        mm += upper;
        vv += upper * upper;
        if (upper < mini) mini = upper;
        if (upper > maxi) maxi = upper;
      }
    }
  }

  /* Normation */

  if (nb > 0)
  {
    mm /= nb;
    vv = vv / nb - mm * mm;
  }
  else
  {
    mm   = 0.;
    vv   = 1.;
    mini = 0.;
    maxi = 1.;
  }

  double delta = maxi - mini;
  if (delta <= 0) delta = ABS(mm) / 10.;
  if (delta <= 0) delta = 1.;
  m2denv.zmean = mm;
  m2denv.zeps  = ABS(mm) / 1.e4;
  m2denv.zstdv = (vv > 0) ? sqrt(vv) : 1.;
  m2denv.zmini = mini - delta * percent;
  m2denv.zmaxi = maxi + delta * percent;

  if (_verbose)
  {
    mestitle(2, "Global Statistics on Raw Elevations (extended by %4.2lf)", percent);
    message("Statistics are derived from compiling bounds (when defined)\n");
    message("Number of valid bounds = %d\n", static_cast<Id>(nb));
    message("Mean                   = %lf\n", m2denv.zmean);
    message("St. Deviation          = %lf\n", m2denv.zstdv);
    message("Tolerance              = %lf\n", m2denv.zeps);
    message("Minimum                = %lf\n", m2denv.zmini);
    message("Maximum                = %lf\n", m2denv.zmaxi);
    message("Range                  = %lf\n", m2denv.zmaxi - m2denv.zmini);
    if (nb <= 0) message("(Range of values are fixed arbitrarily as no bound is found)\n");
  }
}

/****************************************************************************/
/*!
 **  Update global statistics on the raw information
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  percent     Tolerance
 **
 *****************************************************************************/
void MGibbs::_statsUpdate(M2D_Environ& m2denv, double percent) const
{
  double nb   = 0.;
  double mm   = 0.;
  double vv   = 0.;
  double mini = MAXIMUM_BIG;
  double maxi = MINIMUM_BIG;

  /* Loop on the layers */

  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {

    /* Loop on the samples */

    for (Id iech = 0; iech < _nechc; iech++)
    {
      double zval = _dbc->getZVariable(iech, ilayer);

      nb += 1.;
      mm += zval;
      vv += zval * zval;
      if (zval < mini) mini = zval;
      if (zval > maxi) maxi = zval;
    }
  }

  /* Normation */

  if (nb > 0)
  {
    mm /= nb;
    vv = vv / nb - mm * mm;
  }
  else
  {
    mm   = 0.;
    vv   = 1.;
    mini = -0.5;
    maxi = 0.5;
  }

  double delta = maxi - mini;
  if (delta <= 0.) delta = ABS(mm) / 10.;
  if (delta <= 0) delta = 1.;
  m2denv.zmean = mm;
  m2denv.zeps  = ABS(mm) / 1.e4;
  m2denv.zstdv = (vv > 0) ? sqrt(vv) : 1.;
  m2denv.zmini = mini - delta * percent;
  m2denv.zmaxi = maxi + delta * percent;

  if (_verbose)
  {
    mestitle(2, "Global Statistics on Centered Elevations");
    message("Statistics are compiled from initial values within bounds\n");
    message("Number of values = %d\n", static_cast<Id>(nb));
    message("Mean             = %lf\n", m2denv.zmean);
    message("St. Deviation    = %lf\n", m2denv.zstdv);
    message("Tolerance        = %lf\n", m2denv.zeps);
    message("Minimum          = %lf\n", m2denv.zmini);
    message("Maximum          = %lf\n", m2denv.zmaxi);
    message("Range            = %lf\n", m2denv.zmaxi - m2denv.zmini);
  }
}

/****************************************************************************/
/*!
 **  Set the initial elevations at the constraining information
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  njter_max   Maximum number of iterations
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 ** \remarks This function also add the attributes to 'dbin' per layer:
 ** \remarks - the initial value (ELoc::Z)
 **
 *****************************************************************************/
Id MGibbs::_initializeElevations(M2D_Environ& m2denv,
                                 VectorDouble& work,
                                 Id njter_max) const
{
  double eps   = m2denv.zeps;
  Id flag_jter = 0;

  /* Loop on the samples */

  for (Id iech = 0; iech < _nechc; iech++)
  {

    /* Define the values at sample as unconstrained information */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      double zmin  = _dbc->getLocVariable(ELoc::L, iech, ilayer);
      double zmax  = _dbc->getLocVariable(ELoc::U, iech, ilayer);
      work[ilayer] = _drawElevation(m2denv, zmin, zmax);
    }

    /* Loop on iterations for ordering the values */

    for (Id jter = 0; jter < njter_max; jter++)
    {
      flag_jter = 0;

      /* Loop on the layers */

      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {

        /* Determine the bounds at data locations */

        double zmin = _dbc->getLocVariable(ELoc::L, iech, ilayer);
        double zmax = _dbc->getLocVariable(ELoc::U, iech, ilayer);

        /* Loop on the other layers */

        for (Id jlayer = 0; jlayer < _nlayer; jlayer++)
        {
          if (ilayer == jlayer) continue;
          double zval = work[jlayer];

          if (jlayer < ilayer)
          {

            // Comparing with a layer located shallower than the current one

            if (!FFFF(zmax) && zval > zmax)
              flag_jter = 1;
            else
            {
              if (FFFF(zmin))
                zmin = zval + eps;
              else
                zmin = MAX(zmin, zval);
            }
          }
          else
          {

            // Comparing with a layer located deeper than the current one

            if (!FFFF(zmin) && zval < zmin)
              flag_jter = 1;
            else
            {
              if (FFFF(zmax))
                zmax = zval - eps;
              else
                zmax = MIN(zmax, zval);
            }
          }
        }

        // Update target value according to constraints

        work[ilayer] = _drawElevation(m2denv, zmin, zmax);
      }

      /* Interrupt iterations */

      if (!flag_jter) break;
    }

    // Run abort in case of lack of convergence

    if (flag_jter)
    {
      messerr("At constraining sample #%d/%d, correct interval ordering",
              iech + 1, _nechc);
      messerr("has not been reached after %d iterations. Run is aborted",
              njter_max);
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {
        double zmin = _dbc->getLocVariable(ELoc::L, iech, ilayer);
        double zmax = _dbc->getLocVariable(ELoc::U, iech, ilayer);
        _printConstraintsPerPoint(ilayer, iech, work[ilayer],
                                  TEST,
                                  TEST, zmin, zmax);
      }
      messerr("\n");
      messerr(">>> You should check the ordering of your bound variables");
      messerr(">>>in the Well File");
      messerr("\n");
      return (1);
    }

    /* Store the resulting values */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      _dbc->setLocVariable(ELoc::Z, iech, ilayer, work[ilayer]);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Fit the coefficients of the trend terms for each layer
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  percent     Tolerance
 **
 ** \param[out] iatt_f      Pointer in dbin to the added variables ELoc::F
 **
 ** \remarks This function also add the attributes to 'dbin' per layer:
 ** \remarks - the external drift values (ELoc::F)
 **
 *****************************************************************************/
Id MGibbs::_manageDrift(M2D_Environ& m2denv, Id* iatt_f, double percent)
{
  (*iatt_f) = -1;

  /* Core allocation */

  VectorDouble dval;
  if (_flagED)
  {
    dval.resize(_nechin);
  }

  /* Add attributes to 'dbin' */
  /* - the external drift value at data points (optional) */
  /* - the initial value at data points */

  if (_flagED)
  {
    if (db_locator_attribute_add(_dbin, ELoc::F, _nlayer, 0, TEST, iatt_f))
      return 1;
  }

  /* Loop on the layers */

  Id nb = 0;
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {

    /* Export External drift from 'dbout' to 'dbin' (optional) */

    if (_flagED)
    {
      Id colrank = _dbout->getColIdxByLocator(ELoc::F, ilayer);

      // Migrate the information from Grid to Wells

      migrateByAttribute(_dbout, _dbin, {colrank}, 0, VectorDouble(), false, false);

      // Calculate the statistics of the external drift on the grid

      for (Id iech = 0; iech < _dbout->getNSample(); iech++)
      {
        if (!_dbout->isActive(iech)) continue;
        double value = _dbout->getLocVariable(ELoc::F, iech, ilayer);
        if (FFFF(value)) continue;
        nb++;
        if (FFFF(m2denv.dmini) || value < m2denv.dmini) m2denv.dmini = value;
        if (FFFF(m2denv.dmaxi) || value > m2denv.dmaxi) m2denv.dmaxi = value;
      }
    }

    /* Loop on the samples */

    for (Id iech = 0; iech < _nechin; iech++)
    {
      if (!_dbin->isActive(iech)) continue;
      if (_flagED)
      {
        if (FFFF(dval[iech])) continue;
        _dbin->setLocVariable(ELoc::F, iech, ilayer, dval[iech]);
      }
    }
  }

  /* Patch the statistics on drift if no external drift */

  if (!_flagED)
  {
    m2denv.dmini = 0.;
    m2denv.dmaxi = 1.;
  }
  else
  {
    double delta = m2denv.dmaxi - m2denv.dmini;
    m2denv.dmini -= delta * percent;
    m2denv.dmaxi += delta * percent;
  }

  if (_verbose)
  {
    mestitle(2, "Global Statistics on Trends (extended by %4.2lf)", percent);
    message("Statistics are derived from compiling drift at grid nodes\n");
    message("Number of valid nodes  = %d\n", static_cast<Id>(nb));
    message("Minimum Drift          = %lf\n", m2denv.dmini);
    message("Maximum Drift          = %lf\n", m2denv.dmaxi);
    message("Range of Drift         = %lf\n", m2denv.dmaxi - m2denv.dmini);
  }

  return 0.;
}

/****************************************************************************/
/*!
 **  Print the details of the constraints
 **
 ** \param[in]  nech        Number of hard data
 ** \param[in]  ilayer      Rank of the target layer
 **
 *****************************************************************************/
void MGibbs::_printDetails(Id nech, Id ilayer) const
{
  Id nvalue = 0;
  Id nbdmin = 0;
  Id nbdmax = 0;
  for (Id iech = 0; iech < nech; iech++)
  {
    double value = _dbc->getZVariable(iech, ilayer);
    if (!FFFF(value)) nvalue++;
    double lower = _dbc->getLocVariable(ELoc::L, iech, ilayer);
    double upper = _dbc->getLocVariable(ELoc::U, iech, ilayer);
    if (!FFFF(lower)) nbdmin++;
    if (!FFFF(upper)) nbdmax++;
  }

  // Printout

  message("  . Number of hard data    = %d\n", nvalue);
  message("  . Number of lower limits = %d\n", nbdmin);
  message("  . Number of upper limits = %d\n", nbdmax);
}

/****************************************************************************/
/*!
 **  Fit the coefficients of the trend terms for each layer
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  number_hard Number of hard data used to fit the drift
 **
 ** \remarks The drift is only established on data where lower and upper bounds
 ** \remarks are both defined. The drift coefficients are assumed to be the same
 ** \remarks for all layers
 ** \remarks The impact of areal constraints being to important, it has been
 ** \remarks chosen to base the drift fitting only on the first 'number_hard'
 ** \remarks samples (which correspond to constraints coming from 'dbin'.
 **
 *****************************************************************************/
Id MGibbs::_driftFitting(M2D_Environ& m2denv, Id number_hard) const
{
  Id nech = MIN(number_hard, _dbc->getNSample());
  Id nbfl = 1;
  m2denv.dcoef.resize(_nlayer);
  VectorDouble a(nbfl * nbfl);
  VectorDouble b(nbfl);

  /* Loop on the layers */

  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {

    /* Initializations */

    Id numb       = 0;
    double mean   = 0.;
    double ffmean = 0.;
    double stdv   = 0.;
    double mini   = MAXIMUM_BIG;
    double maxi   = MINIMUM_BIG;
    double ffmini = MAXIMUM_BIG;
    double ffmaxi = MINIMUM_BIG;
    for (Id i = 0; i < nbfl; i++)
      b[i] = 0.;
    for (Id i = 0; i < nbfl * nbfl; i++)
      a[i] = 0.;

    /* Loop on the samples */

    for (Id iech = 0; iech < nech; iech++)
    {

      /* Get the values at the data point */

      double epais = _dbc->getZVariable(iech, ilayer);
      if (ilayer > 0)
        epais -= _dbc->getZVariable(iech, ilayer - 1);
      else
        epais -= m2denv.zmini;

      /* Set the drift vector at data point */

      double ff = 1.;
      if (_flagED)
        ff = _incrementED(m2denv, _dbc, ilayer, iech);

      /* Update statistics */

      numb += 1;
      mean += epais;
      stdv += epais * epais;
      ffmean += ff;
      if (epais < mini) mini = epais;
      if (epais > maxi) maxi = epais;
      if (ff < ffmini) ffmini = ff;
      if (ff > ffmaxi) ffmaxi = ff;

      /* Fill the linear system */

      b[0] = ff * epais;
      a[0] = ff * ff;
    }

    /* Save the results */

    m2denv.dcoef[ilayer] = b[0] / a[0];

    /* Normalize statistics */

    if (numb > 0)
    {
      mean /= numb;
      ffmean /= numb;
      stdv = stdv / numb - mean * mean;
      stdv = (stdv > 0) ? sqrt(stdv) : 0.;
    }

    /* Print statistics (optional) */

    if (_verbose)
    {
      message("\nLayer #%d\n", ilayer + 1);
      message("- Number of Constraints = %d \n", numb);
      _printDetails(nech, ilayer);
      message("- Drift:\n");
      if (_flagED)
      {
        message("  . Mean          = %lf\n", ffmean);
        message("  . Minimum       = %lf\n", ffmini);
        message("  . Maximum       = %lf\n", ffmaxi);
      }
      message("  . Coefficient   = %lg\n", m2denv.dcoef[ilayer]);
      message("- Residual:\n");
      message("  . Mean          = %lf\n", mean);
      message("  . St. Deviation = %lf\n", stdv);
      message("  . Minimum       = %lf\n", mini);
      message("  . Maximum       = %lf\n", maxi);
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Save the drift at the grid nodes
 **
 ** \param[in]  m2denv      M2D_Environ structure
 **
 ** \param[out] gwork       Working array
 **
 ** \remarks The drift returned as a surface uses directly the coefficients
 ** \remarks of the linear combinaison.
 **
 *****************************************************************************/
void MGibbs::_driftSave(M2D_Environ& m2denv, VectorDouble& gwork)
{
  for (Id igrid = 0; igrid < _nechout; igrid++)
  {
    if (!_dbout->isActive(igrid)) continue;
    double drift = m2denv.zmini;

    /* Loop no the layers */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      double value = _getDrift(m2denv, _dbout, ilayer, igrid);
      if (FFFF(value))
        drift = TEST;
      else
        drift += value;
      GWORK(ilayer, igrid) = drift;
    }
  }
}

/****************************************************************************/
/*!
 **  Record a new active point
 **
 ** \param[in]  db          Db structure
 ** \param[in]  iech        Sample rank in 'db'
 ** \param[in]  natt        Number of attributes
 **
 ** \param[in,out] number   Number of samples
 ** \param[in,out] tab      Array of samples
 **
 ** \return 1 if:
 ** \return - the sample is masked off
 ** \return - if no bound is defined for current sample (any layer)
 **
 *****************************************************************************/
bool MGibbs::_recordSample(Db* db,
                           Id iech,
                           Id natt,
                           Id& number,
                           VectorDouble& tab) const
{
  if (!db->isActive(iech)) return false;

  // Perform the different assignments

  Id ecr = number * natt;

  // Set the rank

  tab[ecr++] = static_cast<double>(number) + 1;

  // Set the coordinates

  for (Id idim = 0; idim < _ndim; idim++)
  {
    double coor = db->getCoordinate(iech, idim);
    tab[ecr++]  = coor;
  }

  // For each layer, set the bounds and the initial value

  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    double value = db->getLocVariable(ELoc::Z, iech, ilayer);
    double lower = db->getLocVariable(ELoc::L, iech, ilayer);
    double upper = db->getLocVariable(ELoc::U, iech, ilayer);
    // Check that the bounds (if defined) are correctly ordered
    if (!FFFF(lower) && !FFFF(upper) && lower > upper)
    {
      messerr("For sample %d at layer %d, Lower bound (%lf) > Upper bound (%lf)",
              number, ilayer + 1, lower, upper);
      messerr("Bounds constraints are discarded");
      lower = TEST;
      upper = TEST;
    };
    // If hard data is defined, convert into Lower and Upper bounds
    if (!FFFF(value))
    {
      lower = value;
      upper = value;
    }
    // Check that at least one bound is defined
    if (FFFF(lower) && FFFF(upper)) return false;
    tab[ecr++] = lower;
    tab[ecr++] = upper;
    tab[ecr++] = value;
  }

  // For each layer, set the External Drift value (optional)

  if (_flagED)
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      tab[ecr++] = db->getLocVariable(ELoc::F, iech, ilayer);

  /* Increment the number of records by 1 */

  number++;
  return true;
}

/****************************************************************************/
/*!
 **  Define the locators on the newly created Db
 **
 ** \param[in]  db          Db constraints structure
 ** \param[in]  nvarloc     Number of variables
 **
 *****************************************************************************/
void MGibbs::_defineLocators(Db* db, Id nvarloc) const
{
  Id rank = 1;
  db->setLocatorsByUID(_ndim, rank, ELoc::X, 0);
  rank += _ndim;
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    db->setLocatorByUID(rank++, ELoc::L, ilayer);
    db->setLocatorByUID(rank++, ELoc::U, ilayer);
    if (ilayer < nvarloc) db->setLocatorByUID(rank, ELoc::Z, ilayer);
    rank++;
  }
  if (_flagED) db->setLocatorsByUID(_nlayer, rank, ELoc::F, 0);
}

/****************************************************************************/
/*!
 **  Print the Environnement
 **
 *****************************************************************************/
void MGibbs::_printEnviron(const char* title, M2D_Environ& m2denv) const
{
  mestitle(1, title);

  if (_flagED)
    message("Use of External Drift\n");
  else
    message("No External Drift\n");
  message("Z Minimum               = %lf\n", m2denv.zmini);
  message("Z Maximum               = %lf\n", m2denv.zmaxi);
  message("Z Mean                  = %lf\n", m2denv.zmean);
  message("Z St. Deviation         = %lf\n", m2denv.zstdv);
  message("Z Tolerance             = %lf\n", m2denv.zeps);
  message("Drift Minimum           = %lf\n", m2denv.dmini);
  message("Drift Maximum           = %lf\n", m2denv.dmaxi);
  message("Y St. Deviation         = %lf\n", m2denv.ystdv);
}

/****************************************************************************/
/*!
 **  Create a Db containing all the constraining information
 **
 ** \return  Error return code
 **
 ** \param[out] number_hard Number of hard data which will serve for
 **                         seting the optimal drift
 **
 ** \remark Note that the file constructed here contains as many samples
 ** \remark as the number of ACTIVE samples of the input Db
 **
 *****************************************************************************/
Id MGibbs::_createConstraints(Id* number_hard)
{
  Id nechtot = _nechin + _nechout;
  Id natt    = 1;               // Rank
  natt += _ndim;                // Coordinates
  natt += 3 * _nlayer;          // LowBound, UppBound and Variable per layer
  if (_flagED) natt += _nlayer; // External Drift

  /* Core allocation */

  VectorDouble tab(nechtot * natt);

  /* Load information from 'dbin' */

  Id number = 0;
  for (Id iech = 0; iech < _nechin; iech++)
  {
    if (!_recordSample(_dbin, iech, natt, number, tab)) continue;
  }
  *number_hard = number;

  /* Load information from 'dbout' */

  for (Id iech = 0; iech < _nechout; iech++)
  {
    if (!_recordSample(_dbout, iech, natt, number, tab)) continue;
  }

  /* When forcing, the first active sample is used */

  if (number <= 0)
  {
    for (Id iech = 0; iech < _nechin; iech++)
    {
      if (!_recordSample(_dbin, iech, natt, number, tab)) continue;
      if (number > 0) break;
    }
  }

  /* Core reallocation */

  if (number < nechtot) tab.resize(number * natt);

  // Creating names for the variables

  VectorString names;
  names.push_back("Rank");
  for (Id idim = 0; idim < _ndim; idim++)
  {
    (void)gslSPrintf(string_encode, "x-%d", idim + 1);
    names.push_back(string_encode);
  }
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    (void)gslSPrintf(string_encode, "Lower-%d", ilayer + 1);
    names.push_back(string_encode);
    (void)gslSPrintf(string_encode, "Upper-%d", ilayer + 1);
    names.push_back(string_encode);
    (void)gslSPrintf(string_encode, "z-%d", ilayer + 1);
    names.push_back(string_encode);
  }
  if (_flagED)
  {
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      (void)gslSPrintf(string_encode, "Drift%d", ilayer + 1);
      names.push_back(string_encode);
    }
  }

  /* Create the output Db */

  _dbc = Db::createFromSamples(number, ELoadBy::SAMPLE, tab, names,
                               VectorString(), 0);
  if (_dbc == nullptr) return 1;

  _nechc = _dbc->getNSample();
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the inverse of (s2 * Q + B %*% Bt) and store it
 **  into a new QChol object
 **
 ** \return  The calculated sparse matrix
 **
 ** \param[in]  s2          Nugget effect value
 ** \param[in]  Qc          Qc structure (already existing)
 ** \param[in]  Matelem     Matelem structure
 **
 *****************************************************************************/
QChol* MGibbs::_deriveQC(double s2, QChol* Qc, SPDE_Matelem& Matelem)
{
  Id error         = 1;
  MatrixSparse* Bt = nullptr;
  MatrixSparse* B2 = nullptr;
  MatrixSparse* Q  = Matelem.QC->Q;
  MatrixSparse* B  = Matelem.Aproj;

  // Clean the previous Qc (if it exists)

  if (Qc != nullptr) Qc = _manageQchol(-1, Qc);

  // Calculate: Q + t(B) %*% B

  if (DEBUG)
    message("Building Q (Size:%d) with additional nugget effect (%lf)\n",
            Q->getNCols(), s2);
  Bt = B->transpose();
  if (Bt == nullptr) goto label_end;
  B2 = MatrixFactory::prodMatMat<MatrixSparse>(Bt, B);
  if (B2 == nullptr) goto label_end;

  Qc = _manageQchol(1, NULL);
  if (Qc == nullptr) goto label_end;
  Qc->Q = MatrixSparse::addMatMat(Q, B2, s2, 1.);
  if (Qc->Q == nullptr) goto label_end;

  // Perform the Cholesky transform

  error = _qcholCholesky(false, Qc);

  // Free memory

label_end:
  delete Bt;
  delete B2;
  if (error) Qc = _manageQchol(-1, Qc);
  return (Qc);
}

/****************************************************************************/
/*!
 **  Draw a Z-value within bounds
 **
 ** \return Z-Value in the working domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  iter        Rank of the iteration
 ** \param[in]  ilayer      Rank of the layer
 ** \param[in]  iech        Rank of the sample
 ** \param[in]  Zval        Input value
 ** \param[in]  Zcum        Cumulated Z value of layers above
 ** \param[in]  Zmin        Lower bound in Z
 ** \param[in]  Zmax        Upper bound in Z
 ** \param[in]  Ymean       Mean of the Y Law
 ** \param[in]  Ysigma      Standard deviation of the Y Law
 **
 *****************************************************************************/
double MGibbs::_drawGaussian(M2D_Environ& m2denv,
                             bool verbose,
                             Id iter,
                             Id ilayer,
                             Id iech,
                             double Zval,
                             double Zcum,
                             double Zmin,
                             double Zmax,
                             double Ymean,
                             double Ysigma)
{
  double Ymin, Ymax;
  bool flagVerif = true;

  /* Initializations */

  double M = _getM(m2denv, _dbc, 1, ilayer, iech);
  double S = _getS(m2denv, _dbc, 1, ilayer, iech);
  if (!_checkValidityMS(_dbc, ilayer, iech, true, true, M, S))
    messageAbort("- Impossible to have M or S undefined");

  if (verbose)
  {
    message("Input Z elevation=%lf", Zval);
    _printConcatenateInterval(NULL, Zmin, Zmax, 1);
  }

  /* Centering in Z */

  double Zminc = Zmin;
  if (!FFFF(Zminc)) Zminc -= Zcum;
  if (Zminc < 0) Zminc = 0.;
  double Zmaxc = Zmax;
  if (!FFFF(Zmaxc)) Zmaxc -= Zcum;
  if (Zmaxc < 0) Zmaxc = 0.;
  if (verbose) _printConcatenateInterval("Z thickness", Zminc, Zmaxc, 1);

  /* Converting from Z to Y */

  if (FFFF(Zminc))
    Ymin = TEST;
  else
    Ymin = (Zminc == 0) ? TEST : (S * S / 2. + log(Zminc / M)) / S;
  if (FFFF(Zmaxc))
    Ymax = TEST;
  else
    Ymax = (Zmaxc == 0) ? TEST : (S * S / 2. + log(Zmaxc / M)) / S;
  if (verbose) _printConcatenateInterval("Y gaussian", Ymin, Ymax, 1);

  /* Centering in Y */

  if (!FFFF(Ymin)) Ymin = (Ymin - Ymean) / Ysigma;
  if (!FFFF(Ymax)) Ymax = (Ymax - Ymean) / Ysigma;
  if (verbose) _printConcatenateInterval("Y centered", Ymin, Ymax, 1);

  double Yval = law_gaussian_between_bounds(Ymin, Ymax);
  // Two next lines are there for robustification: they should not be removed
  if (!FFFF(Ymin) && Yval < Ymin) Yval = Ymin;
  if (!FFFF(Ymax) && Yval > Ymax) Yval = Ymax;

  Yval = Ymean + Ysigma * Yval;
  Zval = Zcum + M * exp(S * Yval - S * S / 2.);

  if (flagVerif)
  {
    if (std::isinf(Zval))
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message("  Zval=Inf");
      _printConcatenateInterval(NULL, Zmin, Zmax, 1);
      messageAbort("Strange output value for Zval");
    }
    if (!FFFF(Zmin) && Zval < Zmin - m2denv.zeps)
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message(" Zval=%lf", Zval);
      _printConcatenateInterval(NULL, Zmin, Zmax, 1);
      message(" Yval=%lf", Yval);
      _printConcatenateInterval(NULL, Ymin, Ymax, 1);
      messageAbort("Zval should not be smaller than Zmin");
    }
    if (!FFFF(Zmax) && Zval > Zmax + m2denv.zeps)
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message(" Zval=%lf", Zval);
      _printConcatenateInterval(NULL, Zmin, Zmax, 1);
      message(" Yval=%lf", Yval);
      _printConcatenateInterval(NULL, Ymin, Ymax, 1);
      messageAbort("Zval should not be larger than Zmax");
    }
  }

  if (verbose)
  {
    message("Output Z elevation=%lf in", Zval);
    _printConcatenateInterval(NULL, Zmin, Zmax, 1);
  }

  return (Zval);
}

/****************************************************************************/
/*!
 **  Convert a layer-pile at a datum from the working to the true domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  iech        Rank of the sample
 ** \param[in,out] tab      Input/Output array of Z-values (Dimension: nlayer)
 **
 *****************************************************************************/
void MGibbs::_converZ2Y(M2D_Environ& m2denv,
                        Id type,
                        Id iech,
                        VectorDouble& tab) const
{
  double Yval, Zval, Zcur;

  Id flag_undef = 0;
  double Zcum   = m2denv.zmini;
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    double M = _getM(m2denv, _dbc, type, ilayer, iech);
    double S = _getS(m2denv, _dbc, type, ilayer, iech);
    if (!_checkValidityMS(_dbc, ilayer, iech, true, true, M, S) || flag_undef)
    {
      flag_undef = 1;
      Yval       = TEST;
    }
    else
    {
      Zcur = tab[ilayer];
      Zval = Zcur - Zcum;
      if (Zval <= 0)
      {
        flag_undef = 1;
        Yval       = TEST;
      }
      else
        Yval = (S * S / 2. + log(Zval / M)) / S;
      Zcum = Zcur;
    }
    tab[ilayer] = Yval;
  }
}

/****************************************************************************/
/*!
 **  Convert a layer-pile at a datum from the true to the working domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  iech        Rank of the sample
 ** \param[in,out] tab      Input/Ouput array of Y-values (Dimension: nlayer)
 **
 *****************************************************************************/
void MGibbs::_convertY2Z(M2D_Environ& m2denv,
                         Db* db,
                         Id type,
                         Id iech,
                         VectorDouble& tab) const
{
  double Zval, Yval;

  Id flag_undef = 0;
  double Zcur   = m2denv.zmini;
  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    double M = _getM(m2denv, db, type, ilayer, iech);
    double S = _getS(m2denv, db, type, ilayer, iech);
    if (!_checkValidityMS(db, ilayer, iech, false, false, M, S) || flag_undef)
    {
      flag_undef = 1;
      Zcur       = TEST;
    }
    else
    {
      Yval = tab[ilayer];
      Zval = M * exp(S * Yval - S * S / 2.);
      Zcur += Zval;
    }
    tab[ilayer] = Zcur;
  }
}

/****************************************************************************/
/*!
 **  Print the values at a sample location
 **
 ** \param[in]  title       Title
 ** \param[in]  iech        Sample rank
 ** \param[in]  work        Array of values (defined in Z)
 **
 *****************************************************************************/
void MGibbs::_printSample(const char* title,
                          Id iech,
                          VectorDouble& work) const
{
  message("%s - Sample #%d/%d\n", title, iech + 1, _nechc);

  for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
  {
    double zmin = _dbc->getLocVariable(ELoc::L, iech, ilayer);
    double zmax = _dbc->getLocVariable(ELoc::U, iech, ilayer);
    message("Z(%d)=%lf in", ilayer, work[ilayer]);
    _printConcatenateInterval(NULL, zmin, zmax, 1);
  }
}

/****************************************************************************/
/*!
 **  Perform the Gibbs iterations
 **
 ** \return Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  iter        Rank of the iteration
 ** \param[in]  sigma       Standard deviation of the nugget value
 ** \param[in]  ymean       Array of mean values at constraints
 ** \param[in,out] ydat     Array of values at constraints samples
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 ** \remarks The [in,out] argument 'ydat' is expressed in working domain
 ** \remarks It needs to be locally transformed in real domain in order to
 ** \remarks be compared to the bounds.
 **
 *****************************************************************************/
Id MGibbs::_gibbs(M2D_Environ& m2denv,
                  Id iter,
                  double sigma,
                  VectorDouble& ymean,
                  VectorDouble& ydat,
                  VectorDouble& work)
{
  bool local_verbose = false;

  // Loop on the samples

  for (Id iech = 0; iech < _nechc; iech++)
  {

    // Set the initial values

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      work[ilayer] = YDAT(_nechc, ilayer, iech);
    _convertY2Z(m2denv, _dbc, 1, iech, work);
    if (local_verbose)
      _printSample("Entering in Gibbs", iech, work);

    // Loop on the layers

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {

      // Get the elevation of the previous layer

      double zcum = m2denv.zmini;

      // Getting the elevation and the bounds for the current layer

      double zmin = _dbc->getLocVariable(ELoc::L, iech, ilayer);
      double zmax = _dbc->getLocVariable(ELoc::U, iech, ilayer);
      if (local_verbose)
      {
        message("ilayer=%d", ilayer);
        _printConcatenateInterval(NULL, zmin, zmax, 1);
      }

      // Loop on the other layers

      for (Id jlayer = 0; jlayer < _nlayer; jlayer++)
      {
        if (ilayer == jlayer) continue;
        double zval = work[jlayer];
        if (local_verbose)
          message("Constrained by jlayer=%d zval=%lf\n", jlayer, zval);

        if (jlayer < ilayer)
        {

          // Comparing with a layer located shallower than the current one

          if (FFFF(zmin))
            zmin = zval;
          else
            zmin = MAX(zmin, zval);
          zcum = zval;
        }
        else
        {

          // Comparing with a layer located deeper than the current one

          if (FFFF(zmax))
            zmax = zval;
          else
            zmax = MIN(zmax, zval);
        }
      }

      // Drawing plausible values according to constraints

      work[ilayer] = _drawGaussian(m2denv, local_verbose, iter, ilayer,
                                   iech, work[ilayer], zcum, zmin, zmax,
                                   YMEAN(ilayer, iech), sigma);
    }

    // Load the new values

    if (local_verbose)
      _printSample("Exiting Gibbs", iech, work);
    _converZ2Y(m2denv, 1, iech, work);

    /* Store in the extracted vector */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      YDAT(_nechc, ilayer, iech) = work[ilayer];
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check that the Gibbs constraints are fullfilled at datum locations
 **
 ** \return  Error return code
 **
 ** \param[in]  title       Title for the printout (if error)
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  ydat        Array of simulations on the data
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 *****************************************************************************/
bool MGibbs::_checkGibbsData(const char* title,
                             M2D_Environ& m2denv,
                             VectorDouble& ydat,
                             VectorDouble& work)
{
  Id error = 0;

  // Loop on the constraints samples

  for (Id iech = 0; iech < _nechc; iech++)
  {
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      work[ilayer] = YDAT(_nechc, ilayer, iech);
    _convertY2Z(m2denv, _dbc, 1, iech, work);

    // Loop on the layers

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      double depth = work[ilayer];

      // Getting the elevation and the bounds for the current layer

      double zmin = _dbc->getLocVariable(ELoc::L, iech, ilayer);
      double zmax = _dbc->getLocVariable(ELoc::U, iech, ilayer);

      // Check consistency

      if (!FFFF(zmin))
      {
        if (depth < zmin - m2denv.zeps)
        {
          messerr("%s: Sample(%d/%d) of Layer(%d/%d): Depth(%lf) < Lower(%lf)",
                  title, iech + 1, _nechc, ilayer + 1, _nlayer, depth, zmin);
          error++;
        }
      }
      if (!FFFF(zmax))
      {
        if (depth > zmax + m2denv.zeps)
        {
          messerr("%s: Sample(%d/%d) of Layer(%d/%d): Depth(%lf) > Upper(%lf)",
                  title, iech + 1, _nechc, ilayer + 1, _nlayer, depth, zmax);
          error++;
        }
      }
    }
  }

  if (_verbose && error > 0)
  {
    message("%s: %d error(s) found\n", title, error);
  }
  return (error == 0);
}

/****************************************************************************/
/*!
 **  Manage the M2D_Environ structure
 **
 ** \return  Pointer to the M2D_Environ structure
 **
 ** \param[in]  ystdv       Stamdard deviation of the Gaussian Transformed
 ** \param[in]  m2denv      M2D_Environ structure (if NULL, a new one is created)
 **
 *****************************************************************************/
void MGibbs::_manageM2denv(double ystdv, M2D_Environ& m2denv)
{
  m2denv.iatt_fd = -1;
  m2denv.iatt_fg = -1;
  m2denv.zmean   = 0.;
  m2denv.zstdv   = 1.;
  m2denv.zeps    = 0.;
  m2denv.zmini   = TEST;
  m2denv.zmaxi   = TEST;
  m2denv.dmini   = TEST;
  m2denv.dmaxi   = TEST;
  m2denv.ystdv   = ystdv;
  m2denv.dcoef.clear();
}

/****************************************************************************/
/*!
 **  Extract a vector containing the constraints
 **
 ** \param[in]  m2denv      M2D_Environ structure
 **
 ** \param[out] ydat        Array of values at constraints samples
 **                         (Dimension: nech * nlayer)
 ** \param[out] lwork       Array of tentative values (Dimension: nlayer)
 **
 *****************************************************************************/
void MGibbs::_extractVector(M2D_Environ& m2denv,
                            VectorDouble& ydat,
                            VectorDouble& lwork)
{
  for (Id iech = 0; iech < _nechc; iech++)
  {

    /* Loop on the layers */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      lwork[ilayer] = _dbc->getZVariable(iech, ilayer);

    /* Convert from the depth to thickness */

    _converZ2Y(m2denv, 1, iech, lwork);

    /* Store in the extracted vector */

    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      YDAT(_nechc, ilayer, iech) = lwork[ilayer];
  }
}

/****************************************************************************/
/*!
 **  Print the set of constraints
 **
 ** \param[in]  title       Title
 ** \param[in]  db          Db constraints structure
 ** \param[in]  ydat        Array of gaussian values at constraints (optional)
 ** \param[in]  nprint      Maximum number of lines to be printed
 **
 ** \remarks This function tends to produce verbose outputs.
 ** \remarks This is the reason why it has been conditioned to print only
 ** \remarks the values of the first samples. This is controled by the
 ** \remarks internal parameter 'nprint' which can be ruled by keypair:
 ** \remarks     set.keypair("Print_Data",0)
 ** \remarks The default number of samples i s 0 (no printout)
 **
 *****************************************************************************/
void MGibbs::_printDbConstraints(const char* title,
                                 Db* db,
                                 const VectorDouble& ydat,
                                 Id nprint) const
{
  if (!_verbose || nprint == 0) return;

  // Printout

  mestitle(1, title);
  Id nech = db->getNSample();
  if (nprint > 0) nech = MIN(nech, nprint);
  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      double lower = db->getLocVariable(ELoc::L, iech, ilayer);
      double upper = db->getLocVariable(ELoc::U, iech, ilayer);
      double value = db->getZVariable(iech, ilayer);
      double drift = db->getLocVariable(ELoc::F, iech, ilayer);
      double vgaus = (!ydat.empty()) ? YDAT(nech, ilayer, iech) : TEST;
      _printConstraintsPerPoint(ilayer, iech, value, drift, vgaus, lower, upper);
    }
  }
}

bool MGibbs::_checkArguments()
{
  // Check Model
  if (_model == nullptr)
  {
    messerr("Model must be defined");
    return false;
  }
  for (Id icov = 0; icov < _model->getNCov(); icov++)
  {
    CovAniso* cova = _model->getCovAniso(icov);
    if (cova->getType() == ECov::MATERN)
    {
      continue;
    }
    if (cova->getType() == ECov::EXPONENTIAL)
    {
      _convertExponentialToMatern(cova);
      continue;
    }
    if (cova->getType() == ECov::NUGGET)
    {
      if (_model->getSill(icov, 0, 0) > 0)
        _setFilterNugget(_model->getCovAnisoList()->isFiltered(icov));
    }
    else
    {
      messerr("SPDE Model can only support:");
      messerr("- Matern basic structures");
      messerr("- Exponential basic structures");
      messerr("- A complementary Nugget Effect");
      return false;
    }
  }
  if (_getNcovaWithoutNugget() <= 0)
  {
    messerr("The SPDE procedure requires at least one Bessel structure");
    return false;
  }

  _ndim = static_cast<Id>(_model->getNDim());
  if (_ndim != 2)
  {
    messerr("This application is restricted to the 2-D case (ndim=%d)", _ndim);
    return false;
  }
  _nvar = _model->getNVar();
  if (_nvar != 1)
  {
    messerr("This function should be called in the case of a single Model");
    messerr("In your case: %d\n", _nvar);
    return false;
  }

  // Check Input Db
  if (_dbin == nullptr)
  {
    messerr("Input Db must be defined");
    return false;
  }
  if (_dbin->getNDim() != _ndim)
  {
    messerr("Model(%d) and input Db(%d) must have same space dimension",
            _ndim, _dbin->getNDim());
    return false;
  }
  _nechin = _dbin->getNSample();

  // Check Output Db
  if (_dbout == nullptr)
  {
    messerr("Output Db must be defined");
    return false;
  }
  if (!_dbout->isGrid())
  {
    messerr("This application is restricted to a Grid output Db");
    return false;
  }
  if (_dbout->getNDim() != _ndim)
  {
    messerr("Model(%d) and output Db(%d) must have same space dimension",
            _ndim, _dbout->getNDim());
    return false;
  }
  _nechout = _dbout->getNSample();

  return true;
}

bool MGibbs::_checkValidOption()
{
  if (_nlayer <= 0)
  {
    messerr("This application requires the Number of Layers to be positive");
    return false;
  }
  if (_dbin->getNInterval() > _nlayer)
  {
    messerr("This application requires Lower and Upper variables (nint=%d)",
            _dbin->getNInterval());
    messerr("to be defined in the Input Db for each layer (nlayer=%d)", _nlayer);
    return false;
  }
  if (_flagED && _nlayer > _dbout->getNLoc(ELoc::F))
  {
    messerr("External Drifts are used for Drift definition");
    messerr("- Count of F-variables (%d) must match Count of layers (%d)",
            _dbout->getNLoc(ELoc::F), _nlayer);
    return false;
  }
  if (_nbsimu <= 0)
  {
    if (!_flagDrift)
    {
      messerr("When 'nbsimu=0', the option 'flag.drift' is set to TRUE");
      messerr("Then the Optimal Drift is calculated only");
    }
    _flagDrift = true;
  }
  if (_icolPinch >= 0)
  {
    if (!_checkValidPinchout()) return false;
  }

  return true;
}

Id MGibbs::_perform(M2D_Environ& m2denv, Id iatt_out, VectorDouble& gwork)
{
  SPDE_Matelem& Matelem = _getCurrentMatelem(0);
  Id nvertex            = _getNvertex(0);
  MatrixSparse* Bproj   = nullptr;
  QChol* QC             = nullptr;
  double ysigma         = 0.;

  /* Core allocation */

  VectorDouble lwork(_nlayer);

  VectorDouble ydat(_nechc * _nlayer, 0);
  VectorDouble ydatM(_nechc * _nlayer, 0);
  VectorDouble ygrid(nvertex * _nlayer, 0);

  VectorDouble ydat_loc(_nechc);
  VectorDouble ydatM_loc(_nechc);
  VectorDouble ygrid_loc(nvertex);

  VectorDouble rhs(nvertex, 0);
  VectorDouble zkrig(nvertex, 0);
  VectorDouble vwork(nvertex, 0);

  /* Extract the vector of current data */

  if (_verbose) message("\n==> Extracting the Initial Values at Wells\n");
  _extractVector(m2denv, ydat, lwork);
  _printDbConstraints("List of Constraining Data at Wells", _dbc, ydat);
  if (DEBUG) VH::dumpStats("Data vector (initial)", ydat);

  /* Print environment just before entering in iterative process */

  if (_verbose)
    _printEnviron("Environment before Simulations", m2denv);

  /* Loop on the simulations */

  for (Id isimu = 0; isimu < _nbsimu; isimu++)
  {
    message("Simulation #%d/%d\n", isimu + 1, _nbsimu);

    /* Loop on Gibbs iterations */

    double nugget = _model->getTotalSill(0, 0);
    for (Id iter = 0; iter < _niter; iter++)
    {
      if (_verbose) message(">>>> Iteration #%d/%d\n", iter + 1, _niter);

      // Update the Cholesky matrix

      if (iter == 0 && isimu == 0)
      {
        nugget /= 100.;
        ysigma = sqrt(nugget);
        QC     = _deriveQC(nugget, QC, Matelem);
        if (QC == nullptr) return 1;
      }

      // Perform the conditional simulation at meshing vartices

      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {
        VH::extractInPlace(ydat, ydat_loc, ilayer * _nechc);
        VH::extractInPlace(ydatM, ydatM_loc, ilayer * _nechc);
        VH::extractInPlace(ygrid, ygrid_loc, ilayer * nvertex);

        // Non-conditional simulation

        for (Id ip = 0; ip < nvertex; ip++)
          vwork[ip] = law_gaussian();
        _cholSimulate(QC, ygrid_loc, vwork);
        VH::multiplyConstant(ygrid_loc, ysigma);

        // Conditional simulation

        Matelem.Aproj->prodMatVecInPlace(ydat_loc, rhs, true);
        zkrig.fill(0.);
        _cholInvert(QC, zkrig, rhs);
        VH::addInPlace(ygrid_loc, zkrig);

        // Project the Simulation from the vertices onto the Data

        Matelem.Aproj->prodVecMatInPlace(ygrid_loc, ydatM_loc, true);

        // TODO: addition to be checked
        VH::mergeInPlace(ygrid_loc, ygrid, ilayer * nvertex);
        VH::mergeInPlace(ydatM_loc, ydatM, ilayer * _nechc);
      }

      // Perform a Gibbs iteration on the constraints

      if (_gibbs(m2denv, iter, ysigma, ydatM, ydat, lwork)) return 1;
      if (DEBUG) VH::dumpStats("Data vector after Gibbs", ydat);
    }

    /* Check that the Constraints on the Wells are honored */

    if (!_checkGibbsData("Checking Constraints at Wells", m2denv, ydat, lwork)) return 1;

    /* Store the conditional simulation on the grid */

    Bproj = dynamic_cast<MatrixSparse*>(Matelem.amesh->createProjMatrix(_dbout, -1, false));
    if (Bproj == nullptr) return 1;

    /* Project from vertices to grid nodes */

    VectorDouble gout(_nechout);
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      VH::extractInPlace(ygrid, ygrid_loc, ilayer * nvertex);
      Bproj->prodVecMatInPlace(ygrid_loc, gout, false);
      for (Id igrid = 0; igrid < _nechout; igrid++)
        GWORK(ilayer, igrid) = gout[igrid];
    }
    if (DEBUG) VH::dumpStats("Depth on grid (Gaussian)", gwork);

    /* Convert from Gaussian to Depth */

    for (Id igrid = 0; igrid < _nechout; igrid++)
    {
      if (!_dbout->isActive(igrid)) continue;
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
        lwork[ilayer] = GWORK(ilayer, igrid);
      _convertY2Z(m2denv, _dbout, 2, igrid, lwork);
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
        GWORK(ilayer, igrid) = lwork[ilayer];
    }

    if (DEBUG) VH::dumpStats("Depth on grid", gwork);
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      _dbout->setColumnByUIDOldStyle(&GWORK(ilayer, 0),
                                     iatt_out + isimu * _nlayer + ilayer);
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Perform Gibbs on a multilayer setup
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id MGibbs::run()
{
  MatrixSparse* Bproj = nullptr;
  M2D_Environ m2denv;
  SPDE_Option s_option;

  /* Initializations */

  Id error       = 1;
  Id iatt_f      = -1;
  Id iptr_ce     = -1;
  Id iptr_cstd   = -1;
  Id number_hard = 0;
  QChol* Qc      = nullptr;

  law_set_random_seed(_seed);
  VectorDouble gwork(_nechout * _nlayer);
  VectorDouble lwork(_nlayer);

  /* Preliminary checks */
  if (!_checkValidOption()) return 1;

  /* Prepare the M2D_Environ structure */

  double vartot = _model->getTotalSill(0, 0);
  _manageM2denv(sqrt(vartot), m2denv);

  /* Preparing the variables in 'dbout' */

  Id nfois    = (_flagDrift) ? 1 : _nbsimu;
  Id iatt_out = _dbout->addColumnsByConstant(_nlayer * nfois, TEST);
  if (iatt_out < 0) goto label_end;

  /* Global statistics on Raw elevations */

  _statsInit(m2denv);

  /* Manage the Drift: define External Drift on input and output Db */

  if (_verbose)
    message("\n==> Migrating Drift Information from Grid to Wells\n");
  if (_manageDrift(m2denv, &iatt_f)) goto label_end;
  _printDbConstraints("List of Initial Constraining Data", _dbin, VectorDouble());

  /* Constitute the new Db containing all the inequality constraints */
  /* whether they belong to 'dbin' or to 'dbout' */

  if (_verbose)
    message("\n==> Creating a Temporary Data Base with all constraints\n");
  if (_createConstraints(&number_hard)) goto label_end;
  if (_verbose) _dbc->display();

  /* Check SPDE environment */
  // At the first call, only one variable is Z_locatorized in order to
  // let the checks be performed on a mono-variate case (as all variables
  // will share the same Q matrix)
  // Then the environment is set to the multivariate case
  if (_verbose) message("\n==> Checking SPDE Environment\n");
  _defineLocators(_dbc, 1);
  if (!_checkSPDE()) goto label_end;
  _defineLocators(_dbc, _nlayer);

  /* Define initial values at constraints and set in Db */

  if (_verbose) message("\n==> Creating Initial Value within bounds at Wells\n");
  if (_initializeElevations(m2denv, lwork)) goto label_end;

  /* Global statistics on Centered Elevations */

  _statsUpdate(m2denv);

  /* Fitting the coefficients of the drift (external or not) */

  if (_verbose) message("\n==> Fitting the optimal Drift(s)\n");
  if (_driftFitting(m2denv, number_hard)) goto label_end;

  /* Save the drift only (optional) */

  if (_flagDrift)
  {
    _driftSave(m2denv, gwork);
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
    {
      _dbout->setColumnByUIDOldStyle(&GWORK(ilayer, 0), iatt_out + ilayer);
      (void)gslSPrintf(string_encode, "Drift%d", ilayer + 1);
      _dbout->setNameByUID(iatt_out + ilayer, string_encode);
    }
    error = 0;
    goto label_end;
  }

  /**********************************************************************/
  /* From now on, the information is stored as drift increment          */
  /**********************************************************************/

  /* Manage Drift (corrected by pinch-out) is stored in 'dbc' and 'dbout' */

  if (_verbose) message("\n==> Transforming Drift information as Thickness\n");
  if (_manageDriftIncrement(m2denv, 1)) goto label_end;

  /* Prepare all material */

  if (_verbose) message("\n==> Preparing SPDE\n");
  s_option.options = std::vector<SPDE_SS_Option>();

  if (_prepar()) goto label_end;

  if (_perform(m2denv, iatt_out, gwork)) goto label_end;

  // Renaming the simulation outcomes

  for (Id isimu = 0, ecr = 0; isimu < _nbsimu; isimu++)
  {
    for (Id ilayer = 0; ilayer < _nlayer; ilayer++, ecr++)
    {
      (void)gslSPrintf(string_encode, "Layer-%d_Simu-%d", ilayer + 1, isimu + 1);
      _dbout->setNameByUID(iatt_out + ecr, string_encode);
    }
  }

  /* Convert the simulations into the mean and variance */

  if (_flagCE || _flagCStd)
  {
    // Modify the locator to ELoc::GAUSFAC before grouping to CE estimation

    _dbout->setLocatorsByUID(_nbsimu * _nlayer, iatt_out, ELoc::GAUSFAC, 0);

    if (db_simulations_to_ce(_dbout, ELoc::GAUSFAC, _nbsimu, _nlayer, &iptr_ce,
                             &iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    if (!_flagCE)
    {
      _dbout->deleteColumnsByUIDRange(iptr_ce, _nlayer);
      iptr_ce = -1;
    }
    if (!_flagCStd)
    {
      _dbout->deleteColumnsByUIDRange(iptr_cstd, _nlayer);
      iptr_cstd = -1;
    }
    _dbout->deleteColumnsByLocator(ELoc::GAUSFAC);

    // Renaming the resulting variables

    if (iptr_ce >= 0)
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {
        (void)gslSPrintf(string_encode, "Layer-%d_CE", ilayer + 1);
        _dbout->setNameByUID(iptr_ce + ilayer, string_encode);
      }
    if (iptr_cstd >= 0)
      for (Id ilayer = 0; ilayer < _nlayer; ilayer++)
      {
        (void)gslSPrintf(string_encode, "Layer-%d_CStd", ilayer + 1);
        _dbout->setNameByUID(iptr_cstd + ilayer, string_encode);
      }
  }

  /* Set the error code */
  error = 0;

label_end:
  (void)_manageDriftIncrement(m2denv, -1);
  _manageQchol(-1, Qc);
  delete Bproj;
  if (iatt_f >= 0) _dbin->deleteColumnsByUIDRange(iatt_f, _nlayer);
  if (error && iatt_out >= 0)
    _dbout->deleteColumnsByUIDRange(iatt_out, _nlayer);
  return (error);
}

void MGibbs::setOptions(Id nlayer,
                        Id niter,
                        Id seed,
                        Id nbsimu,
                        Id icol_pinch,
                        bool flag_ed,
                        bool flag_drift,
                        bool flag_ce,
                        bool flag_cstd,
                        bool verbose)
{
  _nlayer    = nlayer;
  _niter     = niter;
  _seed      = seed;
  _nbsimu    = nbsimu;
  _icolPinch = icol_pinch;
  _flagED    = flag_ed;
  _flagDrift = flag_drift;
  _flagCE    = flag_ce;
  _flagCStd  = flag_cstd;
  _verbose   = verbose;
}

/****************************************************************************/
/*!
 **  Perform Gibbs on a multilayer setup
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 ** \param[in]  model       Model structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  niter       Number of iterations
 ** \param[in]  nbsimu      Number of simulaations
 ** \param[in]  icol_pinch  Address of the variable containing the pinchout
 ** \param[in]  flag_ed     1 if External Drit is used
 ** \param[in]  flag_drift  True to return the drift only
 **                         False the simulations
 ** \param[in]  flag_ce     True if the conditional expectation
 **                         should be returned instead of simulations
 ** \param[in]  flag_cstd   True if the conditional standard deviation
 **                         should be returned instead of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  verbose     Verbose option
 **
 ** \remarks In 'dbin':
 ** \remarks - the lower and upper bounds must be defined for each datum
 ** \remarks   (set to the locator ELoc::L and ELoc::U
 ** \remarks In 'dbout':
 ** \remarks - the trend (if flag_ed is 1) must be defined and set to
 ** \remarks   the locator ELoc::F
 ** \remarks When defined, the pinchout should be defined as a grid variable
 ** \remarks with values ranging between 0 and 1 (FFFF are admitted).
 ** \remarks It will serve as a multiplier to the Mean thickness maps.
 **
 *****************************************************************************/
Id MLayers_spde(Db* dbin,
                Db* dbout,
                Model* model,
                Id nlayer,
                Id niter,
                Id nbsimu,
                Id icol_pinch,
                bool flag_ed,
                bool flag_drift,
                bool flag_ce,
                bool flag_cstd,
                Id seed,
                bool verbose)
{
  MGibbs mgibbs(dbin, dbout, model);
  if (!mgibbs.isValid()) return 1;

  mgibbs.setOptions(nlayer, niter, seed, nbsimu, icol_pinch,
                    flag_ed, flag_drift, flag_ce, flag_cstd, verbose);

  return mgibbs.run();
}

} // namespace gstlrn