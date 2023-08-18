/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Basic/VectorHelper.hpp"
#include "Stats/PCA.hpp"
#include "Stats/PCAStringFormat.hpp"
#include "Stats/Classical.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

#include <math.h>

PCA::PCA(int nvar)
  : AStringable(),
    _nVar(0),
    _mean(),
    _sigma(),
    _eigval(),
    _eigvec(),
    _c0(),
    _gh(),
    _Z2F(),
    _F2Z()
{
  init(nvar);
}

PCA::PCA(const PCA &m)
    : AStringable(m),
      _nVar(m._nVar),
      _mean(m._mean),
      _sigma(m._sigma),
      _eigval(m._eigval),
      _eigvec(m._eigvec),
      _c0(m._c0),
      _gh(m._gh),
      _Z2F(m._Z2F),
      _F2Z(m._F2Z)
{

}

PCA& PCA::operator=(const PCA &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nVar = m._nVar;
    _mean = m._mean;
    _sigma = m._sigma;
    _eigval = m._eigval;
    _eigvec = m._eigvec;
    _c0 = m._c0;
    _gh = m._gh;
    _Z2F = m._Z2F;
    _F2Z = m._F2Z;
  }
  return *this;
}

PCA::~PCA()
{

}

void PCA::init(int nvar)
{
  _nVar = nvar;
  _mean.resize(nvar);
  _sigma.resize(nvar,0);
  _eigval.resize(nvar,0);
  _eigvec.resize(nvar * nvar,0);
  _c0.resize(nvar * nvar,0);
  _gh.resize(nvar * nvar,0);
  _Z2F.resize(nvar * nvar,0);
  _F2Z.resize(nvar * nvar,0);
}

void PCA::_pcaFunctions(bool verbose)
{
  int nvar = _nVar;

  // Transpose for getting the F2Z function from Z2F

  matrix_transpose(nvar,  nvar, _eigvec.data(), _F2Z.data());

  // Construct Z2F

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      setZ2F(ivar, jvar, getEigVec(ivar, jvar) / sqrt(_eigval[ivar]));

  // Construct F2Z

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      setF2Z(ivar, jvar, getF2Z(ivar, jvar) * sqrt(_eigval[jvar]));

  // Printout of the transition matrices (optional)

  if (verbose)
  {
    print_matrix("PCA Z->F", 0, 1, nvar, nvar, NULL,_Z2F.data());
    print_matrix("PCA F->Z", 0, 1, nvar, nvar, NULL,_F2Z.data());
  }
}

void PCA::_mafFunctions(bool verbose)
{
  int nvar = _nVar;

  // Construct Z2F

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      setZ2F(ivar, jvar, getEigVec(ivar, jvar));

  // Construct F2Z

  MatrixSquareGeneral A(nvar);
  A.setValues(_Z2F);
  A.invert();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      setF2Z(ivar, jvar, A(ivar, jvar));

  // Printout of the transition matrices (optional)

  if (verbose)
  {
    print_matrix("MAF Z->F", 0, 1, nvar, nvar, NULL,_Z2F.data());
    print_matrix("MAF F->Z", 0, 1, nvar, nvar, NULL,_F2Z.data());
  }
}

int PCA::_calculateEigen(bool verbose)
{
  int nvar = _nVar;

  // Eigen decomposition

  if (matrix_eigen(_c0.data(), nvar, _eigval.data(), _eigvec.data())) return 1;

  // Printout of the eigen results (optional)

  if (verbose)
  {
    print_matrix("Eigen values", 0, 1, 1, nvar, NULL, _eigval.data());
    print_matrix("Eigen Vectors", 0, 1, nvar, nvar, NULL,_eigvec.data());
  }
  return 0;
}

int PCA::_calculateGEigen(bool verbose)
{
  int nvar = _nVar;

  // Eigen decomposition

  if (matrix_geigen(_gh.data(), _c0.data(), nvar, _eigval.data(), _eigvec.data())) return 1;

  // Printout of the eigen results (optional)

  if (verbose)
  {
    print_matrix("GEigen values", 0, 1, 1, nvar, NULL, _eigval.data());
    print_matrix("GEigen Vectors", 0, 1, nvar, nvar, NULL,_eigvec.data());
  }
  return 0;
}

String PCA::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const PCAStringFormat* pcafmt = dynamic_cast<const PCAStringFormat*>(strfmt);
  PCAStringFormat dsf;
  if (pcafmt != nullptr) dsf = *pcafmt;
  if (_nVar <= 0) return sstr.str();

  sstr << toTitle(1, "PCA Contents");

  if (dsf.getflagCenter())
  {
    sstr << toMatrix("Means", VectorString(), VectorString(), true, 1, _nVar,
                    _mean);
  }
  if (dsf.getflagStats())
  {
    sstr << toMatrix("Covariance Matrix", VectorString(), VectorString(), true, _nVar,
                     _nVar, _c0);
    if (_gh.size() > 0)
    {
      sstr << toMatrix("Variogram Matrix at lag h", VectorString(), VectorString(), true, _nVar,
                       _nVar, _gh);
    }

    sstr << toMatrix("Matrix MZ2F to transform standardized Variables Z into Factors F",
                     VectorString(), VectorString(), true,
                     _nVar, _nVar, _Z2F);
    sstr << "Y = (Z - m) * MZ2F)" << std::endl;
    sstr << toMatrix("Matrix MF2Z to back-transform Factors F into standardized Variables Z",
                     VectorString(), VectorString(), true,
                     _nVar, _nVar, _F2Z);
    sstr << "Z = m + Y * MF2Z" << std::endl;
  }
  return sstr.str();
}

int PCA::dbZ2F(Db* db,
               bool verbose,
               const NamingConvention& namconv)
{
  if (db == nullptr)
  {
    messerr("You must define 'Db'");
    return 1;
  }
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar != _nVar)
  {
    messerr("The number of Z variables (%d) does not match the number of variables in PCA (%d)",
            nvar,_nVar);
    return 1;
  }

  /* Allocate new variables */

  int iptr = db->addColumnsByConstant(nvar, TEST);
  if (iptr < 0) return 1;

  // Optional title

  if (verbose) mestitle(0,"Transform from Z to Factors");

  /* Rotate the factors into data in the PCA system */

  VectorBool isoFlag = _getVectorIsotopic(db);
  _pcaZ2F(iptr, db, isoFlag, _mean, _sigma);

  /* Optional printout */

  if (verbose)
  {
    VectorInt cols(nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      cols[ivar] = iptr + ivar;
    dbStatisticsPrintByUID(db, cols, {}, 1, 1, "Statistics on Factors","Factor");
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, -1, db, iptr);
  return 0;
}

int PCA::dbF2Z(Db* db,
               bool verbose,
               const NamingConvention& namconv)
{
  if (db == nullptr)
  {
    messerr("You must define 'Db'");
    return 1;
  }
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar != _nVar)
  {
    messerr("The number of Z variables (%d) does not match the number of variables in PCA (%d)",
            nvar,_nVar);
    return 1;
  }

  /* Allocate new variables */

  int iptr = db->addColumnsByConstant(nvar, TEST);
  if (iptr < 0) return 1;

  // Optional title

  if (verbose) mestitle(0,"Transform from Factors to Z");

  /* Rotate the factors into data in the PCA system */

  VectorBool isoFlag = _getVectorIsotopic(db);
  _pcaF2Z(iptr, db, isoFlag, _mean, _sigma);

  /* Optional printout */

  if (verbose)
  {
    VectorInt cols(nvar);
    for (int ivar = 0; ivar < nvar; ivar++) cols[ivar] = iptr + ivar;
    dbStatisticsPrintByUID(db, cols, {}, 1, 1, "Statistics on Variables", "Variable");
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, -1, db, iptr);
  return 0;
}

VectorDouble PCA::getVarianceRatio() const
{
  double total = VectorHelper::cumul(_eigval);
  VectorDouble eignorm = _eigval;
  VectorHelper::divideConstant(eignorm, total);
  return eignorm;
}

/****************************************************************************/
/*!
 **  Fill the mean and variance arrays
 **
 ** \param[in] db          Db descriptor
 ** \param[in] isoFlag     Vector of active samples
 ** \param[in] verbose     Verbose flag
 ** \param[in] flag_nm1    When TRUE, variance is scaled by N-1; otherwise by N
 **
 *****************************************************************************/
void PCA::_calculateNormalization(const Db *db,
                                  const VectorBool &isoFlag,
                                  bool verbose,
                                  bool flag_nm1)
{
  int niso = 0;
  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  VectorDouble data(nvar);

  for (int ivar = 0; ivar < nvar; ivar++)
    _mean[ivar] = _sigma[ivar] = 0.;

  /* Calculate the statistics */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!isoFlag[iech]) continue;
    _loadData(db, iech, data);

    niso++;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      _mean[ivar] += data[ivar];
      _sigma[ivar] += data[ivar] * data[ivar];
    }
  }

  /* Normalization */

  if (niso > 0)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      _mean[ivar] /= niso;
      _sigma[ivar] = (_sigma[ivar] / niso - _mean[ivar] * _mean[ivar]);
      if (flag_nm1) _sigma[ivar] *= (double) (niso / (niso-1.));
      _sigma[ivar] = (_sigma[ivar] > 0) ? sqrt(_sigma[ivar]) : 0.;
    }
  }

  if (verbose)
  {
    message("Number of variables         = %d\n", nvar);
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotropic samples = %d\n", niso);
    print_matrix("Mean", 0, 1, 1, nvar, NULL, _mean.data());
    print_matrix("St. Dev.", 0, 1, 1, nvar, NULL, _sigma.data());
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Calculate the covariance matrix
 **
 ** \param[in]  db          Db descriptor
 ** \param[in]  isoFlag     Vector of active samples
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  flag_nm1    When TRUE, variance is scaled by N-1; otherwise by N
 **
 *****************************************************************************/
void PCA::_covariance0(const Db *db,
                       const VectorBool &isoFlag,
                       bool verbose,
                       bool flag_nm1)
{
  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  int niso = 0;
  VectorDouble data1(nvar);

  /* Calculate the variance-covariance matrix at distance 0 */

  niso = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);
    _center(data1, _mean, _sigma, true, false);

    niso++;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        _c0[ivar * nvar + jvar] += data1[ivar] * data1[jvar];
  }

  /* Normalization */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      if (flag_nm1)
        _c0[ivar * nvar + jvar] /= (niso-1.);
      else
        _c0[ivar * nvar + jvar] /= niso;

  /* Printout of the covariance matrix (optional) */

  if (verbose)
    print_matrix("Variance-Covariance matrix for distance 0", 0, 1, nvar, nvar,
                 NULL,_c0.data());
}

/****************************************************************************/
/*!
 **  Normalize the isotropic array of values
 **
 ** \param[in,out] data      Array of information
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 ** \param[in]  flag_center  True if result must be centered by 'mean'
 ** \param[in]  flag_scale   True if result must be scaled by 'sigma'
 **
 *****************************************************************************/
void PCA::_center(VectorDouble& data,
                  const VectorDouble &mean,
                  const VectorDouble &sigma,
                  bool flag_center,
                  bool flag_scale)
{
  int nvar = (int) mean.size();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (flag_center)
      data[ivar] -= mean[ivar];
    if (flag_scale && sigma[ivar] > 0.)
      data[ivar] /= sigma[ivar];
  }
}

/****************************************************************************/
/*!
 **  Un-normalize the isotropic array of values
 **
 ** \param[in,out] data      Array of information
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 ** \param[in]  flag_center  True if result must be uncentered by 'mean'
 ** \param[in]  flag_scale   True if result must be suncaled by 'sigma'
 **
 *****************************************************************************/
void PCA::_uncenter(VectorDouble& data,
                    const VectorDouble &mean,
                    const VectorDouble &sigma,
                    bool flag_center,
                    bool flag_scale)
{
  int ivar;
  int nvar = (int) mean.size();

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (sigma[ivar] <= 0.) continue;
    if (flag_scale)
      data[ivar] *= sigma[ivar];
    if (flag_center)
      data[ivar] += mean[ivar];
  }
}

/****************************************************************************/
/*!
 **  Procedure for transforming the variables into factors using PCA
 **
 ** \param[in]  iptr         Pointer for storing the result in db
 ** \param[in]  db           Db descriptor
 ** \param[in]  isoFlag      Vector of active samples
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 *****************************************************************************/
void PCA::_pcaZ2F(int iptr,
                  Db *db,
                  const VectorBool isoFlag,
                  const VectorDouble& mean,
                  const VectorDouble& sigma)
{
  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  VectorDouble data1(nvar, 0.);
  VectorDouble data2(nvar, 0.);

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);
    _center(data1, mean, sigma, true, false);

    /* Loop on the factors */

    for (int ifac = 0; ifac < nvar; ifac++)
    {
      double value = 0.;
      for (int ivar = 0; ivar < nvar; ivar++)
        value += getZ2F(ifac, ivar) * data1[ivar];
      data2[ifac] = value;
    }

    // Storage
    for (int ifac = 0; ifac < nvar; ifac++)
      db->setArray(iech, ifac + iptr, data2[ifac]);
  }
}

/****************************************************************************/
/*!
 **  Procedure for transforming the factors into variables using PCA
 **
 ** \param[in]  iptr         Pointer to the storage
 ** \param[in]  db           Db descriptor
 ** \param[in]  isoFlag      Vector of active samples
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 *****************************************************************************/
void PCA::_pcaF2Z(int iptr,
                  Db *db,
                  const VectorBool &isoFlag,
                  const VectorDouble &mean,
                  const VectorDouble &sigma)
{
  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  VectorDouble data1(nvar);
  VectorDouble data2(nvar);

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);

    /* Loop on the factors */

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      // Calculate the projection
      double value = 0.;
      for (int ifac = 0; ifac < nvar; ifac++)
        value += getF2Z(ivar, ifac) * data1[ifac];
      data2[ivar] = value;
    }

    // De-normalize
    _uncenter(data2, mean, sigma, true, false);

    // Storage

    for (int ivar = 0; ivar < nvar; ivar++)
      db->setArray(iech, ivar + iptr, data2[ivar]);
  }
}

int PCA::pca_compute(const Db *db, bool verbose)
{

  /* Initializations */

  if (db == nullptr)
  {
    messerr("You must define the 'Db'");
    return 1;
  }
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar <= 0)
  {
    messerr("You must define 'Db' with some Z-variables");
    return 1;
  }
  init(nvar);

  // Optional title

  if (verbose) mestitle(0,"PCA computation");

  /* Calculate the PCA */

  VectorBool isoFlag = _getVectorIsotopic(db);
  _calculateNormalization(db, isoFlag, verbose, true);

  // Derive the PCA decomposition

  _covariance0(db, isoFlag, verbose, true);

  // Establish the transfer functions

  if (_calculateEigen(verbose)) return 1;

  _pcaFunctions(verbose);

  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the MAF
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  hmin       Lower bound on distance
 ** \param[in]  hmax       Upper bound on distance
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int PCA::maf_compute_interval(Db *db, double hmin, double hmax, bool verbose)
{
  return _mafCompute(db, VarioParam(), -1, -1, hmin, hmax, verbose);
}

/****************************************************************************/
/*!
 **  Evaluate the MAF
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  ilag0      Reference Lag
 ** \param[in]  idir0      Reference direction
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int PCA::maf_compute(Db *db,
                     const VarioParam& varioparam,
                     int ilag0,
                     int idir0,
                     bool verbose)
{
  return _mafCompute(db, varioparam, ilag0, idir0, TEST, TEST, verbose);
}

/****************************************************************************/
/*!
 **  Evaluate the MAF on irregular data
 **
 ** \return  Error returned code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  ilag0      Reference Lag
 ** \param[in]  idir0      Reference direction
 ** \param[in]  hmin       Minimum distance
 ** \param[in]  hmax       Maximum distance
 ** \param[in]  verbose    Verbose flag
 **
 ** \remarks This code functions in two modes:
 ** \remarks - either using varioparam, idir0 and ilag0 (if idir0>=0)
 ** \remarks - or using hmin, hmax
 **
 *****************************************************************************/
int PCA::_mafCompute(Db *db,
                     const VarioParam &varioparam,
                     int ilag0,
                     int idir0,
                     double hmin,
                     double hmax,
                     bool verbose)
{
  /* Initializations */

  if (db == nullptr)
  {
    messerr("You must define 'Db' beforehand");
    return 1;
  }
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar <= 0)
  {
    messerr("You must define 'Db' with some Z-variables");
    return 1;
  }
  init(nvar);

  // Optional title

  if (verbose) mestitle(0,"MAF computation");

  // Calculate the normalization parameters

  VectorBool isoFlag = _getVectorIsotopic(db);
  _calculateNormalization(db, isoFlag, verbose, true);

  /* Calculate the first PCA (centered and normalized) */

  _covariance0(db, isoFlag, verbose, true);

  /* Calculate the variogram matrix at distance [h0-dh,h0+dh] */

  _variogramh(db, varioparam, ilag0, idir0, hmin, hmax, isoFlag, verbose);

  // Derive the MAF decomposition

  if (_calculateGEigen(verbose)) return 1;

  // Establish the transfer functions

  _mafFunctions(verbose);

  return 0;
}

void PCA::_variogramh(Db *db,
                      const VarioParam &varioparam,
                      int ilag0,
                      int idir0,
                      double hmin,
                      double hmax,
                      const VectorBool &isoFlag,
                      bool verbose)
{
  double ps;

  // Initializations

  int nech = db->getSampleNumber();
  int nvar = db->getLocNumber(ELoc::Z);
  int npairs = 0;

  // Core allocations

  VectorDouble data1(nvar);
  VectorDouble data2(nvar);

  /* Loop on samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);

    /* Loop on the second sample */

    for (int jech = 0; jech < iech; jech++)
    {
      if (! isoFlag[jech]) continue;
      _loadData(db, jech, data2);

      /* Should the pair be retained */

      double dist = distance_intra(db, iech, jech, NULL);

      if (idir0 >= 0)
      {
        DirParam dirparam = varioparam.getDirParam(idir0);
        double psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
        double lag = dirparam.getDPas();
        double h0  = ilag0 * lag;
        if (dist < h0 - lag/2 || dist > h0 + lag/2) continue;
        if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                            (int) dirparam.getTolCode())) continue;
        if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                  dirparam.getBench(), dirparam.getCylRad(),
                                  dirparam.getCodirs(), &ps)) continue;
      }
      else
      {
        if (dist < hmin || dist > hmax) continue;
      }

      /* Update the variance-covariance matrix at distance h */

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          double di = data1[ivar] - data2[ivar];
          double dj = data1[jvar] - data2[jvar];
          _gh[ivar * nvar + jvar] += di * dj / 2.;
        }
      npairs++;
    }
  }

  /* Normation */

  if (npairs > 0)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        _gh[ivar * nvar + jvar] /= npairs;

  /* Verbose printout */

  if (verbose)
  {
    message("\n");

    if (idir0 >= 0)
    {
      DirParam dirparam = varioparam.getDirParam(idir0);
      dirparam.display();
      message("Reference Lag               = %d\n", ilag0);
    }
    else
    {
      message("Minimum distance            = %lf\n", hmin);
      message("Maximum distance            = %lf\n", hmax);
    }
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotopic pairs    = %d\n", npairs);
    message("\n");
    print_matrix("Variogram matrix for distance h", 0, 1, nvar, nvar, NULL, _gh.data());
  }
}

VectorBool PCA::_getVectorIsotopic(const Db* db)
{
  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  VectorDouble data(nvar);
  VectorBool isoFlag = VectorBool(nech);

  for (int iech = 0; iech < nech; iech++)
  {
    if (! db->isActive(iech))
      isoFlag[iech] = false;
    else
      isoFlag[iech] = db_is_isotropic(db, iech, data.data());
  }
  return isoFlag;
}

void PCA::_loadData(const Db* db, int iech, VectorDouble& data)
{
  int nvar = (int) db->getLocNumber(ELoc::Z);
  for (int ivar = 0; ivar < nvar; ivar++)
    data[ivar] = db->getLocVariable(ELoc::Z,iech, ivar);
}

VectorDouble PCA::mafOfIndex() const
{
  // Calculate the probability of each interval
  VectorDouble w = _mean;
  int ncut = (int) _mean.size();
  w.push_back(1 - VH::cumul(_mean));
  int nclass = (int) w.size();

  // Normalize the indicator of intervals
  MatrixRectangular i_norm_val(nclass, ncut);
  for (int iclass = 0 ; iclass < ncut; iclass++)
  {
    for (int jclass = 0; jclass < nclass; jclass++)
    {
      double value = (iclass == jclass) ? 1 : 0.;
      value = (value - w[iclass]) / sqrt(w[iclass] * (1. - w[iclass]));
      i_norm_val.setValue(jclass,  iclass, value);
    }
  }

  VectorDouble local(nclass * ncut);
  matrix_product_safe(nclass, ncut, ncut, i_norm_val.getValues().data(), _Z2F.data(), local.data());

  VectorDouble maf_index = VH::concatenate(VH::initVDouble(nclass, 1.), local);

  return maf_index;
}

