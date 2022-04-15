/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Stats/PCA.hpp"
#include "Stats/PCAStringFormat.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixRectangular.hpp"

#include <math.h>

PCA::PCA(int nvar)
  : AStringable(),
    _nVar(0),
    _mean(),
    _sigma(),
    _eigen(),
    _Z2F(),
    _F2Z()
{
  init(nvar);
}

PCA::PCA(const Db *db, bool verbose)
    : AStringable(),
      _nVar(0),
      _mean(),
      _sigma(),
      _eigen(),
      _Z2F(),
      _F2Z()
{
  pca_compute(db, verbose);
}

PCA::PCA(Db *db, double h0, double dh, const DirParam& dirparam, bool verbose)
    : AStringable(),
      _nVar(0),
      _mean(),
      _sigma(),
      _eigen(),
      _Z2F(),
      _F2Z()
{
  maf_compute(db, h0, dh, dirparam, verbose);
}

PCA::PCA(const PCA &m)
    : AStringable(m),
      _nVar(m._nVar),
      _mean(m._mean),
      _sigma(m._sigma),
      _eigen(m._eigen),
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
    _eigen = m._eigen;
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
  _eigen.resize(nvar,0);
  _Z2F.resize(nvar * nvar,0);
  _F2Z.resize(nvar * nvar,0);
}

int PCA::_calculateEigen(VectorDouble& c0)
{
  int nvar = _nVar;
  _eigen.resize (nvar,0);
  _Z2F.resize(nvar * nvar,0);
  _F2Z.resize(nvar * nvar,0);

  /* Eigen decomposition */

  if (matrix_eigen(c0.data(), nvar, _eigen.data(), _Z2F.data())) return (1);
  matrix_transpose(nvar,  nvar, _Z2F.data(), _F2Z.data());

  return(0);
}

String PCA::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const PCAStringFormat* pcafmt = dynamic_cast<const PCAStringFormat*>(strfmt);
  PCAStringFormat dsf;
  if (pcafmt != nullptr) dsf = *pcafmt;

  sstr << toTitle(1, "PCA Contents");

  if (dsf.getflagCenter())
  {
    sstr << toMatrix("Means", VectorString(), VectorString(), true, 1, _nVar,
                    _mean);
    sstr << toMatrix("Standard deviations", VectorString(), VectorString(), true,
                    1, _nVar, _sigma);
  }
  if (dsf.getflagStats())
    sstr << toMatrix("Eigen Values", VectorString(), VectorString(), true, 1,
                    _nVar, _eigen);

  sstr << toMatrix("Matrix M to transform standardized Variables Z into Factors Y",
                   VectorString(), VectorString(), true,
                   _nVar, _nVar, _Z2F);
  sstr << "Y = Z * M (columns  = eigen vectors)" << std::endl;
  sstr << toMatrix("Matrix t(M) to back-transform Factors Y into standardized Variables Z",
                   VectorString(), VectorString(), true,
                   _nVar, _nVar, _F2Z);
  sstr << "Z = Y * t(M) (rows  = eigen vectors)" << std::endl;

  return sstr.str();
}

int PCA::pca_compute(const Db *db, bool verbose)
{

  /* Initializations */

  if (db == nullptr)
  {
    messerr("You must define the 'Db'");
    return 1;
  }
  int nvar = db->getVariableNumber();
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
  if (_pcaCalculate(db, isoFlag, verbose)) return 1;
  return 0;
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
  int nvar = db->getVariableNumber();
  if (nvar != _nVar)
  {
    messerr("The number of Z variables (%d) does not match the number of variables in PCA (%d)",
            nvar,_nVar);
    return 1;
  }

  // Core allocation

  VectorDouble data(nvar);
  VectorDouble mean(nvar);
  VectorDouble sigma(nvar);
  VectorBool isoFlag = _getVectorIsotopic(db);

  /* Allocate new variables */

  int iptr = db->addColumnsByConstant(nvar, TEST);
  if (iptr < 0) return 1;

  // Optional title

  if (verbose) mestitle(0,"Transform from Z to Factors");

  /* Normalization (optional) */

  if (_normalization(db, isoFlag, mean, sigma, verbose)) return 1;

  /* Perform the normalization */

  _pcaZ2F(false, iptr, db, isoFlag, mean, sigma);

  /* Optional printout */

  if (verbose)
  {
    VectorInt cols(nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      cols[ivar] = iptr + ivar;
    db_stats_print(db, cols, VectorString(), 1, 1, "Statistics on Factors","Factor");
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
  int nvar = db->getVariableNumber();
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
  _pcaF2Z(iptr, db, isoFlag);

  /* Optional printout */

  if (verbose)
  {
    VectorInt cols(nvar);
    for (int ivar = 0; ivar < nvar; ivar++) cols[ivar] = iptr + ivar;
    db_stats_print(db, cols, VectorString(), 1, 1, "Statistics on Variables", "Variable");
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, -1, db, iptr);
  return 0;
}

/****************************************************************************/
/*!
 **  Internal function to calculate MAF
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  isoFlag    Vector of active samples
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int PCA::_pcaCalculate(const Db *db,
                       const VectorBool& isoFlag,
                       bool verbose)
{
  int nvar = getNVar();
  VectorDouble mean(nvar);
  VectorDouble sigma(nvar);
  VectorDouble c0(nvar * nvar);
  if (_normalization(db, isoFlag, mean, sigma, verbose)) return 1;
  _covariance0(db, isoFlag, mean, sigma, c0, verbose);
  setMean(mean);
  setSigma(sigma);
  if (_calculateEigen(c0)) return 1;
  return 0;
}

/****************************************************************************/
/*!
 **  Fill the mean and variance arrays
 **

 ** \param[in] db          Db descriptor
 ** \param[in] isoFlag     Vector of active samples
 ** \param[in] verbose     Verbose flag
 **
 ** \param[out] mean       Array of means
 ** \param[out] sigma      Array of standard deviations
 **
 *****************************************************************************/
int PCA::_normalization(const Db *db,
                        const VectorBool& isoFlag,
                        VectorDouble& mean,
                        VectorDouble& sigma,
                        bool verbose)
{
  int niso, nvar, nech;

  // Initializations

  niso = 0;
  nvar = db->getVariableNumber();
  nech = db->getSampleNumber();
  VectorDouble data(nvar);

  for (int ivar = 0; ivar < nvar; ivar++)
    mean[ivar] = sigma[ivar] = 0.;

  /* Calculate the statistics */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!isoFlag[iech]) continue;
    _loadData(db, iech, data);

    niso++;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      mean[ivar] += data[ivar];
      sigma[ivar] += data[ivar] * data[ivar];
    }
  }

  /* Normalization */

  if (niso > 0)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      mean[ivar] /= niso;
      sigma[ivar] = (sigma[ivar] / niso - mean[ivar] * mean[ivar]);
      sigma[ivar] = (sigma[ivar] > 0) ? sqrt(sigma[ivar]) :
                                        0.;
      if (sigma[ivar] <= 0.)
      {
        messerr("Error: Variable (%d) is constant", ivar + 1);
        return 1;
      }
    }
  }

  if (verbose)
  {
    message("Number of variables         = %d\n", nvar);
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotropic samples = %d\n", niso);
    print_matrix("Mean", 0, 1, 1, nvar, NULL, mean.data());
    print_matrix("St. Dev.", 0, 1, 1, nvar, NULL, sigma.data());
    message("\n");
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Internal PCA covariance calculation
 **
 ** \return The array of covariance (or NULL)
 **
 ** \param[in]  db          Db descriptor
 ** \param[in]  isoFlag     Vector of active samples
 ** \param[in]  mean        Array containing the mean
 ** \param[in]  sigma       Array containing the standard deviation
 ** \param[in]  verbose     Verbose flag
 **
 ** \param[out] c0          Vector of covariances at distance 0
 **
 *****************************************************************************/
void PCA::_covariance0(const Db *db,
                       const VectorBool& isoFlag,
                       const VectorDouble& mean,
                       const VectorDouble& sigma,
                       VectorDouble& c0,
                       bool verbose)
{
  int nvar = db->getVariableNumber();
  int nech = db->getSampleNumber();
  int niso = 0;
  VectorDouble data1(nvar);
  for (int i = 0; i < nvar * nvar; i++) c0[i] = 0.;

  /* Calculate the variance-covariance matrix at distance 0 */

  niso = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);
    _center(data1, mean, sigma);

    niso++;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        c0[ivar * nvar + jvar] += data1[ivar] * data1[jvar];
  }

  /* Normalization */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      c0[ivar * nvar + jvar] /= niso;

  /* Printout of the covariance matrix (optional) */

  if (verbose)
    print_matrix("Variance-Covariance matrix for distance 0", 0, 1, nvar, nvar,
                 NULL,c0.data());
}

/****************************************************************************/
/*!
 **  Normalize the isotropic array of values
 **
 ** \param[in,out] data      Array of information
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 *****************************************************************************/
void PCA::_center(VectorDouble& data,
                  const VectorDouble &mean,
                  const VectorDouble &sigma)
{
  int ivar;
  int nvar = (int) mean.size();

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (sigma[ivar] <= 0.) continue;
    data[ivar] = (data[ivar] - mean[ivar]) / sigma[ivar];
  }
}

/****************************************************************************/
/*!
 **  Un-normalize the isotropic array of values
 **
 ** \param[in,out] data      Array of information
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 *****************************************************************************/
void PCA::_uncenter(VectorDouble& data,
                    const VectorDouble &mean,
                    const VectorDouble &sigma)
{
  int ivar;
  int nvar = (int) mean.size();

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (sigma[ivar] <= 0.) continue;
    data[ivar] = data[ivar] * sigma[ivar] + mean[ivar];
  }
}

/****************************************************************************/
/*!
 **  Procedure for transforming the variables into factors using PCA
 **
 ** \param[in]  flag_norm    True if the variables must be rescaled beforehand
 ** \param[in]  iptr         Pointer for storing the result in db
 ** \param[in]  db           Db descriptor
 ** \param[in]  isoFlag      Vector of active samples
 ** \param[in]  mean         Array containing the mean
 ** \param[in]  sigma        Array containing the standard deviation
 **
 ** \remarks The standardization statistics are passed as arguments
 ** \remarks The ones stored in PCA structure are not used.
 **
 *****************************************************************************/
void PCA::_pcaZ2F(bool flag_norm,
                  int iptr,
                  Db *db,
                  const VectorBool isoFlag,
                  const VectorDouble& mean,
                  const VectorDouble& sigma)
{
  int nvar = db->getVariableNumber();
  int nech = db->getSampleNumber();
  VectorDouble data1(nvar);
  VectorDouble data2(nvar);

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (! isoFlag[iech]) continue;
    _loadData(db, iech, data1);
    _center(data1, mean, sigma);

    /* Loop on the factors */

    for (int ifac = 0; ifac < nvar; ifac++)
    {
      double value = 0.;
      for (int ivar = 0; ivar < nvar; ivar++)
        value += getZ2F(ifac, ivar) * data1[ivar];
      if (flag_norm) value /= sqrt(_eigen[ifac]);
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
 **
 *****************************************************************************/
void PCA::_pcaF2Z(int iptr, Db *db, const VectorBool& isoFlag)
{
  int nvar = db->getVariableNumber();
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
    _uncenter(data2, getMean(), getSigma());

    // Storage

    for (int ivar = 0; ivar < nvar; ivar++)
      db->setArray(iech, ivar + iptr, data2[ivar]);
  }
}

/****************************************************************************/
/*!
 **  Evaluate the MAF on irregular data
 **
 ** \return  PCA structure (or nullptr)
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  dirparam   DirParam structure
 ** \param[in]  h0         Reference distance
 ** \param[in]  dh         Tolerance on distance
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int PCA::maf_compute(Db *db,
                     double h0,
                     double dh,
                     const DirParam& dirparam,
                     bool verbose)
{

  /* Initializations */

  if (db == nullptr)
  {
    messerr("You must define 'Db' beforehand");
    return 1;
  }
  int nvar = db->getVariableNumber();
  if (nvar <= 0)
  {
    messerr("You must define 'Db' with some Z-variables");
    return 1;
  }
  init(nvar);
  PCA pca2(nvar);

  /* Allocate new variables */

  VectorInt uidZ = db->getUIDsByLocator(ELoc::Z);
  int iptr = db->addColumnsByConstant(nvar, TEST);
  if (iptr < 0) return 1;

  // Optional title

  if (verbose) mestitle(0,"MAF computation");

  /* Core allocation */

  VectorBool isoFlag = _getVectorIsotopic(db);
  VectorDouble identity(nvar * nvar);
  VectorDouble ch(nvar * nvar);
  for (int i = 0; i < nvar * nvar; i++)
    identity[i] = 0.;

  /* Calculate the first PCA (centered and normalized) */

  if (_pcaCalculate(db, isoFlag, verbose)) return 1;
  if (verbose) display();

  /* Rotate the initial data in the PCA system */

  _pcaZ2F(true, iptr, db, isoFlag, getMean(), getSigma());
  db->setLocatorsByUID(nvar, iptr, ELoc::Z);

  /* Calculate the variance-covariance matrix at distance [h0-dh,h0+dh] */

  if (_covarianceh(db, h0, dh, dirparam, isoFlag, ch, verbose)) return 1;
  if (pca2._calculateEigen(ch)) return 1;

  /* Rotate the initial data in the second PCA system */

  for (int ivar = 0; ivar < nvar; ivar++)
    identity[ivar * nvar + ivar] = 1. / sqrt(getEigen(ivar));
  VectorDouble pcaz2f  = getZ2F();
  VectorDouble pca2z2f = pca2.getZ2F();
  VectorDouble pcaf2z  = getF2Z();
  matrix_product(nvar, nvar, nvar, pcaz2f.data(), identity.data(), pcaz2f.data());
  matrix_product(nvar, nvar, nvar, pcaz2f.data(), pca2z2f.data(), pca2z2f.data());
  setPcaZ2F(pcaz2f);
  pca2.setPcaZ2F(pca2z2f);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int kvar = nvar - ivar - 1;
    setEigen(kvar, pca2.getEigen(ivar));
    for (int jvar = 0; jvar < nvar; jvar++)
      setPcaZ2F(kvar, jvar, pca2.getZ2F(ivar, jvar));
  }
  if (matrix_invert_copy(getZ2F().data(), nvar, pcaf2z.data())) return 1;
  setPcaF2Z(pcaf2z);

  (void) db_attribute_del_mult(db, iptr, nvar);
  db->setLocatorsByUID(uidZ, ELoc::Z);
  return 0;
}

/****************************************************************************/
/*!
 **  Internal PCA translated covariance calculation
 **
 **  Error returned case
 **
 ** \param[in]  db          Db descriptor
 ** \param[in]  h0          Reference distance
 ** \param[in]  dh          Tolerance on distance
 ** \param[in]  dirparam    DirParam structure
 ** \param[in]  isoFlag     Vector of active samples
 ** \param[in]  verbose     Verbose flag
 **
 ** \param[out] ch          Vector of covariances at distance h
 **
 *****************************************************************************/
int PCA::_covarianceh(Db *db,
                      double h0,
                      double dh,
                      const DirParam& dirparam,
                      const VectorBool& isoFlag,
                      VectorDouble& ch,
                      bool verbose)
{
  double ps;

  int nech = db->getSampleNumber();
  int nvar = db->getVariableNumber();
  int npairs = 0;
  double psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());
  VectorDouble data1(nvar);
  VectorDouble data2(nvar);

  // Initialization

  for (int i = 0; i < nvar * nvar; i++)
    ch[i] = 0.;

  /* Core allocation */

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
      if (dist < h0 - dh || dist > h0 + dh) continue;
      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;

      /* Update the variance-covariance matrix at distance h */

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          double di = data1[ivar] - data2[ivar];
          double dj = data1[jvar] - data2[jvar];
          ch[ivar * nvar + jvar] += di * dj;
        }
      npairs++;
    }
  }

  /* Normation */

  if (npairs <= 1)
  {
    messerr("Number of pairs of isotopic samples is smaller than 2");
    return 1;
  }
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      ch[ivar * nvar + jvar] /= npairs;

  /* Verbose printout */

  if (verbose)
  {
    message("\n");
    dirparam.display();
    message("Reference Distance          = %lf\n", h0);
    message("Tolerance on Distance       = %lf\n", dh);
    message("Number of samples in the Db = %d\n", nech);
    message("Number of isotopic pairs    = %d\n", npairs);
    message("\n");
    print_matrix("Variance-Covariance matrix for distance h", 0, 1, nvar, nvar,
                 NULL,ch.data());
  }
  return 0;
}

VectorBool PCA::_getVectorIsotopic(const Db* db)
{
  int nvar = db->getVariableNumber();
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
  int nvar = (int) db->getVariableNumber();
  for (int ivar = 0; ivar < nvar; ivar++)
    data[ivar] = db->getVariable(iech, ivar);
}

VectorDouble PCA::mafOfIndex() const
{
  // Calculate the probability of each interval
  VectorDouble w = _mean;
  int ncut = (int) _mean.size();
  w.push_back(1 - ut_vector_cumul(_mean));
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
  matrix_product(nclass, ncut, ncut, i_norm_val.getValues().data(), _Z2F.data(), local.data());

  VectorDouble maf_index = ut_vector_concatenate(ut_vector_double(nclass, 1.), local);

  return maf_index;
}

