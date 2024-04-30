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
#include "geoslib_old_f.h"

#include "Geometry/GeometryHelper.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCustom.hpp"

#include <string.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <random>

VectorInt VectorHelper::initVInt(int nval, int value)
{
  VectorInt tab(nval, value);
  return tab;
}

VectorDouble VectorHelper::initVDouble(int nval, double value)
{
  VectorDouble tab(nval, value);
  return tab;
}

VectorVectorDouble VectorHelper::initVVDouble(int nval1, int nval2, double value)
{
  VectorVectorDouble tab(nval1, VectorDouble(nval2, value));
  return tab;
}

VectorVectorInt VectorHelper::initVVInt(int nval1, int nval2, int value)
{
  VectorVectorInt tab(nval1, VectorInt(nval2, value));
  return tab;
}

VectorInt VectorHelper::initVInt(const int* values, int number)
{
  if (values == nullptr) return VectorInt();
  VectorInt vec(number);
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorDouble VectorHelper::initVDouble(const double* values, int number)
{
  if (values == nullptr) return VectorDouble();
  VectorDouble vec(number);
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorVectorDouble VectorHelper::initVVDouble(const double* value, int n1, int n2)
{
  if (value == nullptr) return VectorVectorDouble();
  VectorVectorDouble vec;
  vec.resize(n1);
  for (int i1 = 0; i1 < n1; i1++) vec[i1].resize(n2);

  int lec = 0;
  for (int i1 = 0; i1 < n1; i1++)
    for (int i2 = 0; i2 < n2; i2++)
      vec[i1][i2] = value[lec++];
  return vec;
}

void VectorHelper::dump(const String &title, const VectorDouble& tab)
{
  std::stringstream sstr;
  if (!title.empty())
  {
    sstr << title.c_str() << std::endl;
  }
  sstr.precision(20);
  for (int i = 0, n = (int) tab.size(); i < n; i++)
    sstr << std::fixed << tab[i] << std::endl;
  messageFlush(sstr.str());
}

void VectorHelper::display(const String &title, const VectorDouble &vect, bool skipLine)
{
  if (!title.empty())
  {
    message("%s", title.c_str());
    if (skipLine) message("\n");
  }
  messageFlush(VH::toStringAsVD(vect));
}

void VectorHelper::display(const String &title, const VectorString &vect, bool skipLine)
{
  if (!title.empty())
  {
    message("%s", title.c_str());
    if (skipLine) message("\n");
  }
  messageFlush(VH::toStringAsVS(vect));
}

void VectorHelper::display(const String &title, const VectorVectorDouble &vect, bool skipLine)
{
  if (!title.empty())
  {
    message("%s", title.c_str());
    if (skipLine) message("\n");
  }
  messageFlush(VH::toStringAsVVD(vect));
}

void VectorHelper::display(const String &title, const VectorInt &vect, bool skipLine)
{
  if (!title.empty())
  {
    message("%s", title.c_str());
    if (skipLine) message("\n");
  }
  messageFlush(VH::toStringAsVI(vect));
}

String VectorHelper::toStringAsVD(const VectorDouble &vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toStringAsVVD(const VectorVectorDouble &vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toStringAsVS(const VectorString& vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toStringAsVI(const VectorInt &vec)
{
  return toVector(String(), vec);
}

void VectorHelper::displayStats(const String &title, const VectorDouble &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  double mean = 0.;
  double stdv = 0.;
  double mini = 1.e30;
  double maxi = -1.e30;

  for (int i = 0; i < ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    mean += value;
    stdv += value * value;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    mean /= (double) number;
    stdv = stdv / (double) number - mean * mean;
    stdv = (stdv > 0.) ? sqrt(stdv) : 0.;

    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %s\n", toDouble(mini).c_str());
    message("- Maximum  = %s\n", toDouble(maxi).c_str());
    message("- Mean     = %s\n", toDouble(mean).c_str());
    message("- St. Dev. = %s\n", toDouble(stdv).c_str());
  }
  else
  {
    message("No value defined\n");
  }
}

void VectorHelper::displayRange(const String &title, const VectorDouble &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  double mini = 1.e30;
  double maxi = -1.e30;

  for (int i = 0; i < ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %lf\n", mini);
    message("- Maximum  = %lf\n", maxi);
  }
  else
  {
    message("No value defined\n");
  }
}

void VectorHelper::displayRange(const String &title, const VectorInt &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  int mini =  100000000;
  int maxi = -100000000;

  for (int i = 0; i < ntotal; i++)
  {
    int value = vect[i];
    if (FFFF(value)) continue;
    number++;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %d\n", mini);
    message("- Maximum  = %d\n", maxi);
  }
  else
  {
    message("No value defined\n");
  }
}

void VectorHelper::displayNNZ(const String &title, const VectorDouble &vect, int nclass)
{
  int ntotal = (int) vect.size();
  VectorInt total(nclass);
  for (int ic = 0; ic < nclass; ic++) total[ic] = 0.;

  for (int i = 0; i < ntotal; i++)
  {
    double value = ABS(vect[i]);
    double tol = 1.;
    for (int ic = 0; ic < nclass; ic++)
    {
      tol /= 10.;
      if (value > tol) break;
      total[ic] += 1;
    }
  }

  if (!title.empty()) message("%s\n", title.c_str());
  for (int ic = 0; ic < nclass; ic++)
    message("Count below 10.e-%d = %d\n", ic+1, total[ic]);
}

bool VectorHelper::hasUndefined(const VectorDouble& vec)
{
  for (int i = 0, n = (int) vec.size(); i < n; i++)
    if (FFFF(vec[i])) return true;
  return false;
}

int VectorHelper::maximum(const VectorInt &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0;
  int max = -10000000;
  for (auto v : vec)
  {
    if (IFFFF(v)) continue;
    if (flagAbs) v = ABS(v);
    if (v > max) max = v;
  }
  return (max);
}

double VectorHelper::maximum(const VectorVectorDouble& vect, bool flagAbs)
{
  double val = VH::maximum(vect[0]);
  for (int i = 1, n = (int) vect.size(); i < n; i++)
    val = MAX(val,  VH::maximum(vect[i], flagAbs));
  return val;
}

int VectorHelper::minimum(const VectorInt &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0;
  int min = 10000000;
  for (auto v : vec)
  {
    if (IFFFF(v)) continue;
    if (flagAbs) v = ABS(v);
    if (v < min) min = v;
  }
  return (min);
}

/**
 * Calculate the maximum of 'vec' (conditionally to 'aux')
 * @param vec     Target vector
 * @param flagAbs When True, take the absolute value of 'vec' beforehand
 * @param aux     Conditional vector (see remarks)
 * @param mode    Comparison mode
 * @return Minimum value
 *
 * @remarks When 'aux' is defined, the maximum of 'vec' is conditional to 'aux'
 * @remarks - mode=0: statistics calculated only when 'vec' and 'aux' are both defined
 * @remarks - mode>0: statistics calculated only when 'vec' > 'aux'
 * @remarks - mode<0: statistics calculated only when 'vec' < 'aux'
 */
double VectorHelper::maximum(const VectorDouble &vec, bool flagAbs, const VectorDouble& aux, int mode)
{
  if (vec.empty()) return TEST;
  int size = (int) vec.size();
  bool flagAux = (! aux.empty() && (int) aux.size() == size);
  double max = -1.e30;

  if (! flagAux)
  {
    for (auto v : vec)
    {
      if (FFFF(v)) continue;
      if (flagAbs) v = ABS(v);
      if (v > max) max = v;
    }
  }
  else
  {
    const double *ptrv = &vec[0];
    double val_vec;
    const double *ptra = &aux[0];
    double val_aux;

    for (int i = 0; i < size; i++)
    {
      val_vec = (*ptrv);
      val_aux = (*ptra);
      if (! FFFF(val_vec) && ! FFFF(val_aux))
      {
        if (flagAbs) val_vec = ABS(val_vec);

        if (mode > 0)
        {
          if (val_aux > val_vec) continue;
        }
        if (mode < 0)
        {
          if (val_aux < val_vec) continue;
        }
        if (val_vec > max) max = val_vec;
      }
      ptrv++;
      ptra++;
    }
  }
  return (max);
}

/**
 * Calculate the minimum of 'vec' (conditionally to 'aux')
 * @param vec     Target vector
 * @param flagAbs When True, take the absolute value of 'vec' beforehand
 * @param aux     Conditional vector (see remarks)
 * @param mode    Comparison mode
 * @return Minimum value
 *
 * @remarks When 'aux' is defined, the minimum of 'vec' is conditional to 'aux'
 * @remarks - mode=0: statistics calculated only when 'vec' and 'aux' are both defined
 * @remarks - mode>0: statistics calculated only when 'vec' > 'aux'
 * @remarks - mode<0: statistics calculated only when 'vec' < 'aux'
 */
double VectorHelper::minimum(const VectorDouble &vec, bool flagAbs, const VectorDouble& aux, int mode)
{
  if (vec.empty()) return TEST;
  int size = (int) vec.size();
  bool flagAux = (! aux.empty() && (int) aux.size() == size);
  double min = 1.e30;

  if (! flagAux)
  {
    for (auto v : vec)
    {
      if (FFFF(v)) continue;
      if (flagAbs) v = ABS(v);
      if (v < min) min = v;
    }
  }
  else
  {
    const double *ptrv = &vec[0];
    double val_vec;
    const double *ptra = &aux[0];
    double val_aux;

    for (int i = 0; i < size; i++)
    {
      val_vec = (*ptrv);
      val_aux = (*ptra);
      if (! FFFF(val_vec) && ! FFFF(val_aux))
      {
        if (flagAbs) val_vec = ABS(val_vec);

        if (mode > 0)
        {
          if (val_aux > val_vec) continue;
        }
        if (mode < 0)
        {
          if (val_aux < val_vec) continue;
        }
        if (val_vec < min) min = val_vec;
      }
      ptrv++;
      ptra++;
    }
  }
  return (min);
}

double VectorHelper::minimum(const VectorVectorDouble& vect, bool flagAbs)
{
  double val = VH::minimum(vect[0]);
  for (int i = 1, n = (int) vect.size(); i < n; i++)
    val = MAX(val,  VH::minimum(vect[i], flagAbs));
  return val;
}

double VectorHelper::mean(const VectorDouble &vec)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  int number = 0;
  for (auto &v : vec)
  {
    if (FFFF(v)) continue;
    mean += v;
    number++;
  }
  if (number > 0)
    return (mean / (double) number);
  else
    return TEST;
}

int VectorHelper::cumul(const VectorInt& vec)
{
  int total = 0.;
  for (auto &v : vec)
  {
    total += v;
  }
  return total;
}

double VectorHelper::cumul(const VectorDouble& vec)
{
  double total = 0.;
  for (auto &v : vec)
  {
    if (FFFF(v)) continue;
    total += v;
  }
  return total;
}

double VectorHelper::variance(const VectorDouble &vec, bool scaleByN)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  double var = 0.;
  int number = 0;
  for (auto &v : vec)
  {
    if (FFFF(v)) continue;
    var += v * v;
    mean += v;
    number++;
  }
  if (number <= 1) return TEST;
  mean /= (double) number;
  if (scaleByN)
    var = var / (double) number - mean * mean;
  else
    var = (var - (double) number * mean * mean) / (double) (number-1);
  return var;
}

double VectorHelper::correlation(const VectorDouble &veca, const VectorDouble& vecb)
{
  if (veca.size() <= 0 || vecb.size() <= 0 || veca.size() != vecb.size()) return 0.;

  double m1  = 0.;
  double m2  = 0.;
  double v11 = 0.;
  double v22 = 0.;
  double v12 = 0.;
  int number = 0;
  for (int i = 0; i < (int) veca.size(); i++)
  {
    double z1 = veca[i];
    double z2 = vecb[i];
    if (FFFF(z1) || FFFF(z2)) continue;
    v11 += z1 * z1;
    v22 += z2 * z2;
    v12 += z1 * z2;
    m1 += z1;
    m2 += z2;
    number++;
  }
  if (number <= 0) return TEST;
  m1 /= (double) number;
  m2 /= (double) number;
  v11 = v11 / (double) number - m1 * m1;
  v22 = v22 / (double) number - m2 * m2;
  v12 = v12 / (double) number - m1 * m2;
  if (v11 <= 0.) return TEST;
  if (v22 <= 0.) return TEST;
  double corr = v12 / sqrt(v11 * v22);
  return corr;
}

/**
 * Calculate the quantiles
 * @param vec    Array of values
 * @param probas Array of probabilities (sorted by ascending order)
 * @return Vector of data values for the different probabilities
 */
VectorDouble VectorHelper::quantiles(const VectorDouble &vec,
                                     const VectorDouble &probas)
{
  int nproba = (int) probas.size();
  int nech   = (int) vec.size();
  if (nech <= 0 || nproba <= 0) return VectorDouble();

  VectorDouble retval(nproba, TEST);

  // Sort the data in ascending order
  VectorDouble sorted = VH::sort(vec, true);

  for (int ip = 0; ip < nproba; ip++)
  {
    double proba = probas[ip];
    int rank = (int) (proba * (double) nech);

    double value = TEST;
    if (rank < 0)
    {
      value = TEST;
    }
    else if (rank < nech - 1)
    {
      double v1 = sorted[rank];
      double v2 = sorted[rank + 1];
      double p1 = (double) rank / (double) nech;
      double p2 = (double) (1 + rank) / (double) nech;
      value = v1 + (proba - p1) * (v2 - v1) / (p2 - p1);
    }
    else
    {
      value = sorted[nech-1];
    }
    retval[ip] = value;
  }
  return retval;
}

double VectorHelper::stdv(const VectorDouble &vec, bool scaleByN)
{
  double var = variance(vec, scaleByN);
  if (!FFFF(var))
    return (sqrt(var));
  else
    return TEST;
}

double VectorHelper::norm(const VectorDouble &vec)
{
  double ip = innerProduct(vec, vec);
  return sqrt(ip);
}

double VectorHelper::norminf(const VectorDouble &vec)
{
  double norminf = 0.;
  for (int i = 0, nval = (int) vec.size(); i < nval; i++)
  {
    double value = ABS(vec[i]);
    if (value > norminf) norminf = value;
  }
  return norminf;
}

double VectorHelper::median(const VectorDouble &vec)
{
  VectorDouble med;
  for (int i = 0, n = (int) vec.size(); i < n; i++)
    if (! FFFF(vec[i])) med.push_back(vec[i]);

  // Sort the values
  med = sort(med);

  // Return the median value
  int number = (int) med.size();
  if (number <= 0) return TEST;
  if (isOdd(number))
    return med[number / 2];
  else
    return (med[number / 2] + med[number / 2 - 1]) / 2.;
}

double VectorHelper::normDistance(const VectorDouble &veca,
                                  const VectorDouble &vecb)
{
  double prod = 0.;
  double delta = 0.;
  const double *ptra = &veca[0];
  const double *ptrb = &vecb[0];
  for (int i = 0, n = (int) veca.size(); i < n; i++)
  {
    delta = (*ptra) - (*ptrb);
    prod += delta * delta;
    ptra++;
    ptrb++;
  }
  return sqrt(prod);
}

int VectorHelper::product(const VectorInt& vec)
{
  if (vec.empty()) return 0;
  int nprod = 1;
  const int* iptr = &vec[0];
  for (int i = 0, n = (int) vec.size(); i < n; i++)
  {
    nprod *= (*iptr);
    iptr++;
  }
  return nprod;
}

double VectorHelper::product(const VectorDouble& vec)
{
  if (vec.empty()) return 0;
  double nprod = 1.;
  const double* iptr = &vec[0];
  for (int i = 0, n = (int) vec.size(); i < n; i++)
  {
    nprod *= (*iptr);
    iptr++;
  }
  return nprod;
}

void VectorHelper::normalize(VectorDouble &vec)
{
  double ratio = VH::norm(vec);
  if (ratio <= 0.) return;
  for (auto &v : vec)
     v /= ratio;
}

void VectorHelper::normalize(double *tab, int ntab)
{
  int i;
  double norme;

  norme = 0.;
  for (i = 0; i < ntab; i++)
    norme += tab[i] * tab[i];
  norme = sqrt(norme);

  if (norme <= 0.) return;
  for (i = 0; i < ntab; i++)
    tab[i] /= norme;
  return;
}

void VectorHelper::normalizeFromGaussianDistribution(VectorDouble &vec,
                                                     double mini,
                                                     double maxi)
{
  double* iptr = &vec[0];
  for (int i = 0, n = (int) vec.size(); i < n; i++)
  {
    if (! FFFF(*iptr))
      (*iptr) = mini + (maxi - mini) * law_cdf_gaussian(*iptr);
    iptr++;
  }
}

VectorDouble VectorHelper::qnormVec(const VectorDouble& vec)
{
  int number = (int) vec.size();
  VectorDouble retvec(number, TEST);
  for (int i = 0; i < number; i++)
    retvec[i] = law_invcdf_gaussian(vec[i]);
  return retvec;
}

VectorDouble VectorHelper::pnormVec(const VectorDouble& vec)
{
  int number = (int) vec.size();
  VectorDouble retvec(number, TEST);
  for (int i = 0; i < number; i++)
    retvec[i] = law_cdf_gaussian(vec[i]);
  return retvec;
}

VectorDouble VectorHelper::normalScore(const VectorDouble &data,
                                       const VectorDouble &wt)
{
  int nech = (int) data.size();
  VectorDouble vec(nech, TEST);
  if (nech <= 0) return vec;

  // Check dimension of vector
  if (! wt.empty() && nech != (int) wt.size())
  {
    messerr("Arguments 'data' and 'wt' should have the same dimension");
    return VectorDouble();
  }

  // Check that weights of active samples are positive
  double wtotal = 0.;
  double nechtot = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    if (! FFFF(data[iech]))
    {
      double wtloc = 1.;
      if (!wt.empty()) wtloc = wt[iech];
      if (wtloc < 0.)
      {
        messerr("The weight of sample (%d) is negative (%lf)", iech + 1, wtloc);
        return VectorDouble();
      }
      wtotal += wtloc;
      nechtot += 1;
    }
  }
  if (wtotal <= 0.)
  {
    messerr("The sum of weights of active samples is not positive");
    return VectorDouble();
  }
  wtotal *= (1. + nechtot) / nechtot;

  // Get the list of indices sorted by increasing values of data
  VectorInt idx = VH::orderRanks(data);

  // Loop on the samples
  double wpartial = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    int jech = idx[iech];
    if (! FFFF(data[jech]))
    {
      double wtloc = (wt.empty()) ? 1. : wt[jech];
      wpartial += wtloc;
      vec[jech] = law_invcdf_gaussian(wpartial / wtotal);
    }
    else
    {
      vec[jech] = TEST;
    }
  }
  return vec;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (TEST if not defined)
 * @return
 */
bool VectorHelper::isConstant(const VectorDouble& vect, double refval)
{
  if (vect.empty()) return false;
  if (FFFF(refval)) refval = vect[0];
  const double* iptr = &vect[0];
  for (int i = 0, n = (int) vect.size(); i < n; i++)
  {
    if ((*iptr) != refval) return false;
    iptr++;
  }
  return true;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (ITEST if not defined)
 * @return
 */
bool VectorHelper::isConstant(const VectorInt& vect, int refval)
{
  if (vect.empty()) return false;
  if (IFFFF(refval)) refval = vect[0];
  const int* iptr = &vect[0];
  for (int i = 0, n = (int) vect.size(); i < n; i++)
  {
    if ((*iptr) != refval) return false;
    iptr++;
  }
  return true;
}

bool VectorHelper::isSame(const VectorDouble &v1, const VectorDouble &v2, double eps)
{
  if (v1.size() != v2.size()) return false;
  VectorDouble::const_iterator it1(v1.begin());
  VectorDouble::const_iterator it2(v2.begin());
  while (it1 < v1.end())
  {
    if (ABS(*it1 - *it2) > eps) return false;
    it1++;
    it2++;
  }
  return true;
}

bool VectorHelper::isSame(const VectorInt &v1, const VectorInt &v2)
{
  if (v1.size() != v2.size()) return false;
  VectorInt::const_iterator it1(v1.begin());
  VectorInt::const_iterator it2(v2.begin());
  while (it1 < v1.end())
  {
    if (*it1 != *it2) return false;
    it1++;
    it2++;
  }
  return true;
}

void VectorHelper::fill(VectorDouble &vec, double value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void VectorHelper::fill(VectorInt &vec, int value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void VectorHelper::fill(VectorVectorDouble &vec, double value)
{
  for (auto &e : vec)
  {
    VH::fill(e, value);
  }
}

void VectorHelper::fillUndef(VectorDouble& vec, double repl)
{
  VectorDouble::iterator it(vec.begin());
  while (it < vec.end())
  {
    if (FFFF(*it)) *it = repl;
    it++;
  }
}

/**
 * Create an output vector containing the 'number' consecutive numbers starting from 'ideb'
 *
 * @param number  Length of the output vector
 * @param ideb    Index of the first element of the output vector
 * @param step    Step between two consecutive values
 */
VectorInt VectorHelper::sequence(int number, int ideb, int step)
{
  VectorInt vec(number);

  int jdeb = ideb;
  for (int i = 0; i < number; i++)
  {
    vec[i] = jdeb;
    jdeb += step;
  }
  return vec;
}

/**
 * Create an output vector going from 'valFrom' to 'ValTo' by step of 'valStep'
 */

/**
 * Create a vector containing the a sequence of numbers
 * @param valFrom Starting value
 * @param valTo   Ending value
 * @param valStep Step
 * @param ratio   The whole sequence can be ultimately scaled by 'ratio'
 * @return
 */
VectorDouble VectorHelper::sequence(double valFrom,
                                    double valTo,
                                    double valStep,
                                    double ratio)
{
  VectorDouble vec;

  double value = valFrom;
  while (value <= valTo)
  {
    vec.push_back(value / ratio);
    value = value + valStep;
  }
  return vec;
}

VectorDouble VectorHelper::simulateUniform(int n, double mini, double maxi)
{
  VectorDouble vec(n);
  VectorDouble::iterator it(vec.begin());
  while (it < vec.end())
  {
    *it = law_uniform(mini, maxi);
    it++;
  }
  return vec;
}

VectorDouble VectorHelper::simulateBernoulli(int n, double proba, double vone, double velse)
{
  VectorDouble vec(n);
  VectorDouble::iterator it(vec.begin());
  while (it < vec.end())
  {
    double rand = law_uniform(0., 1.);
    if (rand < proba)
      *it= vone;
    else
      *it = velse;
    it++;
  }
  return vec;
}

VectorDouble VectorHelper::simulateGaussian(int n, double mean, double sigma)
{
  VectorDouble vec(n);
  simulateGaussianInPlace(vec,mean,sigma);
  return vec;
}

void VectorHelper::simulateGaussianInPlace(VectorDouble &vec,
                                           double mean,
                                           double sigma)
{
  VectorDouble::iterator it(vec.begin());
  while (it < vec.end())
  {
    *it= mean + sigma * law_gaussian();
    it++;
  }
}

VectorDouble VectorHelper::concatenate(const VectorDouble &veca,
                                       const VectorDouble &vecb)
{
  VectorDouble res = veca;
  for (auto &e: vecb)
    res.push_back(e);
  return res;
}

void VectorHelper::cumulate(VectorDouble &veca,
                            const VectorDouble &vecb,
                            double coeff,
                            double addval)
{
  if (veca.size() != vecb.size())
  {
    messerr("Arguments 'veca' and 'vecb' should have the same dimension. Nothing is done");
    return;
  }

  VectorDouble::iterator ita(veca.begin());
  VectorDouble::const_iterator itb(vecb.begin());
  while (ita < veca.end())
  {
    *ita += coeff * (*itb) + addval;
    ita++;
    itb++;
  }
}

/**
 * Display the first significant values of the input vector.
 * A "significant" value is a value larger than 'tol' in absolute value
 * Values are listed by decreasing importance.
 * @param vec  Input Vector
 * @param tol  Tolerance above which a value is significant (in absolute value)
 * @param nmax Limit on the number of values printed (-1: no limit)
 */
void VectorHelper::getMostSignificant(const VectorDouble& vec, double tol, int nmax)
{
  int nsize = (int) vec.size();
  VectorDouble absval(nsize, 0.);
  int ninvalid = 0;
  for (int i = 0; i < nsize; i++)
  {
    double value = vec[i];
    if (FFFF(value)) continue;
    value = ABS(value);
    if (value <= tol) continue;
    absval[i] = value;
    ninvalid++;
  }

  if (ninvalid <= 0) return;

  VectorInt ranks = orderRanks(absval, false);
  int nend = ninvalid;
  if (nmax > 0) nend = MIN(ninvalid, nmax);
  for (int i = 0; i < nend; i++)
  {
    int j = ranks[i];
    message("Sample %d - Value = %lf\n", j, vec[j]);
  }
  if (nmax > 0 && ninvalid > nmax)
    message("Found %d (out of %d) samples. Print limited to the %d most important ones.\n", ninvalid, nsize, nmax);
}

/**
 * Sample a set of 'ntotal' ranks (unique occurrence)
 * @param ntotal      Dimension to be sampled
 * @param proportion  Proportion of elected samples (in [0,1])
 * @param number      Number of elected samples
 * @param seed        Seed used for the random number generator
 * @param optSort     Sorting: 0 for None; 1 for Ascending; -1 for Descending
 * @return A vector of indices lying between 0 and ntotal-1. No duplicate.
 *
 * @remark If 'proportion' and 'number' are not specified,
 * @remark the output vector has dimension equal to 'ntotal'
 */
VectorInt VectorHelper::sampleRanks(int ntotal,
                                    double proportion,
                                    int number,
                                    int seed,
                                    int optSort)
{
  if (proportion <= 0. && number <= 0) return VectorInt();

  // Find the number of expected values
  int count;
  if (proportion <= 0. && number <= 0)
    count = ntotal;
  else if (proportion > 0.)
    count = (int) (ntotal * proportion);
  else
    count = number;
  count = MIN(ntotal, MAX(1, count));

  VectorInt ranks(ntotal);
  for (int i = 0; i < ntotal; i++) ranks[i] = i;

  shuffle (ranks.begin(), ranks.end(), std::default_random_engine(seed));

  ranks.resize(count);

  // Sort them out
  if (optSort > 0)
    ranks = sort(ranks, true);
  else if (optSort < 0)
    ranks = sort(ranks, false);

  std::vector<int>::iterator it;
  it = std::unique(ranks.begin(), ranks.end());
  ranks.resize(distance(ranks.begin(),it));

  return ranks;
}

VectorDouble VectorHelper::add(const VectorDouble &veca, const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
  {
    messerr("Arguments 'veca' and 'vecb' should have the same dimension. Nothing is done");
    return veca;
  }

  VectorDouble res(veca.size());
  VectorDouble::iterator it(res.begin());
  VectorDouble::const_iterator ita(veca.begin());
  VectorDouble::const_iterator itb(vecb.begin());

  while (it < res.end())
  {
    *it = *ita + *itb;
    it++;
    ita++;
    itb++;
  }
  return res;
}

/**
 * Performs: veca += vecb
 * @param dest Input/Output vector
 * @param src Auxiliary vector
 */
void VectorHelper::addInPlace(VectorDouble &dest, const VectorDouble &src)
{
  if (dest.size() != src.size())
  {
    messerr("Arguments 'dest' and 'src' should have the same dimension. Nothing is done");
    return;
  }

  VectorDouble::iterator itd(dest.begin());
  VectorDouble::const_iterator its(src.begin());
  while (itd < dest.end())
  {
    *itd += *its;
    itd++;
    its++;
  }
}

void VectorHelper::addInPlace(const VectorDouble &veca,
                              const VectorDouble &vecb,
                              VectorDouble &res,
                              int size)
{
  if (size <= 0) size = (int) veca.size();
  if (size != (int) vecb.size())
    my_throw("Wrong size");
  if ((int) res.size() != size) res.resize(size);

  const double* iptra = &veca[0];
  const double* iptrb = &vecb[0];
  double* iptrv = &res[0];
  for (int i = 0; i < size; i++)
  {
    (*iptrv) = (*iptra) + (*iptrb);
    iptrv++;
    iptra++;
    iptrb++;
  }
}

void VectorHelper::addInPlace(const double* veca,
                              const double* vecb,
                              double* res,
                              int size)
{
  for (int i = 0; i < size; i++)
    res[i] = veca[i] + vecb[i];
}

void VectorHelper::addInPlace(const VectorVectorDouble &in1,
                              const VectorVectorDouble &in2,
                              VectorVectorDouble &outv)
{
  for (int is = 0, ns = (int) in1.size(); is < ns; is++)
  {
    for (int i = 0, n = (int) in1[is].size(); i < n; i++)
    {
      outv[is][i] = in2[is][i] + in1[is][i];
    }
  }
}

/**
 * Return a vector containing vecb - veca
 * @param veca Input Vector
 * @param vecb Input Vector
 * @return
 */
VectorDouble VectorHelper::subtract(const VectorDouble &veca,
                                    const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
    my_throw("Wrong size");

  VectorDouble res(veca.size());
  VectorDouble::iterator it(res.begin());
  VectorDouble::const_iterator ita(veca.begin());
  VectorDouble::const_iterator itb(vecb.begin());
  while (it < res.end())
  {
    *it = *itb - *ita;
    it++;
    ita++;
    itb++;
  }
  return res;
}

/**
 * Return a vector containing vecb - veca (integer version)
 * @param veca Input Vector
 * @param vecb Input Vector
 * @return
 */
VectorInt VectorHelper::subtract(const VectorInt &veca,
                                 const VectorInt &vecb)
{
  if (veca.size() != vecb.size())
    my_throw("Wrong size");

  VectorInt res(veca.size());
  VectorInt::iterator it(res.begin());
  VectorInt::const_iterator ita(veca.begin());
  VectorInt::const_iterator itb(vecb.begin());
  while (it < res.end())
  {
    *it = *itb - *ita;
    it++;
    ita++;
    itb++;
  }
  return res;
}

/**
 * Performs: veca -= vecb
 * @param dest Input/Output vector
 * @param src Auxiliary vector
 */
void VectorHelper::subtractInPlace(VectorDouble &dest, const VectorDouble &src)
{
  VectorDouble res;
  if (dest.size() != src.size())
    my_throw("Wrong size");

  VectorDouble::iterator itd(dest.begin());
  VectorDouble::const_iterator its(src.begin());
  while (itd < dest.end())
  {
    *itd -= *its;
    itd++;
    its++;
  }
}

void VectorHelper::subtractInPlace(VectorInt &dest, const VectorInt &src)
{
  VectorInt res;
  if (dest.size() != src.size())
    my_throw("Wrong size");

  VectorInt::iterator itd(dest.begin());
  VectorInt::const_iterator its(src.begin());
  while (itd < dest.end())
  {
    *itd -= *its;
    itd++;
    its++;
  }
}

void VectorHelper::subtractInPlace(const VectorVectorDouble &in1,
                                   const VectorVectorDouble &in2,
                                   VectorVectorDouble &outv)
{
  for (int is = 0, ns = (int) in1.size(); is < ns; is++)
  {
    for (int i = 0, n = (int) in1[is].size(); i < n; i++)
    {
      outv[is][i] = in2[is][i] - in1[is][i];
    }
  }
}


void VectorHelper::multiplyInPlace(VectorDouble &vec, const VectorDouble &v)
{
  if (vec.size() != v.size())
    my_throw("Arguments 'vec' and 'v' should have same dimension");

  VectorDouble::iterator it(vec.begin());
  VectorDouble::const_iterator itv(v.begin());
  while (it < vec.end())
  {
    *it *= (*itv);
    it++;
    itv++;
  }
}

void VectorHelper::divideInPlace(VectorDouble &vec, const VectorDouble &v)
{
  if (vec.size() != v.size())
    my_throw("Arguments 'vec' and 'v' should have same dimension");

  VectorDouble::iterator it(vec.begin());
  VectorDouble::const_iterator itv(v.begin());
  while (it < vec.end())
  {
    if (ABS(*itv) >= EPSILON20)
      *it /= (*itv);
    it++;
    itv++;
  }
}

void VectorHelper::multiplyConstant(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d *= v;});
}

void VectorHelper::multiplyConstantInPlace(const VectorDouble &vecin, double v, VectorDouble& vecout)
{
  VectorDouble::iterator itout(vecout.begin());
  VectorDouble::const_iterator itin(vecin.begin());
  while (itin < vecin.end())
  {
    *itout = (*itin) * v;
    itin++;
    itout++;
  }
}

void VectorHelper::addMultiplyConstantInPlace(double val1,
                                              const VectorVectorDouble &in1,
                                              VectorVectorDouble &outv)
{
  for (int is = 0, ns = (int) in1.size(); is < ns; is++)
  {
    for (int i = 0, n = (int) in1[is].size(); i < n; i++)
    {
      outv[is][i] += val1 * in1[is][i];
    }
  }
}

void VectorHelper::divideConstant(VectorDouble &vec, double v)
{
  if (isZero(v))
  my_throw("division by 0");
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d /= v; });
}

void VectorHelper::copy(const VectorDouble &vecin, VectorDouble &vecout, int size)
{
  if (size < 0) size = (int) vecin.size();
  if (size > (int) vecout.size())
    my_throw("Wrong size");

  VectorDouble::iterator itout(vecout.begin());
  VectorDouble::const_iterator itin(vecin.begin());
  for (int i = 0; i < size; i++)
  {
    (*itout) = (*itin);
    itin++;
    itout++;
  }
}

void VectorHelper::copy(const VectorInt &vecin, VectorInt &vecout, int size)
{
  if (size < 0) size = (int) vecin.size();
  if (size > (int) vecout.size())
    my_throw("Wrong size");

  VectorInt::iterator itout(vecout.begin());
  VectorInt::const_iterator itin(vecin.begin());
  for (int i = 0; i < size; i++)
  {
    (*itout) = (*itin);
    itin++;
    itout++;
  }
}

void VectorHelper::copy(const VectorVectorDouble &inv, VectorVectorDouble &outv)
{
  for (int is = 0, ns = (int) inv.size(); is < ns; is++)
  {
    for (int i = 0, n = (int) inv[is].size(); i < n; i++)
    {
      outv[is][i] = inv[is][i];
    }
  }
}

void VectorHelper::addConstant(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d += v;});
}

void VectorHelper::addConstant(VectorInt &vec, int v)
{
  std::for_each(vec.begin(), vec.end(), [v](int &d)
  { d += v;});
}

VectorDouble VectorHelper::power(const VectorDouble &vec, double power)
{
  VectorDouble res(vec.size());
  VectorDouble::iterator it(res.begin());
  VectorDouble::const_iterator itv(vec.begin());
  while (it < res.end())
  {
    *it = pow(*itv, power);
    it++;
    itv++;
  }
  return res;
}

VectorDouble VectorHelper::inverse(const VectorDouble& vec)
{
  VectorDouble inv(vec.size());
  VectorDouble::iterator it(inv.begin());
  VectorDouble::const_iterator itv(vec.begin());
  while (it < inv.end())
  {
    *it = 1. / *itv;
    it++;
    itv++;
  }
  return inv;
}

int VectorHelper::countUndefined(const VectorDouble &vec)
{
  int count = 0;
  VectorDouble::const_iterator it(vec.begin());
  while (it < vec.end())
  {
    if (FFFF(*it)) count++;
    it++;
  }
  return count;
}

int VectorHelper::countDefined(const VectorDouble &vec)
{
  int count = 0;
  VectorDouble::const_iterator it(vec.begin());
  while (it < vec.end())
  {
    if (! FFFF(*it)) count++;
    it++;
  }
  return count;
}

/**
 * Calculate the diagonal of the box extension
 * @param mini Array of lower coordinates of the box
 * @param maxi Array of upper coordinates of the box
 * @return
 * @remark If one coordinate is undefined, TEST is returned.
 */
double VectorHelper::extensionDiagonal(const VectorDouble &mini,
                                       const VectorDouble &maxi)
{
  double diag = 0.;
  VectorDouble delta = VH::subtract(mini, maxi);
  int ndim = (int) delta.size();
  for (int idim = 0; idim < ndim; idim++)
  {
    double dval = delta[idim];
    if (FFFF(dval)) return TEST;
    diag += dval * dval;
  }
  diag = sqrt(diag);
  return diag;
}

VectorInt VectorHelper::unique(const VectorInt& vecin, int size)
{
  if (size < 0) size = (int) vecin.size();

  VectorInt vecout = vecin;
  vecout.resize(size);
  std::sort(vecout.begin(), vecout.end());
  auto last = std::unique(vecout.begin(), vecout.end());
  vecout.erase(last, vecout.end());
  return vecout;
}

VectorDouble VectorHelper::unique(const VectorDouble& vecin, int size)
{
  if (size < 0) size = (int) vecin.size();

  VectorDouble vecout = vecin;
  vecout.resize(size);
  std::sort(vecout.begin(), vecout.end());
  auto last = std::unique(vecout.begin(), vecout.end());
  vecout.erase(last, vecout.end());
  return vecout;
}

bool VectorHelper::isInList(const VectorInt& vec, int item)
{
  return std::count(vec.begin(), vec.end(), item);
}

VectorInt VectorHelper::sort(const VectorInt& vecin, bool ascending, int size)
{
  if (vecin.empty()) return VectorInt();
  if (size < 0) size = (int) vecin.size();
  VectorInt vecout = vecin;
  vecout.resize(size);
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

VectorDouble VectorHelper::sort(const VectorDouble& vecin, bool ascending, int size)
{
  if (vecin.empty()) return VectorDouble();
  if (size < 0) size = (int) vecin.size();
  VectorDouble vecout = vecin;
  vecout.resize(size);
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

void VectorHelper::sortInPlace(VectorInt& vecin, bool ascending, int size)
{
  if (vecin.empty()) return;
  VectorInt vecout = sort(vecin, ascending, size);
  copy(vecout, vecin, size);
}

void VectorHelper::sortInPlace(VectorDouble& vecin, bool ascending, int size)
{
  if (vecin.empty()) return;
  VectorDouble vecout = sort(vecin, ascending, size);
  copy(vecout, vecin, size);
}

bool VectorHelper::isSorted(const VectorDouble& vec, bool ascending)
{
  int nval = (int) vec.size();

  if (ascending)
  {
    // Ascending order
    for (int i = 1; i < nval; i++)
    {
      if (vec[i] > vec[i - 1]) continue;
      return false;
    }
  }
  else
  {
    // Descending order
    for (int i = 1; i < nval; i++)
    {
      if (vec[i] < vec[i - 1]) continue;
      return false;
    }
  }
  return true;
}

/**
 * From an input list, filter out all the elements which do no lie within [vmin, vmax[,
 * suppress double occurrences and sort them out (ascending or descending)
 * @param vecin Input array (integer)
 * @param vmin  lower bound included (or ITEST)
 * @param vmax  upper bound excluded (or ITEST)
 * @param ascending True for ascending order; False for descending order
 * @return Output array (integers)
 */
VectorInt VectorHelper::filter(const VectorInt &vecin,
                               int vmin,
                               int vmax,
                               bool ascending)
{
  VectorInt vecout = vecin;

  // Sort the vector
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());

  // Unique occurrence
  std::vector<int>::iterator it;
  it = std::unique(vecout.begin(), vecout.end());
  vecout.resize(distance(vecout.begin(),it));

  // Filter out the irrelevant values
  int nech = (int) vecout.size();
  for (int j = 0; j < nech; j++)
  {
    int i = nech - j - 1;
    if (!IFFFF(vmin))
    {
      if (vecout[i] < vmin)
      {
        vecout.erase(vecout.begin()+i);
        continue;
      }
    }
    if (!IFFFF(vmax))
    {
      if (vecout[i] >= vmax)
      {
        vecout.erase(vecout.begin()+i);
        continue;
      }
    }
  }
  return vecout;
}

/**
 * Returns the list complementary to 'sel' within 'vecin'
 * @param vec Initial list
 * @param sel Vector of forbidden elements
 * @return Complementary list
 */
VectorInt VectorHelper::complement(const VectorInt& vec, const VectorInt& sel)
{
  VectorInt rest;
  if (vec.empty()) return rest;
  if (sel.empty()) return vec;

  // Sort

  VectorInt allVec = vec;
  std::sort(allVec.begin(), allVec.end());

  VectorInt offVec = sel;
  std::sort(offVec.begin(), offVec.end());

  int j, k, idx;
  int nvec = (int) allVec.size();
  int noff = (int) offVec.size();
  for (int i = 0; i < nvec; i++)
  {
    j = allVec.at(i);

    // I go through offVec as long as element is strictly less than j
    k = 0;
    idx = offVec.at(k);
    while (idx < j && k < noff)
    {
        idx = offVec.at(k++);
    }

    if (idx != j) // idx not in offElemsVec
    {
      rest.push_back(j);
    }
  }
  return rest;
}

/**
 * Returns the permutation which rearranges the input vector into any order
 * @param vecin Input vector
 * @param ascending True for ascending order; False otherwise
 * @param size Optional dimension of the input vector
 * @return Vector of orders
 */
VectorInt VectorHelper::orderRanks(const VectorInt& vecin, bool ascending, int size)
{
  if (vecin.empty()) return VectorInt();

  if (size < 0) size = (int) vecin.size();
  VectorInt idx(size);
  for (int i = 0; i < size; i++) idx[i] = i;

  // sort indexes based on comparing values in v using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings when v contains elements of equal values
  VectorT<int>::iterator last = idx.begin() + size;
  if (ascending)
  {
    stable_sort(idx.begin(), last,
                [&vecin](size_t i1, size_t i2) {return vecin[i1] < vecin[i2];});
  }
  else
  {
    stable_sort(idx.begin(), last,
                [&vecin](size_t i1, size_t i2) {return vecin[i1] > vecin[i2];});
  }

  return idx;
}

/**
 * Returns the permutation which rearranges the input vector into any order
 * @param vecin Input vector
 * @param ascending True for ascending order; False otherwise
 * @param size Optional vector dimension
 * @return Vector of orders
 */
VectorInt VectorHelper::orderRanks(const VectorDouble& vecin, bool ascending, int size)
{
  if (vecin.empty()) return VectorInt();
  if (size < 0) size = (int) vecin.size();
  VectorInt idx(size);
  for (int i = 0; i < size; i++) idx[i] = i;

  // sort indexes based on comparing values in v using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings when v contains elements of equal values
  VectorT<int>::iterator last = idx.begin() + size;
  if (ascending)
  {
    stable_sort(idx.begin(), last,
                [&vecin](size_t i1, size_t i2) {return vecin[i1] < vecin[i2];});
  }
  else
  {
    stable_sort(idx.begin(), last,
                [&vecin](size_t i1, size_t i2) {return vecin[i1] > vecin[i2];});
  }

  return idx;
}

VectorInt VectorHelper::sortRanks(const VectorDouble& vecin, bool ascending, int size)
{
  if (vecin.empty()) return VectorInt();
  if (size < 0) size = (int) vecin.size();
  VectorInt order = orderRanks(vecin, ascending, size);
  VectorInt idx(size);
  for (int i = 0; i < size; i++) idx[order[i]] = i;

  return idx;
}

VectorInt VectorHelper::reorder(const VectorInt& vecin, const VectorInt& order, int size)
{
  if (size < 0) size = (int) vecin.size();
  VectorInt vecout(size);
  for (int i = 0; i< size; i++)
    vecout[i] = vecin[order[i]];
  return vecout;
}

VectorDouble VectorHelper::reorder(const VectorDouble& vecin, const VectorInt& order, int size)
{
  if (size < 0) size = (int) vecin.size();
  VectorDouble vecout(size);
  for (int i = 0; i< size; i++)
    vecout[i] = vecin[order[i]];
  return vecout;
}

/*****************************************************************************/
/*!
 **  Sorts the (double) array value() and the array ranks() if provided
 **
 ** \param[in]  safe   1 if the value array if preserved
 **                    0 if the value array is also sorted
 ** \param[in]  ascending Sorting order
 ** \param[in]  size   Optional vector dimension
 **
 ** \param[out] ranks  input and output int array
 ** \param[out] values input and output double array
 **
 ** \remark  If ranks = NULL, ranks is ignored
 ** \remark  When using 'size', the remaining part of arrays is unchanged
 **
 *****************************************************************************/
void VectorHelper::arrangeInPlace(int safe,
                                  VectorInt &ranks,
                                  VectorDouble &values,
                                  bool ascending,
                                  int size)
{
  VectorInt order = orderRanks(values, ascending, size);

  if (! ranks.empty())
  {
    VectorInt newranks = reorder(ranks, order, size);
    copy(newranks, ranks, size);
  }

  if (safe == 0)
  {
    VectorDouble newvals = reorder(values, order, size);
    copy(newvals, values, size);
  }
}

/*****************************************************************************/
/*!
 **  Sorts the (int) array value() and the array ranks() if provided
 **
 ** \param[in]  safe   1 if the value array if preserved
 **                    0 if the value array is also sorted
 ** \param[in] ascending True for ascending order; False for descending order
 ** \param[in] size    Optional size
 **
 ** \param[out] ranks  intput and output int array
 ** \param[out] values input and output int array
 **
 ** \remark  If ranks = NULL, ranks is ignored
 ** \remark  When using 'size', the remaining part of arrays is unchanged
 **
 *****************************************************************************/
void VectorHelper::arrangeInPlace(int safe,
                                  VectorInt &ranks,
                                  VectorInt &values,
                                  bool ascending,
                                  int size)
{
  VectorInt order = orderRanks(values, ascending, size);

  if (! ranks.empty())
  {
    VectorInt newranks = reorder(ranks, order, size);
    ranks = newranks;
  }

  if (safe == 0)
  {
    VectorInt newval = reorder(values, order, size);
    values = newval;
  }
}

std::pair<double,double> VectorHelper::rangeVals(const VectorDouble& vec)
{
  std::pair<double,double> res(vec[0],vec[0]);
  for (int i = 1; i < (int)vec.size(); i++)
  {
    res.first  = MIN(res.first ,vec[i]);
    res.second = MAX(res.second,vec[i]);
  }
  return res;
}

double VectorHelper::innerProduct(const VectorDouble &veca,
                                  const VectorDouble &vecb,
                                  int size)
{
  if (size < 0) size = (int) veca.size();
  if (size > (int) veca.size() || size > (int) vecb.size())
    my_throw("Incompatible sizes");

  return innerProduct(veca.data(), vecb.data(), size);
}

double VectorHelper::innerProduct(const double* veca,
                                  const double* vecb,
                                  int size)
{
  double prod = 0.;
  const double *ptra = &veca[0];
  const double *ptrb = &vecb[0];
  for (int i = 0; i < size; i++)
  {
    prod += (*ptra) * (*ptrb);
    ptra++;
    ptrb++;
  }
  return prod;
}

double VectorHelper::innerProduct(const VectorVectorDouble &x,
                                  const VectorVectorDouble &y)
{
  double s = 0.;
  for (int i = 0, n = (int) x.size(); i < n; i++)
    s += VH::innerProduct(x[i], y[i]);
  return s;
}

/**
 * Cross product (limited to 3D)
 * @param veca First vector
 * @param vecb Second Vector
 * @return
 */
VectorDouble VectorHelper::crossProduct3D(const VectorDouble &veca,
                                          const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
    my_throw("Wrong size");
  VectorDouble res(3);
  crossProduct3DInPlace(veca.data(), vecb.data(), res.data());
  return res;
}

void VectorHelper::crossProduct3DInPlace(const double *a, const double *b, double *v)
{
  v[0] = a[1] * b[2] - a[2] * b[1];
  v[1] = a[2] * b[0] - a[0] * b[2];
  v[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * Method which flattens a VectorVectorDouble into a VectorDouble
 * @param vvd Input VectorVectorDouble
 * @return Returned VectorDouble
 */
VectorDouble VectorHelper::flatten(const VectorVectorDouble& vvd)
{
  VectorDouble vd;

  for (int i = 0; i < (int) vvd.size(); i++)
    for (int j = 0; j < (int) vvd[i].size(); j++)
      vd.push_back(vvd[i][j]);

  return vd;
}

void VectorHelper::flattenInPlace(const VectorVectorDouble& vvd, VectorDouble& vd)
{
  int ecr = 0;
  for (int i = 0; i < (int) vvd.size(); i++)
    for (int j = 0; j < (int) vvd[i].size(); j++)
      vd[ecr++] = (vvd[i][j]);
}

VectorVectorDouble VectorHelper::unflatten(const VectorDouble& vd, const VectorInt& sizes)
{
  VectorVectorDouble vvd;

  int lec = 0;
  for (int i = 0, n = (int) sizes.size(); i < n; i++)
  {
    int lng = sizes[i];
    VectorDouble local(lng);
    for (int j = 0; j < lng; j++)
      local[j] = vd[lec++];
    vvd.push_back(local);
  }
  return vvd;
}

void VectorHelper::unflattenInPlace(const VectorDouble& vd, VectorVectorDouble& vvd)
{
  int lec = 0;
  for (int i = 0, n = (int) vvd.size(); i < n; i++)
    for (int j = 0; j < (int) vvd[i].size(); j++)
      vvd[i][j] = vd[lec++];
}

VectorDouble VectorHelper::suppressTest(const VectorDouble& vecin)
{
  VectorDouble vecout;
  for (int i = 0, n = (int) vecin.size(); i < n; i++)
  {
    if (! FFFF(vecin[i])) vecout.push_back(vecin[i]);
  }
  return vecout;
}

void VectorHelper::linearCombinationInPlace(double val1,
                                            const VectorDouble &vd1,
                                            double val2,
                                            const VectorDouble &vd2,
                                            VectorDouble &outv)
{
  if (vd1.empty() || vd2.empty()) return;
  for (int i = 0, n = (int) vd1.size(); i < n; i++)
  {
    double value = 0.;
    if (val1 != 0. && !vd1.empty()) value += val1 * vd1[i];
    if (val2 != 0. && !vd2.empty()) value += val2 * vd2[i];
    outv[i] = value;
  }
}

void VectorHelper::linearCombinationVVDInPlace(double val1,
                                               const VectorVectorDouble &vvd1,
                                               double val2,
                                               const VectorVectorDouble &vvd2,
                                               VectorVectorDouble &outv)
{
  if (vvd1.empty() || vvd2.empty()) return;

  for (int is = 0, ns = (int) vvd1.size(); is < ns; is++)
  {
    for (int i = 0, n = (int) vvd1[is].size(); i < n; i++)
    {
      double value = 0.;
      if (val1 != 0. && ! vvd1.empty()) value += val1 * vvd1[is][i];
      if (val2 != 0. && ! vvd2.empty()) value += val2 * vvd2[is][i];
      outv[is][i] = value;
    }
  }
}

/**
 * Extract the part of a vector 'vecin' (dimensioned to dimension of 'vecout')
 * starting at address 'istart' and copy it into 'vecout'
 * @param vecin  Initial vector
 * @param vecout Resulting vector (already allocated)
 * @param start  Starting address (within 'vecin')
 */
void VectorHelper::extractInPlace(const VectorDouble& vecin, VectorDouble& vecout, int start)
{
  std::copy(vecin.begin() + start, vecin.begin() + start + vecout.size(), vecout.begin());
}

/**
 * Merge 'vecin' into 'vecout' starting at address 'istart'
 * @param vecin  Initial vector
 * @param vecout Vector where 'vecin' should be copied
 * @param start  Starting address (in 'vecout')
 */
void VectorHelper::mergeInPlace(const VectorDouble& vecin, VectorDouble& vecout, int start)
{
  std::copy(vecin.begin(), vecin.end(), vecout.begin() + start);
}

/**
 * Transform a vector of double values as follows
 * @param tab  Vector of double values
 * @param oper_choice Operation on the diagonal term (see Utilities::operate_XXX)
 * @return
 */
void VectorHelper::transformVD(VectorDouble& tab, int oper_choice)
{
  operate_function oper_func = operate_Identify(oper_choice);
  int number = (int) tab.size();
  for (int i = 0; i < number; i++)
    tab[i] = oper_func(tab[i]);
}

/****************************************************************************/
/*!
 **  Fix plausible values for the Direction coefficients.
 **  They must be defined and with norm equal to 1
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in,out]  codir Input/Output Direction coefficients
 **
 *****************************************************************************/
void VectorHelper::normalizeCodir(int ndim, VectorDouble &codir)
{
  double norme;

  if (codir.empty()) return;
  norme = VH::innerProduct(codir, codir, ndim);
  if (norme <= 0.)
  {
    for (int idim = 0; idim < ndim; idim++)
      codir[idim] = 0.;
    codir[0] = 1.;
  }
  else
  {
    norme = sqrt(norme);
    for (int idim = 0; idim < ndim; idim++)
      codir[idim] /= norme;
  }
}

/**
 * Operate the squeeze-and-stretch algorithm forward (see remarks)
 * @param vecin  Input vector (in structural system)
 * @param vecout Output vector (in sugar box system)
 * @param origin Origin of the vertical axis (structural system)
 * @param mesh   Mesh of the vertical axis (structural system)
 * @param top    Elevation of the Top surface
 * @param bot    Elevation of the Bottom surface
 *
 * @remarks The information is contained in 'vecin' which is defined on a regular 1D grid
 * @remarks in the structural system. The purpose is to sample the relevant sub-information
 * @remarks (between 'top' and 'bot') densely in 'vecout'
 */
void VectorHelper::squeezeAndStretchInPlaceForward(const VectorDouble &vecin,
                                                   VectorDouble &vecout,
                                                   double origin,
                                                   double mesh,
                                                   double top,
                                                   double bot)
{
  int nzin  = (int) vecin.size();
  int nzout = (int) vecout.size();
  double thick = top - bot;
  double ratio = thick / nzout;

  // Loop on the positions of the pile in the sugar box system
  for (int iz = 0; iz < nzout; iz++)
  {
    // Corresponding coordinate of the sample in the structural system
    double zzin = bot + (double) iz * ratio;

    // Find the index in the input vector
    int izin = (int) ((zzin - origin) / mesh);
    if (izin < 0 || izin >= nzin) continue;

    // Assign the value
    vecout[iz] = vecin[izin];
  }
}

/**
 * Operate the squeeze-and-stretch algorithm backward (see remarks)
 * @param vecin  Input vector (in sugar box system)
 * @param vecout Output vector (in structural system)
 * @param origin Origin of the vertical axis (structural system)
 * @param mesh   Mesh of the vertical axis (structural system)
 * @param top    Elevation of the Top surface
 * @param bot    Elevation of the Bottom surface
 *
 * @remarks The information is contained in 'vecin' which is defined on a regular 1D grid
 * @remarks (characterized by 'base' and 'mesh')
 * @remarks Extend the relevant information, lying between 'bot' and 'top' in order to fill
 * @remarks the whole vector 'vecout'
 */
void VectorHelper::squeezeAndStretchInPlaceBackward(const VectorDouble &vecin,
                                                    VectorDouble &vecout,
                                                    double origin,
                                                    double mesh,
                                                    double top,
                                                    double bot)
{
  int nzin  = (int) vecin.size();
  int nzout = (int) vecout.size();

  // Blank out the output vector
  vecout.fill(TEST);
  double thick = top - bot;
  if (thick <= 0) return;

  // Get the top and bottom indices in the output vector
  int indbot = floor((bot - origin) / mesh);
  if (indbot < 0) indbot = 0;
  int indtop = ceil((top - origin)  / mesh);
  if (indtop >= nzout) indtop = nzout - 1;

  double ratio = (double) nzin / thick;

  // Loop on the positions of the pile in the structural system
  for (int izout = indbot; izout <= indtop; izout++)
  {
    // Get the location of the sample in the structural system
    double zzout = origin + izout * mesh;

    // Find the index in the input vector (sugar box)
    int izin = ratio * (zzout - bot);
    if (izin < 0 || izin >= nzin) continue;

    // Assign the value
    vecout[izout] = vecin[izin];
  }
}

/*****************************************************************************/
/*!
 **  Find the location of the minimum value within a vector
 **
 ** \return Rank of the minimum value
 **
 ** \param[in]  tab  Vector of values
 **
 *****************************************************************************/
int VectorHelper::whereMinimum(const VectorDouble& tab)
{
  int ibest = -1;
  double vbest = 1.e30;
  for (int i = 0, ntab = (int) tab.size(); i < ntab; i++)
  {
    if (FFFF(tab[i])) continue;
    if (tab[i] > vbest) continue;
    vbest = tab[i];
    ibest = i;
  }
  return ibest;
}

/*****************************************************************************/
/*!
 **  Find the location of the maximum value within a vector
 **
 ** \return Rank of the maximum value
 **
 ** \param[in]  tab  Vector of values
 **
 *****************************************************************************/
int VectorHelper::whereMaximum(const VectorDouble& tab)
{
  int ibest = -1;
  double vbest = -1.e30;
  for (int i = 0, ntab = (int) tab.size(); i < ntab; i++)
  {
    if (FFFF(tab[i])) continue;
    if (tab[i] < vbest) continue;
    vbest = tab[i];
    ibest = i;
  }
  return ibest;
}

VectorDouble VectorHelper::reduceOne(const VectorDouble &vecin, int index)
{
  VectorInt vindex(1);
  vindex[0] = index;
  return reduce(vecin, vindex);
}

VectorDouble VectorHelper::reduce(const VectorDouble &vecin, const VectorInt& vindex)
{
  VectorDouble vecout = vecin;

  // Sort the indices to be removed in ascending order
  VectorInt indexLocal = vindex;
  std::sort(indexLocal.begin(), indexLocal.end());

  int nsel = (int) indexLocal.size();
  for (int j = 0; j < nsel; j++)
  {
    int i = indexLocal[nsel - j - 1];
    vecout.erase(vecout.begin()+i);
  }
  return vecout;
}

void VectorHelper::roundDecimalsInPlace(VectorDouble& vec, int ndec)
{
  for (int i = 0, n = (int) vec.size(); i < n; i++)
  {
    if (FFFF(vec[i])) continue;
    vec[i] = roundDecimals(vec[i], ndec);
  }
}

void VectorHelper::roundDigitsInPlace(VectorDouble& vec, int ndec)
{
  for (int i = 0, n = (int) vec.size(); i < n; i++)
  {
    if (FFFF(vec[i])) continue;
    vec[i] = roundDigits(vec[i], ndec);
  }
}
