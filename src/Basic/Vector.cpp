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

#include "Basic/Vector.hpp"
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

VectorInt ut_vector_int(int nval, int value)
{
  VectorInt tab(nval, value);
  return tab;
}

VectorDouble ut_vector_double(int nval, double value)
{
  VectorDouble tab(nval, value);
  return tab;
}

VectorVectorDouble ut_vector_vector_double(int nval1, int nval2, double value)
{
  VectorVectorDouble tab(nval1, VectorDouble(nval2, value));
  return tab;
}

VectorVectorInt ut_vector_vector_int(int nval1, int nval2, int value)
{
  VectorVectorInt tab(nval1, VectorInt(nval2, value));
  return tab;
}

String ut_vector_string(const VectorDouble &vec)
{
  return toVector(String(), vec);
}
String ut_vector_string(const VectorVectorDouble &vec)
{
  return toVector(String(), vec);
}
String ut_vector_string(const VectorString& vec)
{
  return toVector(String(), vec);
}

String ut_ivector_string(const VectorInt &vec)
{
  return toVector(String(), vec);
}

void ut_vector_display(const String &title, const VectorDouble &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_vector_string(vect));
}

void ut_vector_display(const String &title, const VectorString &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_vector_string(vect));
}

void ut_vector_display(const String &title, const VectorVectorDouble &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_vector_string(vect));
}

void ut_vector_display_stats(const String &title, const VectorDouble &vect)
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

void ut_vector_display_range(const String &title, const VectorDouble &vect)
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

void ut_ivector_display(const String &title, const VectorInt &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_ivector_string(vect));
}

double ut_vector_max(const VectorDouble &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0.;
  double max = -1.e30;
  for (auto v : vec)
  {
    if (flagAbs) v = ABS(v);
    if (v > max) max = v;
  }
  return (max);
}

double ut_vector_min(const VectorDouble &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0.;
  double min = 1.e30;
  for (auto v : vec)
  {
    if (flagAbs) v = ABS(v);
    if (v < min) min = v;
  }
  return (min);
}

int ut_ivector_max(const VectorInt &vec)
{
  if (vec.size() <= 0) return 0;
  int max = -10000000;
  for (auto v : vec)
  {
    if (v > max) max = v;
  }
  return (max);
}

int ut_ivector_min(const VectorInt &vec)
{
  if (vec.size() <= 0) return 0;
  double min = 10000000;
  for (auto v : vec)
  {
    if (v < min) min = v;
  }
  return (min);
}

double ut_vector_mean(const VectorDouble &vec)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  int number = 0;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    mean += v;
    number++;
  }
  if (number > 0)
    mean /= (double) number;
  else
    mean = TEST;
  return (mean);
}

double ut_vector_cumul(const VectorDouble& vec)
{
  double total = 0.;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    total += v;
  }
  return total;
}

double ut_vector_var(const VectorDouble &vec)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  double var = 0.;
  int number = 0;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    var += v * v;
    mean += v;
    number++;
  }
  if (number > 0)
  {
    mean /= (double) number;
    var = var / (double) number - mean * mean;
  }
  else
  {
    mean = TEST;
    var = TEST;
  }
  return (var);
}

double ut_vector_stdv(const VectorDouble &vec)
{
  double var = ut_vector_var(vec);
  if (!FFFF(var))
    return (sqrt(var));
  else
    return TEST;
}

double ut_vector_inner_product(const VectorDouble &vec1,
                               const VectorDouble &vec2)
{
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  double prod = 0.;
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    prod += vec1.at(i) * vec2.at(i);
  return prod;
}

double ut_vector_norm(const VectorDouble &vec)
{
  double ip = static_cast<double>(ut_vector_inner_product(vec, vec));
  return sqrt(ip);
}

/**
 * Cross product (limited to 3D)
 * @param vec1 First vector
 * @param vec2 Second Vector
 * @return
 */
VectorDouble ut_vector_cross_product(const VectorDouble &vec1,
                                     const VectorDouble &vec2)
{
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  VectorDouble res;
  res.push_back(vec1[1] * vec2[2] - vec1[2] * vec2[1]);
  res.push_back(vec1[2] * vec2[0] - vec1[0] * vec2[2]);
  res.push_back(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
  return res;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (TEST if not defined)
 * @return
 */
bool ut_vector_constant(const VectorDouble& vect, double refval)
{
  if (vect.empty()) return false;
  if (FFFF(refval)) refval = vect[0];
  for (int i = 1; i < (int) vect.size(); i++)
    if (vect[i] != refval) return false;
  return true;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (ITEST if not defined)
 * @return
 */
bool ut_ivector_constant(const VectorInt& vect, int refval)
{
  if (vect.empty()) return false;
  if (IFFFF(refval)) refval = vect[0];
  for (int i = 1; i < (int) vect.size(); i++)
    if (vect[i] != refval) return false;
  return true;
}

bool ut_vector_same(const VectorDouble &v1, const VectorDouble &v2, double eps)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int>(v1.size()); i < n; i++)
    if (ABS(v1.at(i) - v2.at(i)) > eps) return false;
  return true;
}

bool ut_ivector_same(const VectorInt &v1, const VectorInt &v2)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int>(v1.size()); i < n; i++)
    if (ABS(v1.at(i) - v2.at(i)) > 0) return false;
  return true;
}

void ut_vector_fill(VectorDouble &vec, double value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void ut_ivector_fill(VectorInt &vec, int value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

VectorDouble ut_vector_concatenate(const VectorDouble& vec1,
                                   const VectorDouble& vec2)
{
  VectorDouble res = vec1;
  for (auto &e: vec2)
    res.push_back(e);
  return res;
}

VectorDouble ut_vector_add(const VectorDouble &vec1, const VectorDouble &vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    res.push_back(vec1.at(i) + vec2.at(i));
  return res;
}

/**
 * Performs: vec1 += vec2
 * @param dest Input/Output vector
 * @param src Auxiliary vector
 */
void ut_vector_add_inplace(VectorDouble &dest, const VectorDouble &src)
{
  VectorDouble res;
  if (dest.size() != src.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(dest.size()); i < n; i++)
    dest[i] += src[i];
}

/**
 * Return a vector containing vec2 - vec1
 * @param vec1 Input Vector
 * @param vec2 Input Vector
 * @return
 */
VectorDouble ut_vector_subtract(const VectorDouble &vec1,
                                const VectorDouble &vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    res.push_back(vec2.at(i) - vec1.at(i));
  return res;
}

VectorDouble ut_vector_power(const VectorDouble &vec, double power)
{
  int size = static_cast<int>(vec.size());
  VectorDouble res(size);
  for (int i = 0; i < size; i++)
    res[i] = pow(vec[i], power);
  return res;
}

VectorDouble ut_vector_simulate_uniform(int n, double mini, double maxi)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = law_uniform(mini, maxi);
  return vec;
}

VectorDouble ut_vector_simulate_bernoulli(int n, double proba, double vone, double velse)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
  {
    double rand = law_uniform(0., 1.);
    if (rand < proba)
      vec[i] = vone;
    else
      vec[i] = velse;
  }
  return vec;
}
/**
 * Performs: vec1 -= vec2
 * @param vec1 Input/Output vector
 * @param vec2 Auxiliary vector
 */
void ut_vector_subtract_inplace(VectorDouble &vec1, const VectorDouble &vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    vec1[i] -= vec2[i];
}

void ut_vector_sum(const VectorDouble &vec1,
                   const VectorDouble &vec2,
                   VectorDouble &res)
{
  if (vec1.size() != vec2.size())
  {
    my_throw("Wrong size");
  }
  int n = (int) vec1.size();
  if ((int) res.size() != n)
  {
    res.resize(n);
  }
  for (int i = 0; i < n; i++)
  {
    res[i] = vec1[i] + vec2[i];
  }
}

GSTLEARN_EXPORT VectorDouble ut_vector_simulate_gaussian(int n, double mean, double sigma)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = mean + sigma * law_gaussian();
  return vec;
}

// random generator function:
int myrandom (int i) { return std::rand()%i;}

/**
 * Sample a set of 'ntotal' consecutive ranks
 * @param ntotal      Dimension to be sampled
 * @param proportion  Proportion of elected samples (in [0,1])
 * @param number      Number of elected samples
 * @param seed        Seed used for the random number generator
 * @return A vector of indices lying between 0 and ntotal-1. No duplicate.
 *
 * @remark If 'proportion' and 'number' are not specified,
 * @remark the output vector has dimension equal to 'ntotal'
 */
VectorInt ut_vector_sample(int ntotal, double proportion, int number, int seed)
{
  // Find the number of expected values
  if (proportion <= 0. && number <= 0) return VectorInt();
  int count;
  if (proportion <= 0. && number <= 0)
    count = ntotal;
  else if (proportion > 0.)
    count = (int) (ntotal * proportion);
  else
    count = number;
  count = MIN(ntotal, MAX(1, count));

  // Define the Seed
  std::srand(seed);

  VectorInt ranks;
  for (int i = 0; i < ntotal; i++) ranks.push_back(i);

  std::random_shuffle ( ranks.begin(), ranks.end(), myrandom);

  ranks.resize(count);
  return ranks;
}

void ut_vector_cumul(VectorDouble &vec1, const VectorDouble &vec2, double coeff)
{
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    vec1[i] += coeff * vec2[i];
}

void ut_vector_copy(VectorDouble &vec1, const VectorDouble &vec2)
{
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(vec1.size()); i < n; i++)
    vec1[i] = vec2[i];
}

void ut_vector_multiply_inplace(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d *= v;});
}

void ut_vector_divide_inplace(VectorDouble &vec, double v)
{
  if (ABS(v) < EPSILON10)
  my_throw("division by 0");
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d /= v;});
}

void ut_vector_addval(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d += v;});
}

void ut_ivector_addval(VectorInt &vec, int v)
{
  std::for_each(vec.begin(), vec.end(), [v](int &d)
  { d += v;});
}

void ut_vector_divide_vec(VectorDouble &vec, const VectorDouble &v)
{
  if (vec.size() != v.size())
  my_throw("Arguments 'vec' and 'v' should have same dimension");
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    if (ABS(v[i]) < EPSILON20)
      my_throw("division by 0");
    vec[i] /= v[i];
  }
}

int ut_vector_count_undefined(const VectorDouble &vec)
{
  int count = 0;
  for (int i = 0; i < (int) vec.size(); i++)
  {
    if (FFFF(vec[i])) count++;
  }
  return count;
}

int ut_vector_prod(const VectorInt& vec)
{
  if (vec.empty()) return 0;
  int nprod = 1;
  for (int i = 0; i < (int) vec.size(); i++)
    nprod *= vec[i];
  return nprod;
}

double ut_vector_prod(const VectorDouble& vec)
{
  if (vec.empty()) return 0;
  double nprod = 1.;
  for (int i = 0; i < (int) vec.size(); i++)
    nprod *= vec[i];
  return nprod;
}

/**
 * Create an output vector containing the 'number' consecutive numbers starting from 'ideb'
 *
 * @param number  Length of the output vector
 * @param ideb    Index of the first element of the output vector
 */
VectorInt ut_ivector_sequence(int number, int ideb)
{
  VectorInt vec(number);

  for (int i = 0; i < number; i++)
    vec[i] = ideb + i;
  return vec;
}

/**
 * Create an output vector going from 'valFrom' to 'ValTo' by step of 'valStep'
 */
VectorDouble ut_vector_sequence(double valFrom, double valTo, double valStep)
{
  VectorDouble vec;

  double value = valFrom;
  while (value <= valTo)
  {
    vec.push_back(value);
    value = value + valStep;
  }
  return vec;
}

int ut_vector_size(const VectorInt &vec)
{
  if (vec.empty()) return 0;
  int size = sizeof(VectorInt) + (sizeof(int) * (int) vec.size());
  return size;
}

int ut_vector_size(const VectorDouble &vec)
{
  if (vec.empty()) return 0;
  int size = sizeof(VectorDouble) + (sizeof(double) * (int) vec.size());
  return size;
}

int ut_vector_size(const VectorVectorInt &vec)
{
  int size = 0;
  if (vec.empty()) return size;
  for (auto i = 0; i != (int) vec.size(); i++)
    size += sizeof(VectorInt) + (sizeof(int) * (int) vec[i].size());
  return size;

}

int ut_vector_size(const VectorVectorDouble &vec)
{
  int size = 0;
  if (vec.empty()) return size;
  for (auto i = 0; i != (int) vec.size(); i++)
    size += sizeof(VectorDouble) + (sizeof(double) * (int) vec[i].size());
  return size;
}

VectorInt ut_ivector_set(int* values, int number)
{
  if (values == nullptr) return VectorInt();
  VectorInt vec(number);
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorDouble ut_vector_set(double* values, int number)
{
  if (values == nullptr) return VectorDouble();
  VectorDouble vec;
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorVectorDouble ut_vector_vector_set(double* value, int n1, int n2)
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

VectorInt ut_ivector_sort(const VectorInt& vecin, bool ascending)
{
  if (vecin.empty()) return VectorInt();

  VectorInt vecout = vecin;
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

VectorDouble ut_vector_sort(const VectorDouble& vecin, bool ascending)
{
  if (vecin.empty()) return VectorDouble();

  VectorDouble vecout = vecin;
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

VectorInt ut_vector_sort_indices(const VectorDouble& vecin)
{
  if (vecin.empty()) return VectorInt();

  VectorInt idx(vecin.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&vecin](size_t i1, size_t i2) {return vecin[i1] < vecin[i2];});

  return idx;
}

/**
 * Calculate the diagonal of the box extension
 * @param mini Array of lower coordinates of the box
 * @param maxi Array of upper coordinates of the box
 * @return
 * @remark If one coordinate is undefined, TEST is returned.
 */
double ut_vector_extension_diagonal(const VectorDouble& mini,
                                    const VectorDouble& maxi)
{
  double diag = 0.;
  VectorDouble delta = ut_vector_subtract(mini, maxi);
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

VectorDouble ut_vector_angle_from_codir(const VectorDouble& codir)
{
  int ndim = (int) codir.size();
  VectorDouble angles(ndim);
  ut_angles_from_codir(ndim, codir, angles);
  return angles;
}

