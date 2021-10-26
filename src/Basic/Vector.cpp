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
#include "math.h"
#include "geoslib_e.h"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include <iomanip>

String ut_vector_string(const VectorDouble& vec)
{
  return toVector(String(),vec);
}

String ut_ivector_string(const VectorInt& vec)
{
  return toVector(String(), vec);
}

void ut_vector_display(const String& title, const VectorDouble& vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_vector_string(vect));
}

void ut_vector_display_stats(const String& title, const VectorDouble& vect)
{
  int ntotal  = (int) vect.size();
  int number  = 0;
  double mean = 0.;
  double stdv = 0.;
  double mini =  1.e30;
  double maxi = -1.e30;

  for (int i=0; i<ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    mean += value;
    stdv += value * value;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (! title.empty()) message("%s\n",title.c_str());
  if (number > 0)
  {
    mean /= (double) number;
    stdv  = stdv / (double) number - mean * mean;
    stdv  = (stdv > 0.) ? sqrt(stdv) : 0.;

    message("- Number of samples = %d / %d\n",number,ntotal);
    message("- Minimum  = %lf\n",mini);
    message("- Maximum  = %lf\n",maxi);
    message("- Mean     = %lf\n",mean);
    message("- St. Dev. = %lf\n",stdv);
  }
  else
  {
    message("No value defined\n");
  }
}

void ut_vector_display_range(const String& title, const VectorDouble& vect)
{
  int ntotal  = (int) vect.size();
  int number  = 0;
  double mini =  1.e30;
  double maxi = -1.e30;

  for (int i=0; i<ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (! title.empty()) message("%s\n",title.c_str());
  if (number > 0)
  {
    message("- Number of samples = %d / %d\n",number,ntotal);
    message("- Minimum  = %lf\n",mini);
    message("- Maximum  = %lf\n",maxi);
  }
  else
  {
    message("No value defined\n");
  }
}

void ut_ivector_display(const String& title,const VectorInt& vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(ut_ivector_string(vect));
}

double ut_vector_max(const VectorDouble& vec, bool flagAbs)
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

double ut_vector_min(const VectorDouble& vec, bool flagAbs)
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

double ut_vector_mean(const VectorDouble& vec)
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

double ut_vector_var(const VectorDouble& vec)
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
    var   = var / (double) number - mean * mean;
  }
  else
  {
    mean = TEST;
    var  = TEST;
  }
  return (var);
}

double ut_vector_stdv(const VectorDouble& vec)
{
  double var = ut_vector_var(vec);
  if (! FFFF(var))
    return (sqrt(var));
  else
    return TEST;
}

double ut_vector_inner_product(const VectorDouble& vec1,
                               const VectorDouble& vec2)
{
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  double prod = 0.;
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    prod += vec1.at(i) * vec2.at(i);
  return prod;
}

double ut_vector_norm(const VectorDouble& vec)
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
VectorDouble ut_vector_cross_product(const VectorDouble& vec1,
                                                 const VectorDouble& vec2)
{
  if (vec1.size() != vec2.size())
  my_throw("Wrong size");
  VectorDouble res;
  res.push_back(vec1[1] * vec2[2] - vec1[2] * vec2[1]);
  res.push_back(vec1[2] * vec2[0] - vec1[0] * vec2[2]);
  res.push_back(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
  return res;
}

bool ut_vector_same(const VectorDouble& v1, const VectorDouble& v2, double eps)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int> (v1.size()); i < n; i++)
    if (abs(v1.at(i) - v2.at(i)) > eps) return false;
  return true;
}

bool ut_ivector_same(const VectorInt& v1, const VectorInt& v2)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int> (v1.size()); i < n; i++)
    if (abs(v1.at(i) - v2.at(i)) > 0) return false;
  return true;
}

void ut_vector_fill(VectorDouble& vec, double value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void ut_ivector_fill(VectorInt& vec, int value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

VectorDouble ut_vector_add(const VectorDouble& vec1, const VectorDouble& vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    res.push_back(vec1.at(i) + vec2.at(i));
  return res;
}

/**
 * Performs: vec1 += vec2
 * @param vec1 Input/Output vector
 * @param vec2 Auxiliary vector
 */
void ut_vector_add_inplace(VectorDouble& vec1, const VectorDouble& vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    vec1[i] += vec2[i];
}

VectorDouble ut_vector_subtract(const VectorDouble& vec1, const VectorDouble& vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    res.push_back(vec2.at(i) - vec1.at(i));
  return res;
}

VectorDouble ut_vector_power(const VectorDouble&vec, double power)
{
  int size = static_cast<int> (vec.size());
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
/**
 * Performs: vec1 -= vec2
 * @param vec1 Input/Output vector
 * @param vec2 Auxiliary vector
 */
void ut_vector_subtract_inplace(VectorDouble& vec1, const VectorDouble& vec2)
{
  VectorDouble res;
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    vec1[i] -= vec2[i];
}

void ut_vector_sum(const VectorDouble& vec1,const VectorDouble& vec2,VectorDouble& res )
{
  if (vec1.size() != vec2.size())
  {
    my_throw("Wrong size");
  }
  int n = (int)vec1.size();
  if ((int)res.size() != n)
  {
    res.resize(n);
  }
  for (int i = 0; i< n; i++)
  {
    res[i] = vec1[i] + vec2[i];
  }
}


VectorDouble ut_vector_simulate_gaussian(int n, double mean, double sigma)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = mean + sigma * law_gaussian();
  return vec;
}

/**
 * Sample a set of 'ntotal' consecutive ranks
 * @param ntotal      Dimension to be sampled
 * @param proportion  Proportion of elected samples (in [0,1])
 * @param seed        Seed used for the random number generator
 * @return A vector of ranks (between 0 and 'ntotal-1'). No duplicate
 */
VectorInt ut_vector_sample(int ntotal, double proportion, int seed)
{
  int number = (int) (ntotal * proportion);
  number = MIN(ntotal, MAX(1, number));
  VectorInt ranks(number);

  int memo = law_get_random_seed();
  if (seed > 0) law_set_random_seed(seed);
  for (int i = 0; i < number; i++)
    ranks[i] = (int) law_uniform(0.,(double) ntotal);

  // Check that the ranks are 'unique'
  std::sort(ranks.begin(), ranks.end());
  auto last = std::unique(ranks.begin(), ranks.end());
  ranks.erase(last, ranks.end());

  law_set_random_seed(memo);
  return ranks;
}

void ut_vector_cumul(VectorDouble& vec1, const VectorDouble& vec2, double coeff)
{
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    vec1[i] += coeff * vec2[i];
}

void ut_vector_copy(VectorDouble& vec1, const VectorDouble& vec2)
{
  if (vec1.size() != vec2.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int> (vec1.size()); i < n; i++)
    vec1[i] = vec2[i];
}

void ut_vector_multiply_inplace(VectorDouble& vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double& d) { d *= v;});
}

void ut_vector_divide_inplace(VectorDouble& vec, double v)
{
  if (ABS(v) < EPSILON10)
    my_throw("division by 0");
  std::for_each(vec.begin(), vec.end(), [v](double& d) { d /= v;});
}

void ut_vector_addval(VectorDouble& vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double& d) { d += v;});
}

void ut_ivector_addval(VectorInt& vec, int v)
{
  std::for_each(vec.begin(), vec.end(), [v](int& d) { d += v;});
}

void ut_vector_divide_vec(VectorDouble& vec, const VectorDouble& v)
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

int ut_vector_count_undefined(const VectorDouble& vec)
{
  int count = 0;
  for (int i = 0; i< (int) vec.size(); i++)
  {
    if (FFFF(vec[i])) count++;
  }
  return count;
}

int ut_ivector_prod(const VectorInt vec)
{
  if (vec.empty()) return 0;
  int nprod = 1;
  for (int i = 0; i < (int) vec.size(); i++) nprod *= vec[i];
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

