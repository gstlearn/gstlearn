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
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"

/**
 * This test is meant to check the generation of random values
 */

void st_do_it(bool style, int seed)
{
  if (style)
    message("\nUsing Old Style Random Number Generator\n");
  else
    message("\nUsing New Style Random Number Generator\n");

  law_set_old_style(style);

  // Setting the seed
  law_set_random_seed(seed);
  message("Getting the seed = %d\n", law_get_random_seed());
}

int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int seed = 432432;
  int number = 10000;
  VectorDouble tab;
  message("All subsequent tests are performed on a set of %d samples\n",number);

  // Sampling Uniform Distribution

  double mini = 2.;
  double maxi = 4.;
  tab.resize(number,0.);
  message("\nUniform Distribution: mini=%lf maxi=%lf\n",mini,maxi);
  law_set_old_style(true);
  for (int i = 0; i < number; i++) tab[i] = law_uniform(mini,maxi);
  VH::displayStats("Old Style", tab);
  law_set_old_style(false);
  for (int i = 0; i < number; i++) tab[i] = law_uniform(mini,maxi);
  VH::displayStats("New Style", tab);

  // Sampling Gaussian distribution

  double mean = 5.;
  double stdv = 2.;
  tab.resize(number,0.);
  message("\nNormal Distribution: mean=%lf sigma=%lf\n",mean,stdv);
  law_set_old_style(true);
  for (int i = 0; i < number; i++) tab[i] = law_gaussian();
  VH::displayStats("Old Style", tab);
  law_set_old_style(false);
  for (int i = 0; i < number; i++) tab[i] = law_gaussian();
  VH::displayStats("New Style", tab);

  // Sampling Exponential distribution

  double param = 2.;
  tab.resize(number,0.);
  message("\nExponential Distribution: param=%lf\n",param);
  law_set_old_style(true);
  for (int i = 0; i < number; i++) tab[i] = law_exponential(param);
  VH::displayStats("Old Style", tab);
  law_set_old_style(false);
  for (int i = 0; i < number; i++) tab[i] = law_exponential(param);
  VH::displayStats("New Style", tab);

  // Sampling Gamma distribution

  double alpha = 3.4;
  double beta  = 1.2;
  tab.resize(number,0.);
  message("\nGamma Distribution: alpha=%lf beta=%lf\n",alpha,beta);
  law_set_old_style(true);
  for (int i = 0; i < number; i++) tab[i] = law_gamma(alpha,beta);
  VH::displayStats("Old Style", tab);
  law_set_old_style(false);
  for (int i = 0; i < number; i++) tab[i] = law_gamma(alpha,beta);
  VH::displayStats("New Style", tab);

  // Sampling Poisson distribution

  double lambda = 4.6;
  number = 100000; // Use larger sample for convergence
  tab.resize(number,0.);
  int ndec = (int) OptCst::query(ECst::NTDEC);
  OptCst::define(ECst::NTDEC, 4);
  message("\nPoisson Distribution: lambda=%lf\n",lambda);
  law_set_old_style(true);
  for (int i = 0; i < number; i++) tab[i] = (double) law_poisson(lambda);
  VH::displayStats("Old Style", tab);
  law_set_old_style(false);
  for (int i = 0; i < number; i++) tab[i] = (double) law_poisson(lambda);
  VH::displayStats("New Style", tab);
  OptCst::define(ECst::NTDEC, ndec);

  // Retrieve seed

  st_do_it(false, seed);
  st_do_it(true, seed);

  // Tossing few Uniform values in both styles

  mini = 2.;
  maxi = 4.;
  number = 5;
  tab.resize(number,0.);
  message("\nUniform Distribution: mini=%lf maxi=%lf\n",mini,maxi);
  law_set_old_style(true);
  law_set_random_seed(seed);
  for (int i = 0; i < number; i++) tab[i] = law_uniform(mini,maxi);
  VH::display("Old Style", tab);
  law_set_old_style(false);
  law_set_random_seed(seed);
  for (int i = 0; i < number; i++) tab[i] = law_uniform(mini,maxi);
  VH::display("New Style", tab);

  return 0;
}
