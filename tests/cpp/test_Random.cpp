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
#include "Basic/Law.hpp"
#include "Basic/AStringable.hpp"

/**
 * This test is meant to check the generation of random values
 */

void st_do_it(bool style, int seed)
{
  int number = 3;

  if (style)
    message("\nUsing Old Style Random Number Generator\n");
  else
    message("\nUsing New Style Random Number Generator\n");

  law_set_old_style(style);

  // Setting the seed
  law_set_random_seed(seed);
  message("Getting the seed = %d\n", law_get_random_seed());

  // Sampling Uniform distribution
  for (int i = 0; i < number; i++)
    message("Uniform=%lf\n",law_uniform());

  // Sampling Gaussian distribution
  for (int i = 0; i < number; i++)
    message("Gaussian=%lf\n",law_gaussian());

  // Sampling Exponential distribution
  for (int i = 0; i < number; i++)
    message("Exponential=%lf\n",law_exponential());

  // Sampling Gamma distribution
  for (int i = 0; i < number; i++)
    message("Gamma=%lf\n",law_gamma(3.4,1.));

  // Sampling Poisson distribution
  for (int i = 0; i < number; i++)
    message("Poisson=%d\n",law_poisson(14.6));

}

int main()
{
  int seed = 432432;

  st_do_it(false, seed);

  st_do_it(true, seed);

  return 0;
}
