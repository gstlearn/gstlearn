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
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixSquare.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This file is means to check the performance of several programming rules
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  Timer timer;
  double result     = 0.;
  double result_ref = 0.;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  law_set_random_seed(1331742);

  // Assigning values in a matrix

  Id nx             = 10000;
  Id naffect        = 5000000;
  VectorDouble vecS = VH::simulateUniform(nx * nx);
  MatrixSquare matS(nx);
  matS.resetFromVD(nx, nx, vecS);

  mestitle(1, "Assigning values into a storage");
  message("Random values are assigned at random locations %d times\n", naffect);
  double bidon = 13431.;

  message("- Assigning value to the vector of dimension %d\n", nx * nx);
  timer.reset();
  for (Id itime = 0; itime < naffect; itime++)
  {
    Id rank    = law_uniform(0, nx * nx);
    vecS[rank] = bidon;
  }
  timer.displayIntervalMilliseconds("Assignment to vector", 120);

  message("- Assigning value to the square matrix of dimension %d x %d\n", nx, nx);
  timer.reset();
  for (Id itime = 0; itime < naffect; itime++)
  {
    Id rank = law_uniform(0, nx);
    matS.setValue(rank, rank, bidon);
  }
  timer.displayIntervalMilliseconds("Assignment to matrix", 75);

  // Comparing several uses of VectorDouble calculations
  // ===================================================
  Id ntimes = 5000;
  Id nsize  = 30000;

  VectorDouble a = VH::simulateGaussian(nsize);
  VectorDouble b = VH::simulateGaussian(nsize);

  mestitle(1, "Comparing various ways of operating Vectors of Double values");
  message("Operations are performed %d times over vectors of size %d\n", ntimes, nsize);

  message("- using: sum_i a[i] * b[i]\n");
  timer.reset();

  for (Id itime = 0; itime < ntimes; itime++)
  {
    result = 0.;
    for (Id ielem = 0; ielem < nsize; ielem++)
    {
      result += a[ielem] * b[ielem];
    }
  }
  timer.displayIntervalMilliseconds("Product of [] terms", 500);
  result_ref = result;

  message("- using iterators\n");
  timer.reset();
  for (Id itime = 0; itime < ntimes; itime++)
  {
    VectorDouble::const_iterator ita(a.begin());
    VectorDouble::const_iterator itb(b.begin());
    result = 0.;
    while (ita < a.end())
    {
      result += (*ita) * (*itb);
      ita++;
      itb++;
    }
  }
  timer.displayIntervalMilliseconds("with iterators", 450);
  if (ABS(result - result_ref) > ABS(result_ref) * EPSILON6)
    message("Results are different: Result = %lf; Ref = %lf\n", result, result_ref);

  message("- using pointers to double\n");
  timer.reset();
  for (Id itime = 0; itime < ntimes; itime++)
  {
    double* pta = a.data();
    double* ptb = b.data();
    result      = 0.;
    for (Id ielem = 0; ielem < nsize; ielem++)
    {
      result += (*pta) * (*ptb);
      pta++;
      ptb++;
    }
  }
  timer.displayIntervalMilliseconds("with pointers", 320);
  if (ABS(result - result_ref) > ABS(result_ref) * EPSILON6)
    message("Results are different: Result = %lf; Ref = %lf\n", result, result_ref);

  message("- using VectorHelper\n");
  timer.reset();
  for (Id itime = 0; itime < ntimes; itime++)
    result = a.innerProduct(b);
  timer.displayIntervalMilliseconds("with VectorHelper", 200);
  if (ABS(result - result_ref) > ABS(result_ref) * EPSILON6)
    message("Results are different: Result = %lf; Ref = %lf\n", result, result_ref);

  message("- using VectorHelper (double)\n");
  timer.reset();
  const double* ptra = a.data();
  const double* ptrb = b.data();
  for (Id itime = 0; itime < ntimes; itime++)
    result = VH::innerProduct(ptra, ptrb, nsize);
  timer.displayIntervalMilliseconds("with VectorHelper (double)", 200);
  if (ABS(result - result_ref) > ABS(result_ref) * EPSILON6)
    message("Results are different: Result = %lf; Ref = %lf\n", result, result_ref);

  message("- using matrix algebra\n");
  MatrixDense mata;
  mata.resetFromVD(1, nsize, a);
  MatrixDense matb;
  matb.resetFromVD(nsize, 1, b);
  MatrixDense res(1, 1);
  timer.reset();
  for (Id itime = 0; itime < ntimes; itime++)
  {
    res.prodMatMatInPlace(&mata, &matb);
    result = res(0, 0);
  }
  timer.displayIntervalMilliseconds("with algebra", 1700);
  if (ABS(result - result_ref) > ABS(result_ref) * EPSILON6)
    message("Results are different: Result = %lf; Ref = %lf\n", result, result_ref);

  /// Sorting the contents of a vector

  mestitle(1, "Testing sorting algorithms");
  Id nech = 10;
  Id size = 7;
  message("We consider a vector of %d values and the corresponding vector of ranks\n", nech);
  message("Only the first %d positions are used\n", size);
  message("This paragraph is not bench-marked as time consumption is too short\n");

  VectorDouble VinVal = VH::simulateUniform(nech);
  VectorInt VinRank   = VH::sequence(nech, 4, 3);
  printVector(VinVal, "Unsorted values", true, true);
  printVector(VinRank, "Unsorted ranks", true, true);

  VectorInt order = VH::orderRanks(VinVal, true, size);
  printVector(order, "Order", true, true);
  VectorDouble VoutVal = VH::sort(VinVal, true, size);
  printVector(VoutVal, "Sorted values", true, true);

  VectorDouble VsortVal = VH::reorder(VinVal, order, size);
  if (!VoutVal.isEqual(VsortVal))
    printVector(VsortVal, "Results are different: Re-ordered values", true, true);

  VectorInt VsortRank = VH::reorder(VinRank, order, size);
  printVector(VsortRank, "Ranks of Sorted values", true, true);

  VH::arrangeInPlace(0, VinRank, VinVal, true, size);
  VinVal.resize(size);
  if (!VoutVal.isEqual(VinVal))
    printVector(VinVal, "Results are different: Re-arranged values", true, true);
  VinRank.resize(size);
  if (!VsortRank.isEqual(VinRank))
    printVector(VinRank, "Re-arranged ranks", true, true);
  return (0);
}
