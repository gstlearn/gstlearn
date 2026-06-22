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
#include "Basic/VectorHelper.hpp"
#include "geoslib_f.h"
#include <chrono>

using namespace gstlrn;

/**
 * Benchmark for VectorHelper::whereElement optimization
 * This demonstrates the performance improvement for large datasets
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  message("\n=== Benchmark: whereElement Performance ===\n");

  // Create a large sorted vector (simulating index1[ivar] from InvNuggetOp)
  Id n = 10000; // Size of the index vector
  VectorInt indices;
  indices.reserve(n);
  for (Id i = 0; i < n; i++)
    indices.push_back(i * 2); // Sorted, with gaps (simulating sparse data)

  // Create search targets (simulating iech values in the loop)
  VectorInt targets;
  targets.reserve(n);
  for (Id i = 0; i < n; i++) targets.push_back(i * 2);

  message("Vector size: %d\n", n);
  message("Number of searches: %d\n", n);
  message("\n");

  // Benchmark 1: Original implementation (searching from beginning each time)
  {
    auto start = std::chrono::high_resolution_clock::now();

    Id count = 0;
    for (Id i = 0; i < n; i++)
    {
      Id pos = VH::whereElement(indices, targets[i]);
      if (pos >= 0) count++;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    message("Original (start=0 always):\n");
    message("  Found %d/%d elements\n", count, n);
    messageNoDiff("  Time: %lld microseconds\n", duration.count());
    message("\n");
  }

  // Benchmark 2: Optimized implementation (searching from last position)
  {
    auto start = std::chrono::high_resolution_clock::now();

    Id count = 0;
    Id lastPos = 0;
    for (Id i = 0; i < n; i++)
    {
      Id pos = VH::whereElement(indices, targets[i], lastPos);
      if (pos >= 0)
      {
        count++;
        lastPos = pos;
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    message("Optimized (start=lastPos):\n");
    message("  Found %d/%d elements\n", count, n);
    messageNoDiff("  Time: %lld microseconds\n", duration.count());
    message("\n");
  }

  message("=== Benchmark Complete ===\n");
  message(
    "Note: The optimized version reduces complexity from O(n²) to O(n)\n");
  message(
    "      For larger datasets (millions of points), the improvement is "
    "dramatic.\n");

  return 0;
}
