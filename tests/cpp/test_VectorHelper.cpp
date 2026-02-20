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

using namespace gstlrn;

/**
 * Test for VectorHelper::whereElement optimization
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Test 1: Basic whereElement functionality (without start parameter)
  message("\n=== Test 1: Basic whereElement ===\n");
  VectorInt vec1 = {10, 20, 30, 40, 50};
  Id pos         = VH::whereElement(vec1, 30);
  if (pos == 2)
    message("PASS: Found 30 at position 2\n");
  else
    message("FAIL: Expected position 2, got %d\n", pos);

  pos = VH::whereElement(vec1, 100);
  if (pos == -1)
    message("PASS: Element 100 not found (returned -1)\n");
  else
    message("FAIL: Expected -1 for missing element, got %d\n", pos);

  // Test 2: Optimized whereElement with start parameter
  message("\n=== Test 2: Optimized whereElement with start ===\n");
  VectorInt vec2 = {5, 15, 25, 35, 45, 55, 65, 75, 85, 95};

  // Search starting from position 0
  pos = VH::whereElement(vec2, 25, 0);
  if (pos == 2)
    message("PASS: Found 25 at position 2 (start=0)\n");
  else
    message("FAIL: Expected position 2, got %d\n", pos);

  // Search starting from position 2 (should still find it)
  pos = VH::whereElement(vec2, 25, 2);
  if (pos == 2)
    message("PASS: Found 25 at position 2 (start=2)\n");
  else
    message("FAIL: Expected position 2, got %d\n", pos);

  // Search starting from position 3 (should not find it as we only search forward)
  pos = VH::whereElement(vec2, 25, 3);
  if (pos == -1)
    message("PASS: Element 25 not found when starting from position 3 (returned -1)\n");
  else
    message("FAIL: Expected -1 when element is before start, got %d\n", pos);

  // Test 3: Sequential search pattern (simulating InvNuggetOp use case)
  message("\n=== Test 3: Sequential search pattern ===\n");
  VectorInt indices = {0, 2, 5, 7, 10, 12, 15, 18, 20, 25};
  VectorInt targets = {0, 2, 5, 7, 10, 12, 15, 18, 20, 25};

  Id lastPos    = 0;
  bool allFound = true;
  for (Id i = 0; i < static_cast<Id>(targets.size()); i++)
  {
    pos = VH::whereElement(indices, targets[i], lastPos);
    if (pos >= 0)
    {
      message("Found %d at position %d (started from %d)\n", targets[i], pos, lastPos);
      lastPos = pos; // Update for next search
    }
    else
    {
      message("FAIL: Could not find %d\n", targets[i]);
      allFound = false;
    }
  }

  if (allFound)
    message("PASS: All elements found in sequential order\n");
  else
    message("FAIL: Some elements not found\n");

  // Test 4: Edge cases
  message("\n=== Test 4: Edge cases ===\n");

  // Empty vector
  VectorInt emptyVec;
  pos = VH::whereElement(emptyVec, 10);
  if (pos == -1)
    message("PASS: Empty vector returns -1\n");
  else
    message("FAIL: Expected -1 for empty vector, got %d\n", pos);

  pos = VH::whereElement(emptyVec, 10, 0);
  if (pos == -1)
    message("PASS: Empty vector with start returns -1\n");
  else
    message("FAIL: Expected -1 for empty vector with start, got %d\n", pos);

  // Single element
  VectorInt singleVec = {42};
  pos                 = VH::whereElement(singleVec, 42, 0);
  if (pos == 0)
    message("PASS: Single element found at position 0\n");
  else
    message("FAIL: Expected position 0 for single element, got %d\n", pos);

  // Start beyond vector size
  VectorInt vec3 = {1, 2, 3};
  pos            = VH::whereElement(vec3, 2, 10);
  if (pos == -1)
    message("PASS: Start beyond vector size returns -1\n");
  else
    message("FAIL: Expected -1 when start > size, got %d\n", pos);

  // New interface
  mestitle(0, "Testing new VectorHelper interface for addition");
  auto nech = 5;
  auto V1   = VH::simulateGaussian(nech, 0., 1.);
  auto V2   = VH::simulateGaussian(nech, 0., 1.);
  VectorDouble V4;

  VectorDouble V3 = VH::add(V1, V2);
  V3.dump("Checking VH::add(V1,V2)");

  VH::add(V4, V1, V2);
  V4.dump("Checking VH::add(v4,V1,V2)");

  VectorDouble V5 = V1 + V2;
  V5.dump("Checking V = V1 + V2");

  VectorDouble V6 = V1;
  V6 += V2;
  V6.dump("Checking V6 += V2");

  message("\n=== All tests completed ===\n");
  return 0;
}
