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
#include "Basic/ASerializable.hpp"
#include "geoslib_define.h"

#include "Transform/YeoJohnson.hpp"
#include "Transform/YeoJohnsonForTest.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  {
    YeoJohnson transform(0.5);
    double x = 3.0;
    double y = transform.inverseTransform(x);
    double z = transform.transform(y);
    std::cout << "Original x: " << x << "\n";
    std::cout << "Inverse Transformed y: " << y << "\n";
    std::cout << "Transformed z: " << z << "\n";
  }

  {
    YeoJohnson transform(0.);
    double x = exp(2.0) - 1;
    double y = transform.inverseTransform(x);
    double z = transform.transform(y);
    std::cout << "Original x: " << x << "\n";
    std::cout << "Inverse Transformed y: " << y << "\n";
    std::cout << "Transformed z: " << z << "\n";
  }

  {
    YeoJohnson transform(0.5);
    double x = -2.0;
    double y = transform.inverseTransform(x);
    double z = transform.transform(y);
    std::cout << "Original x: " << x << "\n";
    std::cout << "Inverse Transformed y: " << y << "\n";
    std::cout << "Transformed z: " << z << "\n";
  }

  {
    YeoJohnson transform(2);
    double x = -exp(2.0) + 1;
    double y = transform.inverseTransform(x);
    double z = transform.transform(y);
    std::cout << "Original x: " << x << "\n";
    std::cout << "Inverse Transformed y: " << y << "\n";
    std::cout << "Transformed z: " << z << "\n";
  }

  YeoJohnsonForTest transformTest(0.5);
  YeoJohnson transformRef(0.5);

  double test_points[] = {-3.0, -1.0, 0.0, 1.0, 2.0, 3.0};
  double lambda[] = {0.0, 0.5, 1.0, 2.0};

  for (double l: lambda)
  {
    transformTest.setLambda(l);
    transformRef.setLambdaValue(l);
    std::cout << "Testing for lambda = " << l << "\n";
    for (double x: test_points)
    {

      double y_test = transformTest.inverseTransform(x);
      double y_ref = transformRef.inverseTransform(x);
      std::cout << "- Test point: " << x << "\n";
      std::cout << "  -Diff inverse: " << std::abs(y_test - y_ref) << "\n";
      double jacobian_test = transformTest.evalJacobian(x);
      double jacobian_ref = transformRef.evalJacobian(x);
      std::cout << "  -Diff Jacobian: "
                << std::abs(jacobian_test - jacobian_ref) << "\n";
    }
  }
  return 0;
}
