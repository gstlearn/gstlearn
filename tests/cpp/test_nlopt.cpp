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

/**
 * This test aims to test the use of the nlopt library for optimization. 
 * Only a quadratic function is minimized in this example.
 */
#include <iostream>
#include "Basic/File.hpp"
#include "geoslib_define.h"
#include <sstream>
#include <nlopt.h>

// Function to minimize
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data = nullptr)
{
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(my_func_data);
    if (grad) {
        grad[0] = 2 * (x[0] - 3);
    }
    return (x[0] - 3) * (x[0] - 3);
}

int main(int argc, char *argv[])
{
    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);
    nlopt_opt opt= nlopt_create(NLOPT_LD_LBFGS, 2);

    // Bounds for each parameter
    nlopt_set_lower_bound(opt, 0, 1);
    nlopt_set_lower_bound(opt, 1, -10);
    nlopt_set_upper_bound(opt, 0, 5);
    nlopt_set_upper_bound(opt, 1, 10);


    nlopt_set_min_objective(opt, myfunc, nullptr);

    // Set the tolerance for the stopping criteria
    nlopt_set_xtol_rel(opt,EPSILON4);

    // Starting point
    double x[2] = {1., 5.0}; 

    double minf; // minimum value of the objective function
    try {
        nlopt_optimize(opt,x, &minf);
        std::cout << "End optimization" << std::endl;
        std::cout << "x = " << x[0] << " " << x[1] << std::endl;
        std::cout << "Minimum value of the objective " << minf << std::endl;
    }
    catch (std::exception &e) {
        std::cerr << "Error in the optimization: " << e.what() << std::endl;
    }

    return 0;
}