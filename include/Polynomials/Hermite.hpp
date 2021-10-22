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
#pragma once

#include "Basic/Vector.hpp"
#include "Matrix/MatrixSGeneral.hpp"

VectorDouble hermitePolynomials(double yc, double r, int nbpoly);
VectorDouble hermiteCoefIndicator(double yc, int nbpoly);
VectorDouble hermiteCoefMetal(double yc, const VectorDouble& phi);
VectorDouble hermiteFunction(double y, int nbpoly);
MatrixSGeneral hermiteIncompleteIntegral(double yc, int nbpoly);
VectorDouble hermiteLognormal(double mean, double sigma, int nbpoly);
double hermiteSeries(const VectorDouble& an, const VectorDouble& hn);

VectorDouble hermiteIndicator(double yc,
                              VectorDouble krigest,
                              VectorDouble krigstd);
double hermiteIndicatorElement(double yc, double krigest, double krigstd);
VectorDouble hermiteIndicatorStd(double yc,
                                 VectorDouble krigest,
                                 VectorDouble krigstd);
double hermiteIndicatorStdElement(double yc, double krigest, double krigstd);
VectorDouble hermiteMetal(double yc,
                          VectorDouble krigest,
                          VectorDouble krigstd,
                          const VectorDouble& phi);
double hermiteMetalElement(double yc,
                           double krigest,
                           double krigstd,
                           const VectorDouble& phi);
VectorDouble hermiteMetalStd(double yc,
                             VectorDouble krigest,
                             VectorDouble krigstd,
                             const VectorDouble& phi);
double hermiteMetalStdElement(double yc,
                              double krigest,
                              double krigstd,
                              const VectorDouble& phi);
VectorDouble hermiteCondExp(VectorDouble krigest,
                            VectorDouble krigstd,
                            const VectorDouble& phi);
double hermiteCondExpElement(double krigest,
                             double krigstd,
                             const VectorDouble& phi);
VectorDouble hermiteCondStd(VectorDouble krigest,
                            VectorDouble krigstd,
                            const VectorDouble& phi);
double hermiteCondStdElement(double krigest,
                             double krigstd,
                             const VectorDouble& phi);
VectorDouble hermiteEvaluateZ2(VectorDouble yk,
                               VectorDouble sk,
                               const VectorDouble& phi);
double hermiteEvaluateZ2(double yk, double sk, const VectorDouble& phi);
