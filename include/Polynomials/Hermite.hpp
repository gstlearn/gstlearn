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

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

GSTLEARN_EXPORT VectorDouble hermitePolynomials(double yc,
                                                double r,
                                                int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteCoefIndicator(double yc, int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteCoefMetal(double yc,
                                              const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteFunction(double y, int nbpoly);
GSTLEARN_EXPORT MatrixSquareGeneral hermiteIncompleteIntegral(double yc,
                                                              int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteLognormal(double mean,
                                              double sigma,
                                              int nbpoly);
GSTLEARN_EXPORT double hermiteSeries(const VectorDouble &an,
                                     const VectorDouble &hn);

GSTLEARN_EXPORT VectorDouble hermiteIndicator(double yc,
                                              VectorDouble krigest,
                                              VectorDouble krigstd);
GSTLEARN_EXPORT double hermiteIndicatorElement(double yc,
                                               double krigest,
                                               double krigstd);
GSTLEARN_EXPORT VectorDouble hermiteIndicatorStd(double yc,
                                                 VectorDouble krigest,
                                                 VectorDouble krigstd);
GSTLEARN_EXPORT double hermiteIndicatorStdElement(double yc,
                                                  double krigest,
                                                  double krigstd);
GSTLEARN_EXPORT VectorDouble hermiteMetal(double yc,
                                          VectorDouble krigest,
                                          VectorDouble krigstd,
                                          const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteMetalElement(double yc,
                                           double krigest,
                                           double krigstd,
                                           const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteMetalStd(double yc,
                                             VectorDouble krigest,
                                             VectorDouble krigstd,
                                             const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteMetalStdElement(double yc,
                                              double krigest,
                                              double krigstd,
                                              const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteCondExp(VectorDouble krigest,
                                            VectorDouble krigstd,
                                            const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteCondExpElement(double krigest,
                                             double krigstd,
                                             const VectorDouble &phi);
VectorDouble hermiteCondStd(VectorDouble krigest,
                            VectorDouble krigstd,
                            const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteCondStdElement(double krigest,
                                             double krigstd,
                                             const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteEvaluateZ2(VectorDouble yk,
                                               VectorDouble sk,
                                               const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteEvaluateZ2(double yk,
                                         double sk,
                                         const VectorDouble &phi);
