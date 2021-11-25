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

#define NBSIMU_DEF 1000

GSTLEARN_EXPORT double integralGaussHermite(double yc,
                                            double r,
                                            const VectorDouble &psi);
GSTLEARN_EXPORT void normalizeResults(int nbsimu,
                                      double &valest,
                                      double &valstd);
GSTLEARN_EXPORT void normalizeResults(int nbsimu, double &valest);

GSTLEARN_EXPORT VectorDouble MCCondExp(VectorDouble krigest,
                                       VectorDouble krigstd,
                                       const VectorDouble &psi,
                                       int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCCondExpElement(double krigest,
                                        double krigstd,
                                        const VectorDouble &psi,
                                        int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT VectorDouble MCCondStd(VectorDouble krigest,
                                       VectorDouble krigstd,
                                       const VectorDouble &psi,
                                       int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCCondStdElement(double krigest,
                                        double krigstd,
                                        const VectorDouble &psi,
                                        int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT VectorDouble MCIndicator(double yc,
                                         VectorDouble krigest,
                                         VectorDouble krigstd,
                                         int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCIndicatorElement(double yc,
                                          double krigest,
                                          double krigstd,
                                          int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT VectorDouble MCIndicatorStd(double yc,
                                            VectorDouble krigest,
                                            VectorDouble krigstd,
                                            int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCIndicatorStdElement(double yc,
                                             double krigest,
                                             double krigstd,
                                             int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT VectorDouble MCMetal(double yc,
                                     VectorDouble krigest,
                                     VectorDouble krigstd,
                                     const VectorDouble &psi,
                                     int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCMetalElement(double yc,
                                      double krigest,
                                      double krigstd,
                                      const VectorDouble &psi,
                                      int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT VectorDouble MCMetalStd(double yc,
                                        VectorDouble krigest,
                                        VectorDouble krigstd,
                                        const VectorDouble &psi,
                                        int nbsimu = NBSIMU_DEF);
GSTLEARN_EXPORT double MCMetalStdElement(double yc,
                                         double krigest,
                                         double krigstd,
                                         const VectorDouble &psi,
                                         int nbsimu = NBSIMU_DEF);
