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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

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
