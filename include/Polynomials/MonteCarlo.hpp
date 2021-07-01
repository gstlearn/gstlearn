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

#define NBSIMU_DEF 1000

double integralGaussHermite(double yc, double r, const VectorDouble& psi);
void normalizeResults(int nbsimu, double& valest, double& valstd);
void normalizeResults(int nbsimu, double& valest);

VectorDouble MCCondExp(VectorDouble krigest,
                       VectorDouble krigstd,
                       const VectorDouble& psi,
                       int nbsimu = NBSIMU_DEF);
double MCCondExpElement(double krigest,
                        double krigstd,
                        const VectorDouble& psi,
                        int nbsimu = NBSIMU_DEF);
VectorDouble MCCondStd(VectorDouble krigest,
                       VectorDouble krigstd,
                       const VectorDouble& psi,
                       int nbsimu = NBSIMU_DEF);
double MCCondStdElement(double krigest,
                        double krigstd,
                        const VectorDouble& psi,
                        int nbsimu = NBSIMU_DEF);
VectorDouble MCIndicator(double yc,
                         VectorDouble krigest,
                         VectorDouble krigstd,
                         int nbsimu = NBSIMU_DEF);
double MCIndicatorElement(double yc,
                          double krigest,
                          double krigstd,
                          int nbsimu = NBSIMU_DEF);
VectorDouble MCIndicatorStd(double yc,
                            VectorDouble krigest,
                            VectorDouble krigstd,
                            int nbsimu = NBSIMU_DEF);
double MCIndicatorStdElement(double yc,
                             double krigest,
                             double krigstd,
                             int nbsimu = NBSIMU_DEF);
VectorDouble MCMetal(double yc,
                     VectorDouble krigest,
                     VectorDouble krigstd,
                     const VectorDouble& psi,
                     int nbsimu = NBSIMU_DEF);
double MCMetalElement(double yc,
                      double krigest,
                      double krigstd,
                      const VectorDouble& psi,
                      int nbsimu = NBSIMU_DEF);
VectorDouble MCMetalStd(double yc,
                        VectorDouble krigest,
                        VectorDouble krigstd,
                        const VectorDouble& psi,
                        int nbsimu = NBSIMU_DEF);
double MCMetalStdElement(double yc,
                         double krigest,
                         double krigstd,
                         const VectorDouble& psi,
                         int nbsimu = NBSIMU_DEF);
