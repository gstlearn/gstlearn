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

#include "Covariances/CovAniso.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Covariances/CovContext.hpp"
#include "Space/SpacePoint.hpp"
#include "gstlearn_export.hpp"


class ACov;
/**
 * \brief
 * This class describes an **elementary covariance**.
 *
 * This covariance is described through the following list of parameters:
 * - the covariance **type**: the list of these types is provided in ECov.hpp
 * - the largest set of parameters for any covariance: **range(s)**, **anisotropy angle(s)**, **third parameter**. Some of these parameters
 * do not make sense, depending on the covariance type: e.g. the range for nugget effect, the third parameter for a spherical
 * structure, ...
 * All these parameters are processed and stored as a **tensor** in order to avoid repetitive calculations.
 * - the **sill**. This comes as a square symmetric matrix whose dimension is equal to the number of variables.
 */
class GSTLEARN_EXPORT CovGneiting: public ACov//, public ICloneable
{
public:
  CovGneiting();  
  CovGneiting(const CovGneiting& r);
  CovGneiting& operator=(const CovGneiting& r);
  virtual ~CovGneiting();


  /// ACov Interface
 
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;
    virtual double eval0(int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode* mode = nullptr) const override;
    virtual void eval0MatInPlace(MatrixSquareGeneral &mat,
                                 const CovCalcMode *mode = nullptr) const override; 
    virtual int getNVariables() const override { return 1; }
    
private:
  CovContext _ctxt;                    /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  CovAniso* _covS;
  CovAniso* _covTemp;

};
