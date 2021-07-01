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

#include "Covariances/ACovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovAniso;

class CovLMC : public ACovAnisoList
{
public:
  CovLMC(const ASpace* space = nullptr);
  CovLMC(const CovLMC &r);
  CovLMC& operator= (const CovLMC &r);
  virtual ~CovLMC();

  virtual IClonable* clone() const override;
};
