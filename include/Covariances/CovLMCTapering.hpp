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
#include "Covariances/CovLMC.hpp"
#include "Covariances/ETape.hpp"

class ASpace;
class SpacePoint;
class CovAniso;
class Model;

typedef struct
{
  String name; // Useless but allows checking the correct match
  int maxNDim;  /* Maximum dimension for validity */
  double (*tapeFunc)(double);
} Def_Tapering;

/* Prototyping the internal covariance functions */
GSTLEARN_EXPORT double _tape_spherical(double);
GSTLEARN_EXPORT double _tape_cubic(double);
GSTLEARN_EXPORT double _tape_triangle(double);
GSTLEARN_EXPORT double _tape_penta(double);
GSTLEARN_EXPORT double _tape_storkey(double);
GSTLEARN_EXPORT double _tape_wendland1(double);
GSTLEARN_EXPORT double _tape_wendland2(double);

GSTLEARN_EXPORT Def_Tapering& D_TAPE(int rank);

class GSTLEARN_EXPORT CovLMCTapering : public CovLMC
{
public:
  CovLMCTapering(const ETape& tapetype,
                 double taperange,
                 const ASpace* space = nullptr);
  CovLMCTapering(const CovLMCTapering &r);
  CovLMCTapering& operator= (const CovLMCTapering &r);
  virtual ~CovLMCTapering();

  virtual IClonable* clone() const override { return new CovLMCTapering(*this); };
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(const ETape& tapetype, double taperange);

  const String& getName() const;
  double getTapeRange() const { return _tapeRange; }
  void setTapeRange(double range) { _tapeRange = range; }

private:
  ETape  _tapeType;
  double _tapeRange;
};
