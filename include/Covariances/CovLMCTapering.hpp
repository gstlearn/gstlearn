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

class ASpace;
class SpacePoint;
class CovAniso;
class Model;

typedef struct
{
  std::string name;
  int maxNDim;  /* Maximum dimension for validity */
  double (*tapeFunc)(double);
} Def_Tapering;

/* Prototyping the internal covariance functions */
double _tape_spherical(double);
double _tape_cubic(double);
double _tape_triangle(double);
double _tape_penta(double);
double _tape_storkey(double);
double _tape_wendland1(double);
double _tape_wendland2(double);

GSTLEARN_EXPORT Def_Tapering& D_TAPE(int rank);

class GSTLEARN_EXPORT CovLMCTapering : public CovLMC
{
public:
  CovLMCTapering(int tapetype = 0,
                 double taperange = 0.,
                 const ASpace* space = nullptr);
  CovLMCTapering(const CovLMCTapering &r);
  CovLMCTapering& operator= (const CovLMCTapering &r);
  virtual ~CovLMCTapering();

  virtual IClonable* clone() const override { return new CovLMCTapering(*this); };
  virtual String toString(int /*level*/) const;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(int tapetype, double taperange);

  int getMaxNDim() const { return _maxNDim; }
  const String& getName() const { return _tapeName; }
  double getTapeRange() const { return _tapeRange; }
  void setTapeRange(double range) { _tapeRange = range; }
  void setTapeType(int type) { _tapeType = type; }

private:
  int _getTapeNumber();

private:
  int    _tapeType;
  double _tapeRange;
  int    _maxNDim;
  String _tapeName;
};
