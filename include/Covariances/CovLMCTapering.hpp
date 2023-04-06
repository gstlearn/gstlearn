/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/ETape.hpp"

#include "gstlearn_export.hpp"
#include "Covariances/CovLMC.hpp"

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

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMCTapering)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar,
                      int jvar,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(const ETape& tapetype, double taperange);

  const String& getName() const;
  double getTapeRange() const { return _tapeRange; }
  void setTapeRange(double range) { _tapeRange = range; }

private:
  ETape  _tapeType;
  double _tapeRange;
};
