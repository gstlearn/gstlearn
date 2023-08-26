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
#include "Enum/ETape.hpp"

#include "Covariances/CovLMCTapering.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"

#include <math.h>

CovLMCTapering::CovLMCTapering(const ETape& tapetype,
                               double taperange,
                               const ASpace* space)
    : CovLMC(space),
      _tapeType(),
      _tapeRange(0)
{
  init(tapetype, taperange);
}

CovLMCTapering::CovLMCTapering(const CovLMCTapering &r)
    : CovLMC(r),
      _tapeType(r._tapeType),
      _tapeRange(r._tapeRange)
{
}

CovLMCTapering& CovLMCTapering::operator=(const CovLMCTapering &r)
{
  if (this != &r)
  {
    CovLMC::operator=(r);
    _tapeType = r._tapeType;
    _tapeRange = r._tapeRange;
  }
  return *this;
}

CovLMCTapering::~CovLMCTapering()
{
}

int CovLMCTapering::init(const ETape& tapetype, double taperange)
{

  /* Preliminary check */

  if (taperange <= 0)
  {
    messerr("The argument 'tape_range' must be strictly positive");
    return 1;
  }

  /* Load the tapering parameters */

  _tapeType    = tapetype;
  _tapeRange   = taperange;

  return 0;
}

Def_Tapering& D_TAPE(int rank)
{
  static Def_Tapering DEF_TAPES[] =
  {
   {"Spherical"    , 3, _tape_spherical },
   {"Cubic"        , 3, _tape_cubic     },
   {"Triangle"     , 1, _tape_triangle  },
   {"Pentamodel"   , 3, _tape_penta     },
   {"Storkey"      , 1, _tape_storkey   },
   {"Wendland1"    , 3, _tape_wendland1 },
   {"Wendland2"    , 3, _tape_wendland2 }
  };
  return DEF_TAPES[rank];
}

double _tape_spherical(double h)
{
  double cov;

  cov = 0.;
  if (h < 1) cov = 1 - 0.5 * h * (3 - h * h);

  return (cov);
}
double _tape_cubic(double h)
{
  double h2, cov;

  cov = 0.;
  h2 = h * h;
  if (h < 1) cov = 1. - h2 * (7. + h * (-8.75 + h2 * (3.5 - 0.75 * h2)));
  cov = MAX(0., cov);

  return (cov);
}

double _tape_triangle(double h)
{
  double cov;

  cov = MAX(0, 1. - h);
  return (cov);
}

double _tape_penta(double h)
{
  double h2, cov;

  cov = 0.;
  h2 = h * h;
  if (h < 1)
    cov = 1.
        - h2 * (22. / 3.
            - h2 * (33.
                - h * (77. / 2.
                    - h2 * (33. / 2. - h2 * (11. / 2. - 5. / 6. * h2)))));

  return (cov);
}

double _tape_storkey(double h)
{
  double cov, pi2;

  cov = 0.;
  pi2 = 2. * GV_PI;
  if (h < 1)
    cov = (2. * (1. - h) * (1. + cos(pi2 * h) / 2.) + 3 / pi2 * sin(pi2 * h))
        / 3.;

  return (cov);
}

double _tape_wendland1(double h)
{
  double h2, cov;

  cov = 0.;
  h2 = h * h;
  if (h < 1) cov = 1 - h2 * (10 - h * (20 - h * (15 - h * 4)));
  return (cov);
}

double _tape_wendland2(double h)
{
  double h2, cov;

  cov = 0.;
  h2 = h * h;
  if (h < 1)
    cov = 1
        - h2 * ((28. / 3.)
            - h2 * (70
                - h * ((448. / 3.) - h * (140 - h * (64 - h * (35. / 3.))))));
  return (cov);
}

String CovLMCTapering::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << ACovAnisoList::toString(strfmt);

  sstr << "Tapering Function     = " << getName() << std::endl;
  sstr << "Tapering Scale        = " << _tapeRange << std::endl;
  return sstr.str();
}

double CovLMCTapering::eval0(int ivar,
                             int jvar,
                             const CovCalcMode* mode) const
{
  double cov0 = CovLMC::eval0(ivar, jvar, mode);
  return cov0;
}

double CovLMCTapering::eval(const SpacePoint& p1,
                            const SpacePoint& p2,
                            int ivar,
                            int jvar,
                            const CovCalcMode* mode) const
{
  // The calculation flag 'as.Vario' must be treated here rather than relying on calculation
  // performed in generic 'eval' method
  double cov = 0.;
  double cov0 = 0.;
  bool asVario = false;
  if (mode == nullptr)
  {
    cov = CovLMC::eval(p1, p2, ivar, jvar);
  }
  else
  {
    CovCalcMode modeloc(*mode);
    asVario = mode->getAsVario();
    modeloc.setAsVario(false);
    cov = CovLMC::eval(p1, p2, ivar, jvar, &modeloc);
    cov0 = CovLMC::eval(p1, p1, ivar, jvar, &modeloc); // or eval0 if stationary
  }

  double h = getSpace()->getDistance(p1, p2) / _tapeRange;
  cov *= D_TAPE(_tapeType.getValue()).tapeFunc(h);

  if (asVario) cov = cov0 - cov;
  return cov;
}

void CovLMCTapering::eval0MatInPlace(MatrixSquareGeneral &mat,
                                     const CovCalcMode *mode) const
{
  // We do not want to call the optimization of ACovAnisoList
  ACov::eval0MatInPlace(mat, mode);
}

void CovLMCTapering::evalMatInPlace(const SpacePoint &p1,
                                    const SpacePoint &p2,
                                    MatrixSquareGeneral &mat,
                                    const CovCalcMode *mode) const
{
  // We do not want to call the optimization of ACovAnisoList
  ACov::evalMatInPlace(p1, p2, mat, mode);
}

const String& CovLMCTapering::getName() const
{
  return _tapeType.getDescr();
}
