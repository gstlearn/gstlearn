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
#include "geoslib_old_f.h"

#include "Enum/ECov.hpp"

#include "Basic/Utilities.hpp"
#include "Model/Tapering.hpp"

#include <math.h>

Tapering::Tapering()
  : AStringable(),
    _type(0)
  , _maxNDim(0)
  , _range(0)
{
}

Tapering::Tapering(const Tapering &m)
    : AStringable(m),
      _type(m._type),
      _maxNDim(m._maxNDim),
      _range(m._range)
{
}

Tapering& Tapering::operator=(const Tapering &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _type = m._type;
    _maxNDim = m._maxNDim;
    _range = m._range;
  }
  return (*this);
}

Tapering::~Tapering()
{
}

int Tapering::init(int tape_type,double tape_range)
{

  /* Preliminary check */

  if (tape_type < 1 || tape_type > getTapeNumber())
  {
    mesArg("Tapering Index",tape_type,getTapeNumber());
    return 1;
  }
  if (tape_range <= 0)
  {
    messerr("The argument 'tape_range' must be strictly positive");
    return 1;
  }

  /* Load the tapering parameters */

  _type    = tape_type - 1;
  _range   = tape_range;
  _name    = D_TAPE(_type).tapeName;
  _maxNDim = D_TAPE(_type).maxNDim;

  return 0;
}

int Tapering::getTapeNumber()
{
  int N_DEF_TAPERING = 7;
  return N_DEF_TAPERING;
}

Def_Tapering& D_TAPE(int rank)
{
  static Def_Tapering DEF_TAPES[] =
  {
   {"Spherical"    , ECov::E_SPHERICAL, 3, _tape_spherical },
   {"Cubic"        , ECov::E_CUBIC,     3, _tape_cubic     },
   {"Triangle"     , ECov::E_TRIANGLE,  1, _tape_triangle  },
   {"Pentamodel"   , ECov::E_PENTA,     3, _tape_penta     },
   {"Storkey"      , ECov::E_STORKEY,   1, _tape_storkey   },
   {"Wendland1"    , ECov::E_WENDLAND1, 3, _tape_wendland1 },
   {"Wendland2"    , ECov::E_WENDLAND2, 3, _tape_wendland2 }
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

String Tapering::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Tapering Function     = " << _name << std::endl;
  sstr << "Tapering Scale        = " << _range << std::endl;
  return sstr.str();
}
