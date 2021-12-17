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
#include "Basic/Utilities.hpp"
#include "Model/Convolution.hpp"
#include "geoslib_enum.h"
#include "geoslib_f.h"
#include <math.h>

Convolution::Convolution()
  : AStringable()
  , _type(0)
  , _dir(0)
  , _discNumber(0)
  , _count(0)
  , _range(0)
  , _delX()
  , _delY()
  , _delZ()
  , _weight()
{
  _convFunc = nullptr;
}

Convolution::Convolution(const Convolution &m)
    : AStringable(m),
      _type(m._type),
      _dir(m._dir),
      _discNumber(m._discNumber),
      _count(m._count),
      _range(m._range),
      _delX(m._delX),
      _delY(m._delY),
      _delZ(m._delZ),
      _weight(m._weight)
{
  _convFunc = m._convFunc;
}

Convolution& Convolution::operator=(const Convolution &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _type = m._type;
    _dir = m._dir,
    _discNumber = m._discNumber,
    _count = m._count;
    _range = m._range;
    _delX = m._delX;
    _delY = m._delY;
    _delZ = m._delZ;
    _weight = m._weight;
    _convFunc = m._convFunc;
  }
  return (*this);
}

Convolution::~Convolution()
{
}

int Convolution::init(int conv_type,
                      int conv_idir,
                      int conv_ndisc,
                      double conv_range)
{
  int navail[3];
  double total[2];

  // Preliminary checks

  if (conv_type < 1 || conv_type > getConvNumber())
  {
    mesArg("Convolution Type", conv_type, getConvNumber());
    return 1;
  }
  if (conv_idir < 1 || conv_idir > 5)
  {
    mesArg("Convolution Direction Index", conv_idir, 5);
    messerr("The argument 'conv_idir' should lie between 1 and 5");
    return 1;
  }
  if (conv_ndisc < 1)
  {
    messerr("The number of discretization points must be larger than 1");
    return 1;
  }
  if (conv_range <= 0.)
  {
    messerr("The convolution range should be strictly positive");
    return 1;
  }

  /* Load the Convolution parameters */

  _type = conv_type;
  _dir = conv_idir;
  _discNumber = conv_ndisc;
  _range = conv_range;
  _convFunc = D_CONV(conv_type).convFunc;
  _convName = D_CONV(conv_type).convName;
  double diameter = _range * D_CONV(conv_type).convScale;

  /* Calculate the discretization points */

  for (int i = 0; i < 3; i++)
    navail[i] = 0;
  switch (conv_idir)
  {
    case CONV_DIRX:
      navail[0] = 1;
      break;

    case CONV_DIRY:
      navail[1] = 1;
      break;

    case CONV_DIRZ:
      navail[2] = 1;
      break;

    case CONV_DIRXY:
      navail[0] = navail[1] = 1;
      break;

    case CONV_DIRXYZ:
      navail[0] = navail[1] = navail[2] = 1;
      break;

    default:
      messerr("This convolution direction (%d) does not exist", _dir);
      return 1;
  }

  /* Calculate the discretization points */

  int number = 1;
  for (int i = 0; i < 3; i++)
    if (navail[i] > 0) number *= 2 * _discNumber + 1;
  _delX.resize(number * 2);
  _delX.resize(number * 2);
  _delX.resize(number * 2);
  _weight.resize(number * 2);
  _count = 2 * number;

  total[0] = total[1] = 0.;
  int ecr = 0;
  for (int iconv = 0; iconv < 2; iconv++)
    for (int ix = -conv_ndisc * navail[0]; ix <= conv_ndisc * navail[0]; ix++)
      for (int iy = -conv_ndisc * navail[0]; iy <= conv_ndisc * navail[0]; iy++)
        for (int iz = -conv_ndisc * navail[0]; iz <= conv_ndisc * navail[0];
            iz++)
        {
          _delX[ecr] = diameter * iconv * ix / (2 * conv_ndisc + 1);
          double vx = (navail[0] > 0) ? _convFunc(iconv, _delX[ecr]) :
                                        1.;
          _delY[ecr] = diameter * iconv * iy / (2 * conv_ndisc + 1);
          double vy = (navail[1] > 0) ? _convFunc(iconv, _delY[ecr]) :
                                        1.;
          _delZ[ecr] = diameter * iconv * iz / (2 * conv_ndisc + 1);
          double vz = (navail[2] > 0) ? _convFunc(iconv, _delZ[ecr]) :
                                        1.;
          _weight[ecr] = vx * vy * vz;
          total[iconv] += _weight[ecr];
          ecr++;
        }

  /* Normalize the weights */

  ecr = 0;
  for (int iconv = 0; iconv < 2; iconv++)
    for (int i = 0; i < number; i++, ecr++)
      _weight[ecr] /= total[iconv];

  return 0;
}

int Convolution::getConvNumber()
{
  int N_DEF_CONV = 4;
  return N_DEF_CONV;
}

/*****************************************************************************/
/*!
**  Generic convolution function
**
** \return  Weight of the convolution function
**
** \param[in]  number number of times the convolution function is applied
** \param[in]  v      normalized distance for convolution calculation
**
*****************************************************************************/
double _conv_uniform(int number, double v)
{
  double dp;

  dp = (number == 1) ? 1. : 2. - v;
  return(dp);
}

double _conv_exponential(int number, double v)
{
  double dp;

  dp = (number == 1) ? exp(-v): (1. + v) * exp(-v);
  return(dp);
}

double _conv_gaussian(int number, double v)
{
  double dp;

  dp = (number == 1) ? exp(-v*v/2.): exp(-v*v/4.);
  return(dp);
}

double _conv_sincard(int number, double v)
{
  double dp,v2,dv;

  v = (number == 1) ? GV_PI * v : 2. * GV_PI * v;
  v2 = v * v;
  if (v < 1.0e-10)
  {
    dp = (number == 1) ? 1. - v2/3: (1. - v2/20.) / 6.;
  }
  else
  {
    dv = sin(v) / v;
    dp = (number == 1) ? dv * dv: (1.-dv)/v2;
  }
  return(dp);
}

String Convolution::toString(int /*level*/) const
{
  std::stringstream sstr;

  sstr << "Convolution type      = " << _convName  << std::endl;
  sstr << "Convolution direction = " << _dir   << std::endl;
  sstr << "Nb. discretization    = " << _discNumber << std::endl;
  sstr << "Convolution Scale     = " << _range << std::endl;

  return sstr.str();
}

Def_Conv& D_CONV(int rank)
{
  static Def_Conv DEF_CONVS[] =
  {
   {"Uniform"     , 1.,       _conv_uniform     },
   {"Exponential" , 2.995732, _conv_exponential },
   {"Gaussian"    , 1.730818, _conv_gaussian    },
   {"Sincard"     , 20.371,   _conv_sincard     }
  };

  return DEF_CONVS[rank];
}
