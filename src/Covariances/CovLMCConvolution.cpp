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
#include "Covariances/CovLMCConvolution.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "geoslib_f.h"

#include <math.h>

CovLMCConvolution::CovLMCConvolution(int conv_type,
                                     int conv_dir,
                                     double conv_range,
                                     int conv_ndisc,
                                     const ASpace* space)
    : CovLMC(space),
      _convType(conv_type),
      _convDir(conv_dir),
      _convRange(conv_range),
      _convDiscNumber(conv_ndisc),
      _convNumber(0),
      _convIncr(),
      _convWeight()
{
  init(conv_type, conv_dir, conv_range, conv_ndisc);
}

CovLMCConvolution::CovLMCConvolution(const CovLMCConvolution &r)
    : CovLMC(r),
      _convType(r._convType),
      _convDir(r._convDir),
      _convRange(r._convRange),
      _convDiscNumber(r._convDiscNumber),
      _convNumber(r._convNumber),
      _convIncr(r._convIncr),
      _convWeight(r._convWeight)
{
}

CovLMCConvolution& CovLMCConvolution::operator=(const CovLMCConvolution &r)
{
  if (this != &r)
  {
    CovLMC::operator=(r);
    _convType = r._convType;
    _convDir = r._convDir;
    _convRange = r._convRange;
    _convDiscNumber = r._convDiscNumber;
    _convNumber = r._convNumber;
    _convIncr = r._convIncr;
    _convWeight = r._convWeight;
  }
{
  }
  return *this;
}

CovLMCConvolution::~CovLMCConvolution()
{
}

int CovLMCConvolution::init(int conv_type,
                            int conv_idir,
                            double conv_range,
                            int conv_ndisc)
{
  if (conv_type < 0 || conv_type >= _getConvNumber())
  {
    mesArg("CovLMCConvolution Type", conv_type, _getConvNumber());
    return 1;
  }
  if (conv_idir < 0 || conv_idir >= 4)
  {
    mesArg("CovLMCConvolution Direction Index", conv_idir, 4);
    messerr("The argument 'conv_idir' should lie between 0 and 4");
    return 1;
  }
  if (conv_ndisc < 1)
  {
    messerr("The number of discretization points must be larger than 1");
    return 1;
  }
  if (conv_range <= 0.)
  {
    messerr("The CovLMCConvolution range should be strictly positive");
    return 1;
  }

  /* Load the CovLMCConvolution parameters */

  int ndim = getNDim();
  _convType = conv_type;
  _convDir  = conv_idir;
  _convDiscNumber = conv_ndisc;
  _convRange = conv_range;
  double diameter = _convRange * D_CONV(conv_type).convScale;

  /* Calculate the discretization points */

  int navail[3];
  for (int i = 0; i < 3; i++) navail[i] = 0;
  switch (conv_idir)
  {
    case 0:
      navail[0] = (ndim >= 1) ? _convDiscNumber : 0;
      break;

    case 1:
      navail[1] = (ndim >= 2) ? _convDiscNumber : 0;
      break;

    case 2:
      navail[2] = (ndim >= 3) ? _convDiscNumber : 0;
      break;

    case 3:
      navail[0] = (ndim >= 1) ? _convDiscNumber : 0;
      navail[1] = (ndim >= 2) ? _convDiscNumber : 0;
      break;

    case 4:
      navail[0] = (ndim >= 1) ? _convDiscNumber : 0;
      navail[1] = (ndim >= 2) ? _convDiscNumber : 0;
      navail[2] = (ndim >= 3) ? _convDiscNumber : 0;
      break;

    default:
      messerr("This CovLMCConvolution direction (%d) does not exist", _convDir);
      return 1;
  }

  /* Calculate the total number of discretization points */

  _convNumber = 1;
  for (int i = 0; i < 3; i++)
    _convNumber *= 2 * navail[i] + 1;

  _convIncr = MatrixRectangular(ndim,_convNumber);
  _convWeight.resize(_convNumber);

  double delta, weight;
  int ecr = 0;
  double total = 0.;
  for (int ix = -navail[0]; ix <= navail[0]; ix++)
    for (int iy = -navail[1]; iy <= navail[1]; iy++)
      for (int iz = -navail[2]; iz <= navail[2]; iz++)
      {
        double local = 1.;
        if (ndim >= 1)
        {
          delta  = diameter * ix / (2 * conv_ndisc + 1);
          weight = D_CONV(conv_type).convFunc(delta);
          _convIncr.setValue(0, ecr, delta);
          local *= weight;
        }
        if (ndim >= 2)
        {
          delta = diameter * iy / (2 * conv_ndisc + 1);
          weight = D_CONV(conv_type).convFunc(delta);
          _convIncr.setValue(1, ecr, delta);
          local *= weight;
        }
        if (ndim >= 3)
        {
          delta = diameter * iz / (2 * conv_ndisc + 1);
          weight = D_CONV(conv_type).convFunc(delta);
          _convIncr.setValue(2, ecr, delta);
          local *= weight;
        }
        _convWeight[ecr] = local;
        total += local;
        ecr++;
      }

  /* Normalize the weights */

  for (int i = 0; i < _convNumber; i++) _convWeight[i] /= total;

  return 0;
}

Def_Convolution& D_CONV(int rank)
{
  static Def_Convolution DEF_CONVS[] =
  {
   {"Uniform"     , 1.,       _conv_uniform     },
   {"Exponential" , 2.995732, _conv_exponential },
   {"Gaussian"    , 1.730818, _conv_gaussian    },
   {"Sincard"     , 20.371,   _conv_sincard     }
  };
  return DEF_CONVS[rank];
}

double _conv_uniform(double /* v */)
{
  double dp = 1.;
  return(dp);
}

double _conv_exponential(double v)
{
  double dp = exp(-v);
  return(dp);
}

double _conv_gaussian(double v)
{
  double dp = exp(-v*v/2.);
  return(dp);
}

double _conv_sincard(double v)
{
  double dp,v2,dv;

  v = GV_PI * v;
  v2 = v * v;
  if (v < 1.0e-10)
  {
    dp = 1. - v2/3;
  }
  else
  {
    dv = sin(v) / v;
    dp = dv * dv;
  }
  return(dp);
}

String CovLMCConvolution::toString(int level) const
{
  std::stringstream sstr;

  sstr << ACovAnisoList::toString(level);

  sstr << "Convolution type      = " << D_CONV(_convType).convName  << std::endl;
  sstr << "Convolution direction = " << _convDir   << std::endl;
  sstr << "Nb. discretization    = " << _convDiscNumber << std::endl;
  sstr << "Convolution Scale     = " << _convRange << std::endl;

  return sstr.str();
}

double CovLMCConvolution::eval0(int ivar,
                                int jvar,
                                const CovCalcMode& mode) const
{
  double cov0 = 0.;
  SpacePoint p11;
  SpacePoint p22;
  for (int i1 = 0; i1 < _convNumber; i1++)
  {
    double w1 = _convWeight[i1];
    p11.move(_convIncr.getColumn(i1));
    for (int i2 = 0; i2 < _convNumber; i2++)
    {
      double w2 = _convWeight[i2];
      p22.move(_convIncr.getColumn(i2));
      cov0 += CovLMC::eval(ivar, jvar, p11, p22, mode) * w1 * w2;
    }
  }
  return cov0;
}

double CovLMCConvolution::eval(int ivar,
                               int jvar,
                               const SpacePoint& p1,
                               const SpacePoint& p2,
                               const CovCalcMode& mode) const
{
  // The calculation flag 'as.Vario' must be treated here rather than relying on calculation
  // performed internally in 'eval' function
  CovCalcMode modeloc(mode);
  bool asVario = mode.getAsVario();
  modeloc.setAsVario(false);

  double cov = 0.;
  SpacePoint p11(p1);
  SpacePoint p22(p2);
  for (int i1 = 0; i1 < _convNumber; i1++)
  {
    double w1 = _convWeight[i1];
    p11.move(_convIncr.getColumn(i1));
    for (int i2 = 0; i2 < _convNumber; i2++)
    {
      double w2 = _convWeight[i2];
      p22.move(_convIncr.getColumn(i2));
      cov += CovLMC::eval(ivar, jvar, p11, p22, modeloc) * w1 * w2;
    }
  }

  if (asVario)
  {
    double cov0 = eval0(ivar,jvar,modeloc);
    cov = cov0 - cov;
  }
  return cov;
}

int CovLMCConvolution::_getConvNumber()
{
  int N_DEF_CONV = 4;
  return N_DEF_CONV;
}
