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
#include "Enum/EConvDir.hpp"
#include "Enum/EConvType.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Matrix/MatrixRectangular.hpp"

#include <math.h>

CovLMCConvolution::CovLMCConvolution(const EConvType& conv_type,
                                     const EConvDir&  conv_dir,
                                     double conv_range,
                                     int conv_ndisc,
                                     const ASpace* space)
    : CovAnisoList(space),
      _convType(conv_type),
      _convDir(conv_dir),
      _convDiscNumber(conv_ndisc),
      _convRange(conv_range),
      _convNumber(0),
      _convIncr(),
      _convWeight()
{
  init(conv_type, conv_dir, conv_range, conv_ndisc);
}

CovLMCConvolution::CovLMCConvolution(const CovLMCConvolution &r)
    : CovAnisoList(r),
      _convType(r._convType),
      _convDir(r._convDir),
      _convDiscNumber(r._convDiscNumber),
      _convRange(r._convRange),
      _convNumber(r._convNumber),
      _convIncr(r._convIncr),
      _convWeight(r._convWeight)
{
}

CovLMCConvolution& CovLMCConvolution::operator=(const CovLMCConvolution &r)
{
  if (this != &r)
  {
    CovAnisoList::operator=(r);
    _convType = r._convType;
    _convDir = r._convDir;
    _convDiscNumber = r._convDiscNumber;
    _convRange = r._convRange;
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

void CovLMCConvolution::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  ACov::_loadAndAddEvalCovMatBiPointInPlace(mat, p1, p2, mode);
}
void CovLMCConvolution::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                                     const SpacePoint &pwork1,
                                                     const SpacePoint &pwork2,
                                                     const CovCalcMode *mode) const
{
  ACov::_addEvalCovMatBiPointInPlace(mat, pwork1, pwork2, mode);
}
int CovLMCConvolution::init(const EConvType& conv_type,
                            const EConvDir&  conv_idir,
                            double conv_range,
                            int conv_ndisc)
{
  for (auto &e: _covAnisos)
  {
    e->setOptimEnabled(false);
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
  double diameter = _convRange * D_CONV(_convType.getValue()).convScale;

  /* Calculate the discretization points */

  int navail[3];
  for (int i = 0; i < 3; i++) navail[i] = 0;
  switch (_convDir.getValue())
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
      messerr("This CovLMCConvolution direction (%d) does not exist",
              _convDir.getValue());
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
          weight = D_CONV(_convType.getValue()).convFunc(delta);
          _convIncr.setValue(0, ecr, delta);
          local *= weight;
        }
        if (ndim >= 2)
        {
          delta = diameter * iy / (2 * conv_ndisc + 1);
          weight = D_CONV(_convType.getValue()).convFunc(delta);
          _convIncr.setValue(1, ecr, delta);
          local *= weight;
        }
        if (ndim >= 3)
        {
          delta = diameter * iz / (2 * conv_ndisc + 1);
          weight = D_CONV(_convType.getValue()).convFunc(delta);
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

String CovLMCConvolution::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << CovAnisoList::toString(strfmt);

  sstr << "Convolution type      = " << _convType.getDescr() << std::endl;
  sstr << "Convolution direction = " << _convDir.getDescr()   << std::endl;
  sstr << "Nb. discretization    = " << _convDiscNumber << std::endl;
  sstr << "Convolution Scale     = " << _convRange << std::endl;

  return sstr.str();
}

double CovLMCConvolution::eval0(int ivar,
                                int jvar,
                                const CovCalcMode* mode) const
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
      cov0 += CovAnisoList::eval(p11, p22, ivar, jvar, mode) * w1 * w2;
    }
  }
  return cov0;
}

double CovLMCConvolution::eval(const SpacePoint& p1,
                               const SpacePoint& p2,
                               int ivar,
                               int jvar,
                               const CovCalcMode* mode) const
{
  SpacePoint p11;
  SpacePoint p22;

  // The calculation flag 'as.Vario' must be treated here rather than relying on calculation
  // performed internally in 'eval' function

  bool asVario = false;
  CovCalcMode modeloc;
  if (mode != nullptr)
  {
    modeloc = *mode;
    asVario = mode->getAsVario();
    modeloc.setAsVario(false);
  }

  double cov = 0.;
  p11 = p1;
  p22 = p2;
  for (int i1 = 0; i1 < _convNumber; i1++)
  {
    double w1 = _convWeight[i1];
    p11.move(_convIncr.getColumn(i1));
    for (int i2 = 0; i2 < _convNumber; i2++)
    {
      double w2 = _convWeight[i2];
      p22.move(_convIncr.getColumn(i2));
      double covloc = 0.;
      if (mode == nullptr)
        covloc = CovAnisoList::eval(p11, p22, ivar, jvar);
      else
        covloc = CovAnisoList::eval(p11, p22, ivar, jvar, &modeloc);
      cov += covloc * w1 * w2;
    }
  }

  if (asVario)
  {
    double cov0 = 0.;
    p11 = p1;
    p22 = p1;
    for (int i1 = 0; i1 < _convNumber; i1++)
    {
      double w1 = _convWeight[i1];
      p11.move(_convIncr.getColumn(i1));
      for (int i2 = 0; i2 < _convNumber; i2++)
      {
        double w2 = _convWeight[i2];
        p22.move(_convIncr.getColumn(i2));
        double covloc = 0.;
        if (mode == nullptr)
          covloc = CovAnisoList::eval(p11, p22, ivar, jvar);
        else
        {
          covloc = CovAnisoList::eval(p11, p22, ivar, jvar, &modeloc);
        }
        cov0 += covloc * w1 * w2;
      }
    }
    cov = cov0 -cov;
  }
  return cov;
}
