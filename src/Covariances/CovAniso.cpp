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
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"

#include "Space/ASpace.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include <math.h>

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int COVWGT[4][5] = { { 2, -2, 0, 0, 0 },
                            { 6, -8, 2, 0, 0 },
                            { 20, -30, 12, -2, 0 },
                            { 70, -112, 56, -16, 2 } };

CovAniso::CovAniso(const ECov& type, const CovContext& ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  _updateFromContext();
}

CovAniso::CovAniso(const String& symbol, const CovContext& ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _cova = CovFactory::createCovFunc(covtype, ctxt);
  _updateFromContext();
}

CovAniso::CovAniso(const ECov& type,
                   double range,
                   double param,
                   double sill,
                   const CovContext& ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This entry is dedicated to the Monovariate case");
    messerr("Additional parameters are ignored and stand constructor is used");
  }
  else
  {
    _updateFromContext();
    _sill.setValue(0, 0, sill);
    setParam(param);
    setRange(range);
  }
}

CovAniso::CovAniso(const CovAniso &r)
    : ACov(r),
      _ctxt(r._ctxt),
      _cova(CovFactory::duplicateCovFunc(*r._cova)),
      _sill(r._sill),
      _aniso(r._aniso)
{
}

CovAniso& CovAniso::operator=(const CovAniso &r)
{
  if (this != &r)
  {
    ACov::operator =(r);
    _ctxt = r._ctxt;
    _cova = CovFactory::duplicateCovFunc(*r._cova);
    _sill = r._sill;
    _aniso = r._aniso;
  }
  return *this;
}

CovAniso::~CovAniso()
{
  delete _cova;
}

void CovAniso::setContext(const CovContext& ctxt)
{
  _ctxt = ctxt;
  _updateFromContext();
}

void CovAniso::setParam(double param)
{
  if (!_cova->hasParam()) return;
  _cova->setParam(param);
}

void CovAniso::setSill(double sill)
{
  /// TODO : Replace throw by error message management
  if (getNVariables() != 1)
  my_throw("Number of provided sill doesn't match number of variables");
  _sill.reset(1, 1, sill);
}

void CovAniso::setSill(const MatrixSquareGeneral& sill)
{
  if (getNVariables() != sill.getNSize())
  my_throw("Number of provided sills doesn't match number of variables");
  _sill = sill;
}

void CovAniso::setSill(const VectorDouble& sill)
{
  int size = static_cast<int> (sill.size());
  int nvar = getNVariables();
  if (size != nvar * nvar)
    my_throw("Number of provided sills doesn't match number of variables");
  _sill.setValues(sill);
}

void CovAniso::setSill(int ivar, int jvar, double sill)
{
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
  if (!_sill.isValid(ivar, jvar))
  my_throw("Invalid indices in setSill");
  _sill.setValue(ivar, jvar, sill);
}

void CovAniso::initSill(double value)
{
  _sill.fill(value);
}

void CovAniso::setRange(double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON20)
  {
    messerr("Range is too small (%lf). It has been replaced by 1.",range);
    range = 1;
  }
  double scadef = _cova->getScadef();
  setScale(range / scadef);
}

void CovAniso::setRanges(const VectorDouble& range)
{
  if (!hasRange()) return;
  if (range.size() != getNDim())
  my_throw("Inconsistency on Space Dimension");
  for (unsigned int i = 0; i < range.size(); i++)
    if (range[i] <= EPSILON20)
    my_throw("A range should not be too small");

  double scadef = _cova->getScadef();
  VectorDouble scale = range;
  ut_vector_divide_inplace(scale, scadef);
  setScales(scale);
}

void CovAniso::setRange(int idim, double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON20)
  my_throw("A range should not be too small");
  double scadef = _cova->getScadef();
  setScale(idim, range / scadef);
}

void CovAniso::setScale(double scale)
{
  if (!hasRange()) return;
  if (scale <= EPSILON20)
  my_throw("A scale should not be too small");
  _aniso.setRadius(scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * scale);
}

void CovAniso::setScales(const VectorDouble& scale)
{
  if (!hasRange()) return;
  for (unsigned int i = 0; i < scale.size(); i++)
    if (scale[i] <= EPSILON20)
    my_throw("A scale should not be too small");
  _aniso.setRadiusVec(scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * ut_vector_max(scale));
}

void CovAniso::setScale(int idim, double scale)
{
  if (scale <= EPSILON20)
  my_throw("A scale should not be too small");
  _aniso.setRadiusDir(idim, scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * ut_vector_max(_aniso.getRadius()));
}

void CovAniso::setAnisoRotation(const Rotation& rot)
{
  if (!hasRange()) return;
  _aniso.setRotation(rot);
}

void CovAniso::setAnisoRotation(const VectorDouble& rot)
{
  if (!hasRange()) return;
  int ndim = getNDim();
  if ((int) rot.size() != ndim * ndim)
    my_throw("Error in dimension of 'rot'");
  Rotation r(ndim);
  r.setMatrixDirectByVector(rot);
  _aniso.setRotation(r);
}

void CovAniso::setAnisoAngles(const VectorDouble& angles)
{
  if (!hasRange()) return;
  _aniso.setRotationAngles(angles);
}

void CovAniso::setAnisoAngle(int idim, double angle)
{
  if (!hasRange()) return;
  _aniso.setRotationAngle(idim, angle);
}

bool CovAniso::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : check something in CovAniso::isConsistent?
  return _cova->isConsistent();
}

double CovAniso::eval0(int ivar, int jvar, const CovCalcMode& mode) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;

  // Calculate unit distance by applying anisotropy
  double cov = _cova->evalCov(0);

  // Converting into a variogram is ignored as the result is obvious
  // and this could be not significant most of the time

  // Scale by the sill
  if (! mode.getUnitary())
  {
    double sill;
    if (mode.getEnvelop() != 0)
    {
      double sign = (mode.getEnvelop() > 0) ? 1 : -1;
      double coef = sqrt(getSill(ivar, ivar) * getSill(jvar, jvar));
      sill = sign * coef;
    }
    else
    {
      sill = getSill(ivar, jvar);
    }
    cov *= sill;
  }

  return (cov);
}

// TODO Replace p1 and p2 by SpaceTarget
double CovAniso::eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  double cov = 0;

  // Calculate unit distance by applying anisotropy
  /// TODO: if composite space : return h1, h2, ... (number of sub-space)
  double h = getSpace()->getDistance(p1, p2, _aniso);

   // Calculation of Standardized Covariance

  int norder = mode.getOrderVario();
  if (norder > 0)
  {

    // Calculate High-order Variogram (only valuable when h != 0)

    for (int iwgt = 1; iwgt < NWGT[norder]; iwgt++)
    {
      double hp = h * (1. + iwgt);
      cov += COVWGT[norder][iwgt] * _cova->evalCov(hp);
    }
    cov /= NORWGT[norder];
  }
  else
  {
    cov = _cova->evalCov(h);

    // Convert into a variogram
    if (mode.getAsVario())
    {
      cov = _cova->evalCov(0) - cov;
    }
  }

  // Scale by the sill
  if (!mode.getUnitary())
  {
    double sill;
    if (mode.getEnvelop() != 0)
    {
      double sign = (mode.getEnvelop() > 0) ? 1 : -1;
      double coef = sqrt(getSill(ivar, ivar) * getSill(jvar, jvar));
      sill = sign * coef;
    }
    else
    {
      sill = getSill(ivar, jvar);
    }
    cov *= sill;
  }

  return (cov);
}

String CovAniso::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  // Covariance Name
  sstr << _cova->toString();

  // Sill - Factor / Slope information
  if (_cova->hasRange() > 0)
  {

    // A sill is defined

    if (getNVariables() > 1)
    {
      sstr << toMatrix("- Sill matrix:",VectorString(),VectorString(),0,
                       getNVariables(),getNVariables(),_sill.getValues());
    }
    else
    {
      sstr << "- Sill         = " << toDouble(getSill(0, 0)) << std::endl;
    }

    // Isotropy vs anisotropy

    if (_aniso.isIsotropic())
    {
      sstr << "- Range        = " << toDouble(getRange(0)) << std::endl;
      if (isAsymptotic())
        sstr << "- Theo. Range  = " << toDouble(getScale(0)) << std::endl;
    }
    else
    {
      sstr << toVector("- Ranges       = ", getRanges());
      if (isAsymptotic())
        sstr << toVector("- Theo. Ranges = ", getScales());
      if (!_aniso.getRotation().isIdentity())
      {
        sstr << toVector("- Angles       = ", getAnisoAngles());
        sstr << toMatrix("- Rotation Matrix",VectorString(),VectorString(),true,
                         getNDim(),getNDim(),getAnisoRotMatVec());
      }
    }
  }
  else if (_cova->hasRange() < 0)
  {
    // The sill is not defined: use slope instead

    if (getNVariables() > 1)
    {
      MatrixSquareSymmetric slopes = _sill;
      double range = getRange(0);
      for (int ivar = 0; ivar < getNVariables(); ivar++)
        for (int jvar = 0; jvar < getNVariables(); jvar++)
          slopes.setValue(ivar, jvar, _sill.getValue(ivar,jvar) / range);
      sstr << toMatrix("- Slope matrix:",VectorString(),VectorString(),0,
                       getNVariables(),getNVariables(),slopes.getValues());
    }
    else
    {
      sstr << "- Slope        = " << toDouble(getSlope(0, 0));
    }

    if (! _aniso.isIsotropic())
    {
      sstr << toVector("- Aniso, Coeff = ", _aniso.getRadius());
      if (!_aniso.getRotation().isIdentity())
      {
        sstr << toVector("- Angles       = ", getAnisoAngles());
        sstr << toMatrix("- Rotation Matrix",VectorString(),VectorString(),true,
                         getNDim(),getNDim(),getAnisoRotMatVec());
      }
    }
  }
  else
  {
    // Only sill is defined

    if (getNVariables() > 1)
     {
       sstr << toMatrix("- Sill matrix:",VectorString(),VectorString(),0,
                        getNVariables(),getNVariables(),_sill.getValues());
     }
     else
     {
       sstr << "- Sill         = " << toDouble(getSill(0, 0)) << std::endl;
     }
  }

  return sstr.str();
}

double CovAniso::getSill(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  if (!_sill.isValid(ivar, jvar))
  my_throw("Wrong variable index while getting the sill");
  return _sill.getValue(ivar, jvar);
}

/**
 * Return the Slope calculated as the sill / range(idim=0)
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @return
 */
double CovAniso::getSlope(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  if (!_sill.isValid(ivar, jvar))
    my_throw("Wrong variable index while getting the sill");
  if (hasRange() == 0) return TEST;
  double range = getRange(0);
  return _sill.getValue(ivar, jvar) / range;
}

VectorDouble CovAniso::getRanges() const
{
  VectorDouble range = getScales();
  double scadef = _cova->getScadef();
  if (!hasRange()) scadef = 0.;
  ut_vector_multiply_inplace(range, scadef);
  return range;
}

void CovAniso::setType(const ECov& type)
{
  delete _cova;
  _cova = CovFactory::createCovFunc(type, _ctxt);
}

/**
 * This function returns the range in the isotropic case
 * In the anisotropic case, it returns the largest range over all directions
 * @return
 */
double CovAniso::getRange() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getRange(0);
  else
    return ut_vector_max(getRanges());
}

double CovAniso::getTheoretical() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getScale(0);
  else
    return ut_vector_max(getScales());
}

const VectorDouble CovAniso::getAnisoCoeffs() const
{
  VectorDouble coef = getRanges();
  double max = ut_vector_max(coef);
  if (max < EPSILON10)
  my_throw("Range is null");
  ut_vector_divide_inplace(coef, max);
  return coef;
}

/**
 * For compatibility, this function returns 0 if the Covariance has no Third Parameter
 *
 * @return Third parameter
 */
double CovAniso::getParam() const
{
  if (!hasParam())
    return 0.;
  else
    return _cova->getParam();
}

void CovAniso::_updateFromContext()
{
  int ndim = getNDim();
  int nvar = getNVariables();
  _sill.reset(nvar, nvar, 1.);
  _aniso.init(ndim);
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  int ndim = getNDim();
  VectorDouble dd(ndim);
  double delta = hmax / ndisc;
  double total = 0.;
  switch (ndim)
  {
    case 1:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
      {
        dd[0] = delta * j1;
        total += delta * eval(0,0,dd,SpacePoint());
      }
      break;

    case 2:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
        {
          dd[0] = delta * j1;
          dd[1] = delta * j2;
          total += delta * delta * eval(0,0,dd,SpacePoint());
        }
      break;

    case 3:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
          for (int j3 = -ndisc; j3 <= ndisc; j3++)
          {
            dd[0] = delta * j1;
            dd[1] = delta * j2;
            dd[2] = delta * j3;
            total += delta * delta * delta * eval(0,0,dd,SpacePoint());
          }
      break;

    default:
      my_throw("Integral Range has been programmed for Space Dimensio  between 1 and 3");
  }
  return total;
}

bool CovAniso::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= getNVariables())
  {
    mesArg("Rank of the Variable",1,getNVariables());
    return false;
  }
  return true;
}

int CovAniso::getGradParamNumber() const
{
  int ndim = getNDim();
  int number = 0;

  // Anisotropy ranges
  if (hasRange())
  {
    // Anisotropy ranges
    number += ndim;

    // Anisotropy Rotation angles
    if (ndim == 2)
      number ++;
    else
      number += ndim;
  }
  return number;
}

double CovAniso::scale2range(const ECov &type, double scale, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double factor = cova->getScadef();
  return scale * factor;
}

double CovAniso::range2scale(const ECov &type, double range, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double factor = cova->getScadef();
  return range / factor;
}

CovAniso* CovAniso::createIsotropic(const CovContext& ctxt,
                                    const ECov& type,
                                    double range,
                                    double sill,
                                    double param)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  return new CovAniso(type, range, param, sill, ctxt);
}

CovAniso* CovAniso::createAnisotropic(const CovContext& ctxt,
                                      const ECov& type,
                                      const VectorDouble& ranges,
                                      double sill,
                                      double param)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CovAniso* cov = new CovAniso(type, ctxt);
  cov->setRanges(ranges);
  cov->setSill(sill);
  cov->setParam(param);
  return cov;
}

CovAniso* CovAniso::createIsotropicMulti(const CovContext& ctxt,
                                         const ECov& type,
                                         double range,
                                         const MatrixSquareGeneral& sills,
                                         double param)
{
  CovAniso* cov = new CovAniso(type, ctxt);
  int nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr("Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
            nvar,ctxt.getNVar());
    return nullptr;
  }
  cov->setRange(range);
  cov->setSill(sills);
  cov->setParam(param);
  return cov;
}

CovAniso* CovAniso::createAnisotropicMulti(const CovContext& ctxt,
                                           const ECov& type,
                                           const VectorDouble& ranges,
                                           const MatrixSquareGeneral& sills,
                                           double param)
{

  int nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr(
        "Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
        nvar, ctxt.getNVar());
    return nullptr;
  }
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CovAniso* cov = new CovAniso(type, ctxt);
  cov->setRanges(ranges);
  cov->setSill(sills);
  cov->setParam(param);
  return cov;
}

void CovAniso::copyCovContext(const CovContext& ctxt)
{
  _ctxt.copyCovContext(ctxt);
  if (_cova != nullptr) _cova->copyCovContext(ctxt);
}

