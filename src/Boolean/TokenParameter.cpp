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
#include "Boolean/TokenParameter.hpp"
#include "Boolean/ETLaw.hpp"
#include "Basic/Law.hpp"

TokenParameter::TokenParameter()
    : AStringable(),
      _law(ETLaw::CONSTANT),
      _valarg()
{
}

TokenParameter::TokenParameter(const TokenParameter &r)
    : AStringable(r),
      _law(r._law),
      _valarg(r._valarg)
{
}

TokenParameter& TokenParameter::operator=(const TokenParameter &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _law = r._law;
    _valarg = r._valarg;
  }
  return *this;
}

TokenParameter::~TokenParameter()
{
}

double TokenParameter::getValue() const
{
  if (_law == ETLaw::CONSTANT)
    return _valarg[0];

  if (_law == ETLaw::UNIFORM)
    return law_uniform(_valarg[0], _valarg[1]);

  if (_law == ETLaw::GAUSSIAN)
    return _valarg[0] + _valarg[1] * law_gaussian();

  if (_law == ETLaw::EXPONENTIAL)
    return _valarg[0] + law_exponential(_valarg[1]);

  if (_law == ETLaw::GAMMA)
    return _valarg[0] + law_gamma(_valarg[1]);

  if (_law == ETLaw::STABLE)
    return law_stable(_valarg[0], _valarg[1], _valarg[2], _valarg[3]);

  if (_law == ETLaw::BETA1)
    return law_beta1(_valarg[0], _valarg[1]);

  if (_law == ETLaw::BETA2)
    return law_beta2(_valarg[0], _valarg[1]);

  return TEST;
}

String TokenParameter::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  switch (_law.toEnum())
  {
    case ETLaw::E_CONSTANT:
      sstr << " Constant=" << _valarg[0] << std::endl;
      break;

    case ETLaw::E_UNIFORM:
      sstr << " Uniform within [" << _valarg[0] << ";" << _valarg[1];
      break;

    case ETLaw::E_GAUSSIAN:
      sstr << " Gaussian - Mean=" << _valarg[0] << "- Stdv=" << _valarg[1];
      break;

    case ETLaw::E_EXPONENTIAL:
      sstr << " Exponential - Mean=" << _valarg[0] << " - Scale=" << _valarg[1];
      break;

    case ETLaw::E_GAMMA:
      sstr << " Gamma - Mean=" << _valarg[0] << " - Scale=" << _valarg[1];
      break;

    case ETLaw::E_STABLE:
      sstr << " Stable - Alpha=" << _valarg[0] << " - Beta=" << _valarg[1]
           << " - Gamma=" << _valarg[2] << " - Delta=" << _valarg[3];
      break;

    case ETLaw::E_BETA1:
      sstr << " Beta1 - Par1=" << _valarg[0] << " - Par2=" << _valarg[1];
      break;

    case ETLaw::E_BETA2:
      sstr << " Beta2 - Par1=" << _valarg[0] << " - Par2=" << _valarg[1];
      break;
    }
  return sstr.str();
}
