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
#include "Fractures/Family.hpp"
#include "Basic/AStringable.hpp"

Family::Family(double orient,
               double dorient,
               double theta0,
               double alpha,
               double ratcst,
               double prop1,
               double prop2,
               double aterm,
               double bterm,
               double range)
    : AStringable(),
      _orient(orient),
      _dorient(dorient),
      _theta0(theta0),
      _alpha(alpha),
      _ratcst(ratcst),
      _prop1(prop1),
      _prop2(prop2),
      _aterm(aterm),
      _bterm(bterm),
      _range(range)
{
}

Family::Family(const Family& r)
    : AStringable(r),
      _orient(r._orient),
      _dorient(r._dorient),
      _theta0(r._theta0),
      _alpha(r._alpha),
      _ratcst(r._ratcst),
      _prop1(r._prop1),
      _prop2(r._prop2),
      _aterm(r._aterm),
      _bterm(r._bterm),
      _range(r._range)
{
}

Family& Family::operator=(const Family& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _orient = r._orient;
    _dorient =r ._dorient;
    _theta0 = r._theta0;
    _alpha = r._alpha;
    _ratcst = r._ratcst;
    _prop1 = r._prop1;
    _prop2 = r._prop2;
    _aterm = r._aterm;
    _bterm = r._bterm;
    _range = r._range;
  }
  return *this;
}

Family::~Family()
{
}

String Family::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Average Fault Orientation       = " << _orient << "(deg)" << std::endl;
  sstr << "Tolerance for Orientation       = " << _dorient << "(deg)" << std::endl;
  sstr << "Reference Poisson Intensity     = " << _theta0 << std::endl;
  sstr << "Intensity from thick. exponent  = " << _alpha << std::endl;
  sstr << "Intensity Constant/Shaped ratio = " << _ratcst << std::endl;
  sstr << "Survival constant probability   = " << _prop1 << std::endl;
  sstr << "Survival length-dependent proba = " << _prop2 << std::endl;
  sstr << "Survival cumul. length exponent = " << _aterm << std::endl;
  sstr << "Survival thickness exponent     = " << _bterm << std::endl;
  sstr << "Fracture repulsion area Range   = " << _range << std::endl;

  return sstr.str();
}

