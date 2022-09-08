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
#include "Model/Option_AutoFit.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

Option_AutoFit::Option_AutoFit()
    : AStringable(),
      _verbose(false),
      _wmode(2),
      _maxiter(1000),
      _flag_intrinsic(0),
      _tolstop(1.e-6),
      _tolred(1.e-6),
      _epsdelta(1.e-5),
      _tolsigma(5.),
      _initdelta(1.)
{
}

Option_AutoFit::Option_AutoFit(const Option_AutoFit &m)
    : AStringable(m),
      _verbose(m._verbose),
      _wmode(m._wmode),
      _maxiter(m._maxiter),
      _flag_intrinsic(m._flag_intrinsic),
      _tolstop(m._tolstop),
      _tolred(m._tolred),
      _epsdelta(m._epsdelta),
      _tolsigma(m._tolsigma),
      _initdelta(m._initdelta)
{

}

Option_AutoFit& Option_AutoFit::operator=(const Option_AutoFit &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _verbose = m._verbose;
    _wmode = m._wmode;
    _maxiter = m._maxiter;
    _flag_intrinsic = m._flag_intrinsic;
    _tolstop = m._tolstop;
    _tolred = m._tolred;
    _epsdelta = m._epsdelta;
    _tolsigma = m._tolsigma;
    _initdelta = m._initdelta;
  }
  return *this;
}

Option_AutoFit::~Option_AutoFit()
{

}

String Option_AutoFit::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(1,"Optimization parameters");

  if (getVerbose()) sstr << "Verbose option is switched ON"    << std::endl;
  sstr << "- Optimization weighting mode       " << getWmode()     << std::endl;
  sstr << "- Maximum number of iterations      " << getMaxiter()   << std::endl;
  sstr << "- Stopping criterion                " << getTolstop()   << std::endl;
  sstr << "- Stopping Criterion (scaled)       " << getTolred()    << std::endl;
  sstr << "- Minimum increment value           " << getEpsdelta()  << std::endl;
  sstr << "- Variance percentage for stripping " << getTolsigma()  << std::endl;
  sstr << "- Initial increment value           " << getInitdelta() << std::endl;
  if (getFlagIntrinsic())
    sstr << "- Resulting Multivariate Model should be Intrinsic"   << std::endl;

  return sstr.str();
}
