/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Simulation/SimuBooleanParam.hpp"

SimuBooleanParam::SimuBooleanParam(int maxiter,
                                   double tmax,
                                   double background,
                                   double facies,
                                   const VectorDouble& dilate)
    : AStringable(),
      _maxiter(maxiter),
      _tmax(tmax),
      _background(background),
      _facies(facies),
      _dilate(dilate)
{
}

SimuBooleanParam::SimuBooleanParam(const SimuBooleanParam &r)
: AStringable(r),
   _maxiter(r._maxiter),
   _tmax(r._tmax),
   _background(r._background),
   _facies(r._facies),
   _dilate(r._dilate)
{
}

SimuBooleanParam& SimuBooleanParam::operator=(const SimuBooleanParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _maxiter = r._maxiter;
    _tmax = r._tmax;
    _background = r._background;
    _facies = r._facies;
    _dilate = r._dilate;
  }
  return *this;
}

SimuBooleanParam::~SimuBooleanParam()
{
}

int    _maxiter;
double _tmax;
double _background;
double _facies;
VectorDouble _dilate;


String SimuBooleanParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "- Maximum Number of iterations = " << _maxiter << std::endl;
  sstr << "- Maximum time available = " << _tmax << std::endl;
  sstr << "- Value assigned to the Background = " << _background << std::endl;
  sstr << "- Value assigned to the Facies     = " << _facies << std::endl;
  if (! _dilate.empty())
    sstr << toVector("Dilation", _dilate);
  return sstr.str();
}

double SimuBooleanParam::getDilate(int idim) const
{
  if (_dilate.empty()) return 0.;
  return _dilate[idim];
}
