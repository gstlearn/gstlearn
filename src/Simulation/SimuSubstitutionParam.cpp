/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Basic/AStringable.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"

#include <math.h>

#define TRANS(i,j)     (_trans[(j) + _nfacies * (i)])

SimuSubstitutionParam::SimuSubstitutionParam(int nfacies,
                                             double intensity,
                                             bool flag_direct,
                                             bool flag_coding,
                                             bool flag_orient)
    : AStringable(),
      _nfacies(nfacies),
      _nstates(0),
      _colfac(-1),
      _flagDirect(flag_direct),
      _flagCoding(flag_coding),
      _flagOrient(flag_orient),
      _flagAuto(true),
      _intensity(intensity),
      _factor(0.),
      _colang(),
      _vector(),
      _trans()
{
  (void) isValid();
}

SimuSubstitutionParam::SimuSubstitutionParam(const SimuSubstitutionParam &r)
    : AStringable(r),
      _nfacies(r._nfacies),
      _nstates(r._nstates),
      _colfac(r._colfac),
      _flagDirect(r._flagDirect),
      _flagCoding(r._flagCoding),
      _flagOrient(r._flagOrient),
      _flagAuto(r._flagAuto),
      _intensity(r._intensity),
      _factor(r._factor),
      _colang(r._colang),
      _vector(r._vector),
      _trans(r._trans)
{
}

SimuSubstitutionParam& SimuSubstitutionParam::operator=(const SimuSubstitutionParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _nfacies = r._nfacies;
    _nstates = r._nstates;
    _colfac = r._colfac;
    _flagDirect = r._flagDirect;
    _flagCoding = r._flagCoding;
    _flagOrient = r._flagOrient;
    _flagAuto = r._flagAuto;
    _intensity = r._intensity;
    _factor = r._factor;
    _colang = r._colang;
    _vector = r._vector;
    _trans = r._trans;
  }
  return *this;
}

SimuSubstitutionParam::~SimuSubstitutionParam()
{
}

String SimuSubstitutionParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Number of Facies = " << _nfacies << std::endl;
  if (! _flagAuto)
    sstr << "Number of States = " << _nstates << std::endl;
  sstr << "Intensity of Poisson Point Process = " << _intensity << std::endl;
  if (_flagDirect)
    sstr << "Direction information performed Internally" << std::endl;
  else
    sstr << "Direction information provided in the Db" << std::endl;
  if (_flagCoding)
    sstr << "Coding process performed internally" << std::endl;
  else
    sstr << "Coding not performed: Result is the Direction information" << std::endl;
  if (_flagOrient)
    sstr << toVector("Vector orthogonal to desorientation layering", _vector);
  sstr << "Factor for desorientation strength (0: isotropic; 1: stratified) = " <<
      _factor << std::endl;
  sstr << toVector("Transition probability matrix", _trans);
  if (_colfac >= 0)
    sstr << "Attribute rank for desorientation factor = " << _colfac << std::endl;
  if (! _colang.empty())
    sstr << toVector("Attribute ranks for Desorientation Vector", _colang);

  return sstr.str();
}

bool SimuSubstitutionParam::isValid(bool verbose)
{
  // Check the desorientation information

  if (isFlagOrient())
  {

    /* Check the (constant) angle */

    isValidOrientation(_vector, verbose);

    /* Check the (constant) desorientation factor */

    if (getColfac() < 0)
    {
      isValidFactor(&_factor, verbose);
    }
  }

  // Check the transition information

  if (_trans.empty())
  {
    _trans = VectorDouble(_nfacies * _nfacies, 1. / _nfacies);
  }
  else
  {
    if (! _isValidTransition(verbose))
      _trans = VectorDouble(_nfacies * _nfacies, 1. / _nfacies);
  }

  // Check the irreductibility

  if (isFlagCoding())
  {
    if (! _isIrreductibility(verbose)) return false;
  }

  return true;
}

/****************************************************************************/
/*!
 **  Check if the transition matrix is irreductible
 **
 ** \return  1 if the transition matrix is not irreductible; 0 otherwise
 **
 ** \param[in]  verbose  Verbose option
 **
 *****************************************************************************/
bool SimuSubstitutionParam::_isIrreductibility(bool verbose)
{

  /* Check that the transition matrix is correct */

  for (int i = 0; i < _nfacies; i++)
  {
    double total = 0.;
    for (int j = 0; j < _nfacies; j++)
    {
      if (TRANS(i,j) < 0. || TRANS(i,j) > 1.) return false;
      total += TRANS(i, j);
    }
    if (total <= 0.) return false;
    for (int j = 0; j < _nfacies; j++)
      TRANS(i,j) /= total;
  }

  /* Check the irreductibility */

  VectorInt flag(_nfacies);
  flag[0] = 0;
  int nend = 0;
  int ndeb = 0;
  for (int i = 1; i < _nfacies; i++)
  {
    flag[i] = 0;
    if (TRANS(i,0) > 0)
    {
      flag[i] = 1;
      nend++;
    }
  }

  while (ndeb != nend)
  {
    for (int i = 0; i < _nfacies; i++)
      if (flag[i])
        for (int j = 0; j < _nfacies; j++)
          if (i != j && TRANS(j,i) > 0) flag[j] = 1;
    ndeb = nend;
    for (int i = nend = 0; i < _nfacies; i++)
      nend += flag[i];
  }
  if (nend != _nfacies) return false;

  /* Printout (conditional) */

  if (verbose)
    print_matrix("Transitions", 0, 1, _nfacies, _nfacies, NULL, _trans.data());

  return true;
}

/*****************************************************************************
 **
 ** Checks the validity of an orientation vector
 **
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
void SimuSubstitutionParam::isValidOrientation(VectorDouble& vector,
                                               bool verbose) const
{
  int ndim = (int) vector.size();
  double total = 0.;
  for (int i = 0; i < ndim; i++)
    total += vector[i] * vector[i];
  if (total <= 0.)
  {
    if (verbose)
    {
      messerr("The desorientation vector should not be zero");
      messerr("It is set to the first Direction Unit vector");
    }
    vector[0] = 1.;
    total = 1.;
  }
  for (int i = 0; i < ndim; i++)
    vector[i] /= sqrt(total);
}

/*****************************************************************************
 **
 ** Checks the validity of an orientation factor
 **
 ** \param[in,out] factor   Input and output factor values
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
void SimuSubstitutionParam::isValidFactor(double* factor, bool verbose) const
{
  if (*factor < 0.)
  {
    if (verbose)
    {
      messerr("The desorientation factor cannot be negative");
      messerr("It is set to 0.");
    }
    *factor = 0.;
  }
  if (*factor > 1.)
  {
    if (verbose)
    {
      messerr("The desorientation factor cannot be larger than 1");
      messerr("It is set to 1.");
    }
    *factor = 1.;
  }
  return;
}

bool SimuSubstitutionParam::isAngleLocal() const
{
  if (_colang.empty()) return false;
  for (int i = 0; i < (int) _colang.size(); i++)
    if (_colang[i] >= 0) return true;
  return false;
}

bool SimuSubstitutionParam::isLocal() const
{
  return isAngleLocal() || _colfac >= 0;
}

bool SimuSubstitutionParam::_isValidTransition(bool verbose,double eps)
{
  if ((int) _trans.size() != _nfacies * _nfacies) return false;

  for (int irow = 0; irow < _nfacies; irow++)
  {
    double total = 0.;
    for (int icol = 0; icol < _nfacies; icol++)
    {
      total += TRANS(irow, icol);
    }
    if (ABS(total - 1.) > eps)
    {
      if (verbose)
        messerr("Transition: Sum of elements of row(%d) must be 1 (%lf)",
                irow+1,total);
      return false;
    }
  }
  return true;
}

int SimuSubstitutionParam::getColang(int idim) const
{
  if (idim < (int) _colang.size())
    return _colang[idim];
  else
    return 0.;
}
