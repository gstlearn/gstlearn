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
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Basic/Vector.hpp"

#define TRANS(i,j)     (_trans[(j) + _nfacies * (i)])

SimuSubstitutionParam::SimuSubstitutionParam()
    : AStringable(),
      _nfacies(0),
      _nstates(0),
      _colfac(0),
      _flagDirect(false),
      _flagCoding(false),
      _flagOrient(false),
      _flagAuto(false),
      _intensity(1.),
      _factor(0.),
      _vector(),
      _colang(),
      _trans()
{
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
      _vector(r._vector),
      _colang(r._colang),
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
    _vector = r._vector;
    _colang = r._colang;
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
  return sstr.str();
}

bool SimuSubstitutionParam::isValid(bool verbose) const
{
  bool flag_angloc = false;
  if (isFlagCoding())
  {
    if (! _isIrreductibility(verbose)) return false;
  }

  /* Check that the validity of the desorientation information */

  if (isFlagOrient())
  {
    for (int i = 0; i < 3; i++)
      if (getColang(i) >= 0) flag_angloc = true;

    /* Check the (constant) angle */

    if (!flag_angloc && st_check_orientation(getVector(), 1))
      goto label_end;

    /* Check the (constant) desorientation factor */

    if (subparam.getColfac() < 0 && st_check_factor(&subparam.getFactor(), 1))
      goto label_end;

    flag_local = (flag_angloc || subparam.getColfac() >= 0);
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
 ** \returns 1 if the vector is not valid; 0 otherwise
 **
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
int SimuSubstitutionParam::_check_orientation(bool verbose)
{
  double total = 0.;
  for (int i = 0; i < 3; i++)
    total += _vector[i] * _vector[i];
  if (total <= 0.)
  {
    if (verbose)
    {
      messerr("The desorientation vector should not be zero");
      return (1);
    }
    _vector[0] = 1.;
    total = 1.;
  }
  for (i = 0; i < 3; i++)
    _vector[i] /= sqrt(total);
  return (0);
}
