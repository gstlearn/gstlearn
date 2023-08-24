/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Boolean/ModelBoolean.hpp"
#include "Boolean/AShape.hpp"
#include "Basic/Law.hpp"

ModelBoolean::ModelBoolean(double thetaCst, bool flagStat)
    : AStringable(),
      _flagStat(flagStat),
      _thetaCst(thetaCst),
      _shapes()
{
}

ModelBoolean::ModelBoolean(const ModelBoolean &r)
    : AStringable(r),
      _flagStat(r._flagStat),
      _thetaCst(r._thetaCst),
      _shapes(r._shapes)
{

}

ModelBoolean& ModelBoolean::operator=(const ModelBoolean &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _flagStat = r._flagStat;
    _thetaCst = r._thetaCst;
    _shapes = r._shapes;
  }
  return *this;
}

ModelBoolean::~ModelBoolean()
{
  for (int itok = 0; itok < (int) _shapes.size(); itok++)
    delete _shapes[itok];
}

void ModelBoolean::addToken(const AShape& token)
{
  _shapes.push_back(dynamic_cast<AShape*>(token.clone()));
}

/****************************************************************************/
/*!
 **  Normalize the proportions
 **
 *****************************************************************************/
void ModelBoolean::normalizeProportions()

{
  int nb_tokens = (int) _shapes.size();
  double total = 0.;
  for (int itok = 0; itok < nb_tokens; itok++)
    total += _shapes[itok]->getProportion();

  if (ABS(total) <= 0.)
  {
    for (int itok = 0; itok < nb_tokens; itok++)
      _shapes[itok]->setProportion(1. / (double) nb_tokens);
  }
  else
  {
    for (int itok = 0; itok < nb_tokens; itok++)
      _shapes[itok]->setProportion( _shapes[itok]->getProportion() / total);
  }
  return;
}

BooleanObject* ModelBoolean::generateObject(int ndim) const
{
  int nb_token = (int) _shapes.size();

  /* Calculate the total probability */

  double total = 0.;
  for (int itok = 0; itok < nb_token; itok++)
    total += _shapes[itok]->getProportion();
  if (total <= 0.) return nullptr;

  /* Find the type of token to be generated */

  double value = total * law_uniform(0., 1.);
  int rank = -1;
  double cumul = 0.;
  for (int itok = 0; itok < nb_token; itok++)
  {
    cumul += _shapes[itok]->getProportion();
    rank = itok;
    if (value < cumul) break;
  }
  if (rank < 0) rank = nb_token - 1;
  return _shapes[rank]->generateObject(ndim);
}

String ModelBoolean::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNbTokens() <= 0) return sstr.str();

  sstr << toTitle(0, "Object Model");
  if (_flagStat)
     sstr << "- Poisson Intensity = "<< _thetaCst << std::endl;
   else
     sstr << "- Variable Poisson Intensity" << std::endl;

  for (int itok = 0; itok < getNbTokens(); itok++)
  {
    sstr << toTitle(1, "Token %d", itok+1);
    sstr << _shapes[itok]->toString(strfmt);
  }
  return sstr.str();
}
