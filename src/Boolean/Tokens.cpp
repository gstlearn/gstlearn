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
#include "Boolean/Tokens.hpp"
#include "Boolean/AToken.hpp"
#include "Basic/Law.hpp"

Tokens::Tokens()
    : _tokens()
{
}

Tokens::Tokens(const Tokens &r)
    : _tokens(r._tokens)
{

}

Tokens& Tokens::operator=(const Tokens &r)
{
  if (this != &r)
  {
    _tokens = r._tokens;
  }
  return *this;
}

Tokens::~Tokens()
{
  for (int itok = 0; itok < (int) _tokens.size(); itok++)
    delete _tokens[itok];
}

void Tokens::addToken(const AToken& token)
{
  _tokens.push_back((AToken*)token.clone());
}

/****************************************************************************/
/*!
 **  Normalize the proportions
 **
 ** \param[in]  tokens Tokens structure
 **
 *****************************************************************************/
void Tokens::normalizeProportions()

{
  int nb_tokens = (int) _tokens.size();
  double total = 0.;
  for (int itok = 0; itok < nb_tokens; itok++)
    total += _tokens[itok]->getProportion();

  if (ABS(total) <= 0.)
  {
    for (int itok = 0; itok < nb_tokens; itok++)
      _tokens[itok]->setProportion(1. / (double) nb_tokens);
  }
  else
  {
    for (int itok = 0; itok < nb_tokens; itok++)
      _tokens[itok]->setProportion( _tokens[itok]->getProportion() / total);
  }
  return;
}

Object* Tokens::generateObject(int ndim) const
{
  int nb_token = (int) _tokens.size();

  /* Calculate the total probability */

  double total = 0.;
  for (int itok = 0; itok < nb_token; itok++)
    total += _tokens[itok]->getProportion();
  if (total <= 0.) return nullptr;

  /* Find the type of token to be generated */

  double value = total * law_uniform(0., 1.);
  int rank = -1;
  double cumul = 0.;
  for (int itok = 0; itok < nb_token; itok++)
  {
    cumul += _tokens[itok]->getProportion();
    rank = itok;
    if (value < cumul) break;
  }
  if (rank < 0) rank = nb_token - 1;
  return _tokens[rank]->generateObject(ndim);
}
