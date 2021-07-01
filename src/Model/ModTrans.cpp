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
#include "geoslib_f.h"
#include "Basic/Utilities.hpp"
#include "Model/ModTrans.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"

ModTrans::ModTrans()
  : _modTransMode(MODEL_PROPERTY_NONE)
  , _conv(nullptr)
  , _tape(nullptr)
  , _anam(nullptr)
  , _anamIClass(0)
  , _anamNClass(0)
  , _anamPointBlock(0)
  , _anamStrCount()
  , _anamMeans()
{
}

ModTrans::ModTrans(const ModTrans &m)
    : _modTransMode(m._modTransMode),
      _conv(m._conv),
      _tape(m._tape),
      _anam(m._anam),
      _anamIClass(m._anamIClass),
      _anamNClass(m._anamNClass),
      _anamPointBlock(m._anamPointBlock),
      _anamStrCount(m._anamStrCount),
      _anamMeans(m._anamMeans)
{
}

ModTrans& ModTrans::operator=(const ModTrans &m)
{
  if (this != &m)
  {
    _modTransMode = m._modTransMode;
    _conv = m._conv;
    _tape = m._tape;
    _anam = m._anam;
    _anamIClass = m._anamIClass;
    _anamNClass = m._anamNClass;
    _anamPointBlock = m._anamPointBlock;
    _anamStrCount = m._anamStrCount;
    _anamMeans = m._anamMeans;
  }
  return (*this);
}

ModTrans::~ModTrans()
{
  if (_conv != nullptr) delete _conv;
  if (_tape != nullptr) delete _tape;
  if (_anam != nullptr) delete _anam;
}

void ModTrans::cancelProperty()
{
  _modTransMode = MODEL_PROPERTY_NONE;
}

int ModTrans::addConvolution(int conv_type,
                             int conv_idir,
                             int conv_ndisc,
                             double conv_range)
{
  _conv = new Convolution;
  if (_conv->init(conv_type, conv_idir, conv_ndisc,conv_range)) return 1;
  _modTransMode = MODEL_PROPERTY_CONV;

  return 0;
 }


int ModTrans::addAnamorphosis(int anam_type,
                              int anam_nclass,
                              int anam_iclass,
                              int anam_var,
                              double anam_coefr,
                              double anam_coefs,
                              VectorDouble& anam_strcnt,
                              VectorDouble& anam_stats)
{
  /* Preliminary checks */

  if (! (anam_iclass == 0 || anam_iclass < anam_nclass))
  {
    messerr("The rank of the active factor (%d) is incorrect",anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)",
            anam_nclass-1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }

  /* Load the parameters */

  _anam->setType(anam_type);
  if (anam_type == ANAM_HERMITIAN)
  {
    _anam = new AnamHermite();
    AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(_anam);
    anam_hermite->setRCoef(anam_coefr);
  }
  else if (anam_type == ANAM_DISCRETE_IR)
  {
    _anam = new AnamDiscreteIR();
    AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(_anam);
    anam_discrete_IR->setNCut(anam_nclass);
    anam_discrete_IR->setRCoef(anam_coefr);
    anam_discrete_IR->setStats(anam_stats);
  }
  else if (anam_type == ANAM_DISCRETE_DD)
  {
    _anam = new AnamDiscreteDD();
    AnamDiscreteDD* anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(_anam);
    anam_discrete_DD->setNCut(anam_nclass);
    anam_discrete_DD->setSCoef(anam_coefs);
    anam_discrete_DD->setStats(anam_stats);
  }
  else
  {
    messerr("Unknown Anamorphosis type int Definition of Model Transformation");
    return 1;
  }

  _anamIClass = anam_iclass;
  _anamNClass = anam_nclass;
  _anamPointBlock = FFFF(anam_var) ? 0 : anam_var - 1;

  if (!anam_strcnt.empty())
  {
    _anamStrCount.resize(anam_nclass - 1);
    for (int i = 0; i < anam_nclass - 1; i++)
      _anamStrCount[i] = anam_strcnt[i];
  }

  // Define the Property

  _modTransMode = MODEL_PROPERTY_ANAM;

  return 0;
}

int ModTrans::addTapering(int tape_type,double tape_range)
{
  _tape = new Tapering();
  if (_tape->init(tape_type, tape_range)) return 1;
  _modTransMode = MODEL_PROPERTY_TAPE;

  return 0;
}

std::string ModTrans::toString(int level) const
{
  std::stringstream sstr;

  if (_modTransMode == MODEL_PROPERTY_NONE) return sstr.str();
  mestitle(1,"Additional Properties");

  switch (_modTransMode)
  {
    case MODEL_PROPERTY_CONV:
      sstr << _conv->toString(level);
      break;

    case MODEL_PROPERTY_TAPE:
      sstr << _tape->toString(level);
      break;

    case MODEL_PROPERTY_ANAM:
      sstr << _anam->toString(level);
      break;
  }
  return sstr.str();
}
