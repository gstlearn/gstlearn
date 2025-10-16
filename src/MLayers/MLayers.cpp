/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "MLayers/MLayers.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Model/Model.hpp"
#include "geoslib_old_f.h"

#define IAD(n, i, j)      ((n) * (i) + (j))
#define INVS(_npar, i, j) (invS[IAD(_npar, i, j)])

namespace gstlrn
{
MLayers::MLayers()
  : AStringable()
  , _flagSame(false)
  , _flagVel(false)
  , _flagCumul(true)
  , _flagExt(false)
  , _flagZ(true)
  , _colrefD(-1)
  , _colrefT(-1)
  , _colrefB(-1)
  , _matchTime(false)
  , _ptime(ELoc::Z)
  , _nlayers(0)
  , _nbfl(0)
  , _nech(0)
  , _npar(0)
  , _neq(0)
  , _dbin()
  , _dbout()
  , _model()
{
}

MLayers::MLayers(const MLayers& m)
  : AStringable(m)
  , _flagSame(m._flagSame)
  , _flagVel(m._flagVel)
  , _flagCumul(m._flagCumul)
  , _flagExt(m._flagExt)
  , _flagZ(m._flagZ)
  , _colrefD(m._colrefD)
  , _colrefT(m._colrefT)
  , _colrefB(m._colrefB)
  , _matchTime(m._matchTime)
  , _ptime(m._ptime)
  , _nlayers(m._nlayers)
  , _nbfl(m._nbfl)
  , _nech(m._nech)
  , _npar(m._npar)
  , _neq(m._neq)
  , _dbin(m._dbin)
  , _dbout(m._dbout)
  , _model(m._model)
{
}

MLayers& MLayers::operator=(const MLayers& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _flagSame  = m._flagSame;
    _flagVel   = m._flagVel;
    _flagCumul = m._flagCumul;
    _flagExt   = m._flagExt;
    _flagZ     = m._flagZ;
    _colrefD   = m._colrefD;
    _colrefT   = m._colrefT;
    _colrefB   = m._colrefB;
    _matchTime = m._matchTime;
    _ptime     = m._ptime;
    _nlayers   = m._nlayers;
    _nbfl      = m._nbfl;
    _nech      = m._nech;
    _npar      = m._npar;
    _neq       = m._neq;
    _dbin      = m._dbin;
    _dbout     = m._dbout;
    _model     = m._model;
  }
  return *this;
}

MLayers::~MLayers()
{
}

String MLayers::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  static const char* NOK[] = {"NO", "YES"};

  std::stringstream sstr;
  sstr << "List of Apices" << std::endl;

  sstr << toTitle(0, "Multi-Layers Environments") << std::endl;
  if (_flagVel)
    sstr << "Working in Velocity" << std::endl;
  else
    sstr << "Working in Depth" << std::endl;
  if (_flagCumul)
    sstr << "Producing estimation in Depth" << std::endl;
  else
    sstr << "Producing estimation in Thickness" << std::endl;
  if (_flagZ)
    sstr << "Results are converted into Depth" << std::endl;

  sstr << "Do the Input and Output Db coincide: " << NOK[_flagSame] << std::endl;
  sstr << "Using External Drift functions: " << NOK[_flagExt] << std::endl;
  sstr << "Is Time used as External Drift: %s" << NOK[_matchTime] << std::endl;
  if (_colrefD >= 0)
    sstr << "Rank of the Reference Depth Map = " << _colrefD << std::endl;
  if (_colrefT >= 0)
    sstr << "Rank of the Reference Time Map = " << _colrefT << std::endl;
  if (_colrefB >= 0)
    sstr << "Rank of the Bottom Depth Map = " << _colrefB << std::endl;

  sstr << "Number of layers = " << _nlayers << std::endl;
  sstr << "Number of drift functions (per layer) = " << _nbfl << std::endl;
  sstr << "Number of active samples (including collocated duplicates) = " << _nech << std::endl;

  return sstr.str();
}

/****************************************************************************/
/*!
 **  Returns the number of drift functions
 **
 ** \return  Number of drift conditions
 **
 *****************************************************************************/
Id MLayers::_getNDrift() const
{
  switch (_irfRank)
  {
    case -1:
      return (0);
      break;

    case 0:
      if (!_flagExt)
        return (1);
      else
        return (2);
      break;

    case 1:
      if (!_flagExt)
        return (3);
      else
        return (4);
      break;

    default:
      messageAbort("_irfRank must be -1, 0 or 1");
      break;
  }
  return (0);
}

void MLayers::setGeneral(Db* dbin,
                         DbGrid* dbout,
                         Model* model,
                         bool flag_std,
                         bool flag_same,
                         bool flag_vel,
                         bool flag_cumul,
                         bool flag_ext,
                         bool flag_Z,
                         bool match_time,
                         Id irf_rank,
                         const VectorDouble& prior_means,
                         const VectorDouble& prior_vars,
                         const String& namerefd,
                         const String& namereft,
                         const String& namerefb)
{
  _dbin      = dbin;
  _dbout     = dbout;
  _model     = model;
  Id colrefd = (namerefd.empty()) ? -1 : _dbout->getColIdx(namerefd);
  Id colreft = (namereft.empty()) ? -1 : _dbout->getColIdx(namereft);
  Id colrefb = (namerefb.empty()) ? -1 : _dbout->getColIdx(namerefb);
  _flagStd   = flag_std;
  _irfRank   = irf_rank;
  _nlayers   = _model->getNVar();
  _flagSame  = flag_same;
  _flagVel   = flag_vel;
  _flagCumul = flag_cumul;
  _flagExt   = flag_ext;
  _flagZ     = flag_Z;
  _colrefD   = colrefd;
  _colrefT   = colreft;
  _colrefB   = colrefb;
  _matchTime = match_time;
  _ptime     = (match_time) ? ELoc::F : ELoc::TIME;
  _nbfl      = _getNDrift();
  _nechmax   = _dbin->getNSample();
  _nech      = 0;
  _npar      = _nbfl * _nlayers;
  _neq       = _nechmax + _npar;
  _priorMean = prior_means;
  _priorVars = prior_vars;
}

/****************************************************************************/
/*!
 **  Returns the absolute grid node absolute index which is the closest to a
 **  given sample of a Db
 **  In the case of same input and output file, simply return 'iech'
 **
 ** \return The rank or ITEST
 **
 ** \param[in]  iech      Rank in the input Db
 **
 *****************************************************************************/
Id MLayers::_locateSampleInDbout(Id iech) const
{
  /* In the case the input and output files coincide, simply return 'iech' */
  if (_flagSame) return iech;

  Id ndim = _dbout->getNDim();
  VectorInt indg(ndim, 0);
  VectorDouble coor(ndim);

  /* The files are different */
  for (Id idim = 0; idim < _dbin->getNDim(); idim++)
    coor[idim] = _dbin->getCoordinate(iech, idim);
  if (point_to_grid(_dbout, coor.data(), 0, indg.data()) != 0) return ITEST;
  return _dbout->indiceToRank(indg);
}

/****************************************************************************/
/*!
 **  Check if the target layer rank is consistent
 **
 ** \param[in]  string    Name of the calling procedure
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \remarks If this target layer rank is not correct, mes_abort() is called
 ** \remarks and the program is interrupted as this must never happen.
 **
 *****************************************************************************/
void MLayers::_isValidLayer(const char* string, Id ilayer0) const
{
  if (ilayer0 >= 1 && ilayer0 <= _nlayers) return;

  messerr("Error when calling function %s", string);
  messerr("- Number of layers         = %d", _nlayers);
  messerr("- Rank of the target layer = %d", ilayer0);
  messageAbort("This error should never happen");
}

/****************************************************************************/
/*!
 **  Fill the proportion vector at a output location, up to the target layer
 **
 ** \returns 1 if the proportion vector cannot be defined; 0 otherwise
 **
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \param[out] props     Working array (Dimension: nlayers)
 **
 *****************************************************************************/
Id MLayers::_getPropsResult(Id iech,
                            Id ilayer0,
                            VectorDouble& props) const
{
  double tlast;
  _isValidLayer("_getPropsResult", ilayer0);
  for (Id ilayer = 0; ilayer < _nlayers; ilayer++)
    props[ilayer] = 0.;

  /* Dispatch */

  if (_flagVel)
  {

    /* Working in velocities */

    double t0 = (_colrefT >= 0) ? _dbout->getArray(iech, _colrefT) : 0.;
    if (FFFF(t0)) return (1);
    double t1 = _dbout->getFromLocator(_ptime, iech, ilayer0 - 1);
    if (FFFF(t1)) return (1);
    tlast = t0;

    /* Loop on the layers */

    for (Id ilayer = 0; ilayer < ilayer0; ilayer++)
    {
      double tt = _dbout->getFromLocator(_ptime, iech, ilayer);
      if (FFFF(tt)) return (1);
      double pval = (tt - tlast) / (t1 - t0);
      if (pval < 0 || pval > 1) return (1);
      tlast         = tt;
      props[ilayer] = pval;
    }
  }
  else
  {

    /* Working in depth */

    for (Id ilayer = 0; ilayer < ilayer0; ilayer++)
      props[ilayer] = 1.;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Fill the proportion vector at a data location, up to the target layer
 **
 ** \returns 1 if the proportion vector cannot be defined; 0 otherwise
 **
 ** \param[in]  iech      Rank of the target sample
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \param[out] props     Working array (Dimension: nlayers)
 **
 *****************************************************************************/
Id MLayers::_getPropsData(Id iech,
                          Id ilayer0,
                          VectorDouble& props) const
{
  _isValidLayer("_getPropsData", ilayer0);
  props.fill(0.);

  // Get the sample rank in the output Db of the sample from the input Db //
  Id igrid = _locateSampleInDbout(iech);
  if (IFFFF(igrid)) return (1);

  // Evaluate the proportion vector //
  if (_getPropsResult(igrid, ilayer0, props)) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Return the external drift value at an output location for a target layer
 **
 ** \returns  The external drift value of TEST
 **
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (Starting from 1)
 **
 *****************************************************************************/
double MLayers::_getDriftResult(Id iech, Id ilayer0)
{
  if (!_flagExt) return (TEST);
  _isValidLayer("_getDriftResult", ilayer0);
  return _dbout->getLocVariable(ELoc::F, iech, ilayer0 - 1);
}

/****************************************************************************/
/*!
 **  Return the external drift value at an input location for a target layer
 **
 ** \returns  The external drift value of TEST
 **
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (Starting from 1)
 **
 *****************************************************************************/
double MLayers::_getDriftData(Id iech, Id ilayer0)
{
  if (!_flagExt) return (TEST);
  _isValidLayer("_getDriftData", ilayer0);
  Id igrid = _locateSampleInDbout(iech);
  if (IFFFF(igrid)) return (TEST);
  return _getDriftResult(igrid, ilayer0);
}

/****************************************************************************/
/*!
 **  Calculate the array of covariances for zero distance
 **
 ** \param[in]  prop1     Working array at first point (Dimension: nlayers)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] c00       Returned array (Dimension = nlayers)
 **
 ** \remarks:  This array depends on the target location through proportions
 **
 *****************************************************************************/
void MLayers::_getC00(const VectorDouble& prop1,
                      MatrixSquare& covtab,
                      VectorDouble& c00)
{
  _model->evaluateMatInPlace(nullptr, VectorDouble(), covtab, true);

  if (_flagCumul)
  {
    for (Id k = 0; k < _nlayers; k++)
    {
      double value        = 0.;
      bool flag_interrupt = false;
      for (Id i = 0; i <= k && !flag_interrupt; i++)
        for (Id j = 0; j <= k && !flag_interrupt; j++)
        {
          if (FFFF(prop1[i]) || FFFF(prop1[j]))
            flag_interrupt = true;
          else
            value += prop1[i] * prop1[j] * covtab.getValue(i, j);
        }
      c00[k] = (flag_interrupt) ? TEST : value;
    }
  }
  else
  {
    for (Id k = 0; k < _nlayers; k++)
      c00[k] = covtab.getValue(k, k);
  }
}

/****************************************************************************/
/*!
 **  Calculate the covariance between data and data
 **
 ** \return  The covariance terms or TEST
 **
 ** \param[in]  ilayer    Layer of interest (first point). Starting from 1
 ** \param[in]  prop1     Working array at first point (Dimension: nlayers)
 ** \param[in]  jlayer    Layer of interest (second point). Starting from 1
 ** \param[in]  prop2     Working array at second point (Dimension: nlayers)
 ** \param[in]  dd        Distance vector (or NULL for zero-distance)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 **
 ** \remarks:  As this function may return TEST, TEST value should be tested
 **
 *****************************************************************************/
double MLayers::_getCIJ(Id ilayer,
                        const VectorDouble& prop1,
                        Id jlayer,
                        const VectorDouble& prop2,
                        const VectorDouble& dd,
                        MatrixSquare& covtab)
{
  VectorDouble d1(2);
  _isValidLayer("_covarianceCIJ", ilayer);
  _isValidLayer("_covarianceCIJ", jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (!dd.empty()) ? dd[0] : 0.;
  d1[1] = (!dd.empty()) ? dd[1] : 0.;
  _model->evaluateMatInPlace(nullptr, d1, covtab, true);

  /* Evaluate the covariance term */

  double value = 0.;
  for (Id i = 0; i < ilayer; i++)
    for (Id j = 0; j < jlayer; j++)
    {
      if (FFFF(prop1[i])) return (TEST);
      if (FFFF(prop2[j])) return (TEST);
      value += prop1[i] * prop2[j] * covtab.getValue(i, j);
    }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the covariance between data and target
 **
 ** \return  The covariance terms or TEST
 **
 ** \param[in]  ilayer    Layer of interest (data point). Starting from 1
 ** \param[in]  prop1     Working array at data point (Dimension: nlayers)
 ** \param[in]  jlayer    Layer of interest (target point). Starting from 1
 ** \param[in]  dd        Distance vector (or NULL for zero-distance)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 **
 ** \remarks:  As this function may return TEST, TEST value should be tested
 **
 *****************************************************************************/
double MLayers::_getCI0(Id ilayer,
                        const VectorDouble& prop1,
                        Id jlayer,
                        const VectorDouble& dd,
                        MatrixSquare& covtab)
{
  VectorDouble d1(2);
  _isValidLayer("_getCI0", ilayer);
  _isValidLayer("_getCI0", jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (!dd.empty()) ? dd[0] : 0.;
  d1[1] = (!dd.empty()) ? dd[1] : 0.;
  _model->evaluateMatInPlace(nullptr, d1, covtab, true);

  /* Evaluate the covariance term */

  double value = 0.;
  for (Id i = 0; i < ilayer; i++)
  {
    if (FFFF(prop1[i])) return (1);
    value += prop1[i] * covtab.getValue(i, jlayer - 1);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the drift terms
 **
 ** \return  Error return code
 **
 ** \param[in]  coor        Array of coordinates
 ** \param[in]  propval     Value for the proportion (used if flag_cumul)
 ** \param[in]  drext       Value of the external drift
 ** \param[in,out] ipos_loc Address for the first drift term.
 **                         On output, address for the next term after the drift
 **
 ** \param[out] b           Array for storing the drift
 **
 *****************************************************************************/
Id MLayers::_fillDrift(const VectorDouble& coor,
                       double propval,
                       double drext,
                       Id* ipos_loc,
                       VectorDouble& b) const
{
  if (_flagExt && FFFF(drext)) return (1);
  Id ipos = *ipos_loc;
  switch (_nbfl)
  {
    case 0:
      break;

    case 1:
      b[ipos++] = propval;
      break;

    case 2:
      b[ipos++] = propval;
      b[ipos++] = propval * drext;
      break;

    case 3:
      b[ipos++] = propval;
      b[ipos++] = propval * coor[0];
      b[ipos++] = propval * coor[1];
      break;

    case 4:
      b[ipos++] = propval;
      b[ipos++] = propval * coor[0];
      b[ipos++] = propval * coor[1];
      b[ipos++] = propval * drext;
      break;
  }
  *ipos_loc = ipos;
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the L.H.S. for one data
 **
 ** \returns:  1 if the L.H.S. vector cannot be calculated; 0 otherwise
 **
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  iech0     Rank of the target sample (used for ext. drift)
 ** \param[in]  ilayer0   Rank of the layer of interest (Starting from 1)
 ** \param[in]  coor      Coordinates of the data
 ** \param[in]  prop0     Working array at first point (Dimension: nlayers)
 **
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b         R.H.S. vector (Dimension = neq)
 **
 *****************************************************************************/
Id MLayers::_computeLhsOne(VectorInt& seltab,
                           Id iech0,
                           Id ilayer0,
                           VectorDouble& coor,
                           VectorDouble& prop0,
                           VectorDouble& prop2,
                           MatrixSquare& covtab,
                           VectorDouble& b)
{
  double coor2[2];
  VectorDouble d1(2);

  /* Initializations */

  b.fill(0);

  /* Covariance part */

  Id jjech = 0;
  for (Id jech = 0; jech < _dbin->getNSample(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = _dbin->getCoordinate(jech, 0);
    coor2[1] = _dbin->getCoordinate(jech, 1);
    for (Id jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
    {
      Id jlayer = (jfois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, jech)) : _nlayers;

      /* Evaluate the proportion vector */

      if (_getPropsData(jech, jlayer, prop2)) return (1);

      /* Calculate the distance vector */

      d1[0]    = coor[0] - coor2[0];
      d1[1]    = coor[1] - coor2[1];
      b[jjech] = _getCIJ(ilayer0, prop0, jlayer, prop2, d1, covtab);
      if (FFFF(b[jjech])) return (1);
    }
  }

  /* Drift part */

  for (Id i = 0; i < ilayer0; i++)
  {
    double drext = _getDriftData(iech0, i + 1);
    if (_fillDrift(coor, prop0[i], drext, &jjech, b)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the Kriging R.H.S.
 **
 ** \returns:  1 if the R.H.S. has not been calculated; 0 otherwise
 **
 ** \param[in]  coor      Coordinates of the target sample
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  iechout   Rank of the target sample
 ** \param[in]  ilayer0   Rank of the layer of interest (Starting from 1)
 **
 ** \param[out] prop0     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b         R.H.S. vector (Dimension = neq)
 **
 *****************************************************************************/
Id MLayers::_computeRHS(VectorDouble& coor,
                        VectorInt& seltab,
                        Id iechout,
                        Id ilayer0,
                        VectorDouble& prop0,
                        VectorDouble& prop2,
                        MatrixSquare& covtab,
                        VectorDouble& b)
{
  double coor2[2];
  VectorDouble d1(2);

  /* Get the coordinates of the target */

  _isValidLayer("_computeRHS", ilayer0);
  coor[0] = _dbout->getCoordinate(iechout, 0);
  coor[1] = _dbout->getCoordinate(iechout, 1);

  /* Initialize the vector with zeroes */

  for (Id i = 0; i < _neq; i++)
    b[i] = 0.;

  /* Covariance part */

  for (Id jech = 0, jjech = 0; jech < _dbin->getNSample(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = _dbin->getCoordinate(jech, 0);
    coor2[1] = _dbin->getCoordinate(jech, 1);
    for (Id ifois = 0; ifois < seltab[jech]; ifois++, jjech++)
    {
      Id jlayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, jech)) : _nlayers;

      /* Evaluate the proportion vector */

      (void)_getPropsData(jech, jlayer, prop2);

      /* Calculate the distance vector */

      d1[0] = coor2[0] - coor[0];
      d1[1] = coor2[1] - coor[1];
      if (_flagCumul)
        b[jjech] = _getCIJ(ilayer0, prop0, jlayer, prop2, d1, covtab);
      else
        b[jjech] = _getCI0(jlayer, prop2, ilayer0, d1, covtab);
      if (FFFF(b[jjech])) return (1);
    }
  }

  /* Drift part */

  Id ideb = (_flagCumul) ? 0 : ilayer0 - 1;
  for (Id i = ideb; i < ilayer0; i++)
  {
    Id ipos        = _nech + _nbfl * i;
    double drext   = _getDriftResult(iechout, i + 1);
    double propval = (_flagCumul) ? prop0[i] : 1.;
    if (_fillDrift(coor, propval, drext, &ipos, b)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the Kriging L.H.S.
 **
 ** \returns:  1 if the L.H.S. vector cannot be calculated; 0 otherwise
 **
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] prop1     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] atot      L.H.S. (square) matrix
 ** \param[out] acov      L.H.S. (square) covariance matrix
 **
 *****************************************************************************/
Id MLayers::_computeLHS(VectorInt& seltab,
                        VectorDouble& prop1,
                        VectorDouble& prop2,
                        MatrixSquare& covtab,
                        MatrixSquare& atot,
                        MatrixSquare& acov)
{
  VectorDouble coor(2);
  VectorDouble b(_neq);
  atot.fill(0);

  /* Loop on the first sample */

  for (Id iech = 0, iiech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = _dbin->getCoordinate(iech, 0);
    coor[1] = _dbin->getCoordinate(iech, 1);
    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech)) : _nlayers;

      /* Evaluate the proportion vector */

      if (_getPropsData(iech, ilayer, prop1)) return (1);

      /* Loop on the second sample */

      if (_computeLhsOne(seltab, iech, ilayer, coor,
                         prop1, prop2, covtab, b)) return (1);

      for (Id i = 0; i < _neq; i++)
        atot.setValue(iiech, i, b[i]);
    }
  }

  /* Symmetrization */
  for (Id iiech = 0; iiech < _neq; iiech++)
    for (Id jjech = 0; jjech <= iiech; jjech++)
      atot.setValue(iiech, jjech, atot.getValue(jjech, iiech));

  /* Extraction of the Covariance matrix */
  for (Id iech = 0; iech < _nech; iech++)
    for (Id jech = 0; jech < _nech; jech++)
      acov.setValue(iech, jech, atot.getValue(iech, jech));

  return (0);
}

/****************************************************************************/
/*!
 **  Establish the vector of data
 **
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] zval      The data vector (Dimension: neq)
 **
 *****************************************************************************/
void MLayers::_getDataVector(VectorInt& seltab, VectorDouble& zval)
{
  double value;

  /* Initialize the vector with zeroes */

  Id igrid = 0;
  zval.fill(0.);

  /* Loop on the samples */

  for (Id iech = 0, iiech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    /* Calculate the grid node index (optional) */

    if (_colrefD >= 0 || _colrefT >= 0 ||
        _colrefB >= 0 || _flagVel)
      igrid = _locateSampleInDbout(iech);

    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech)) : _nlayers;

      if (ifois == 0)
      {

        /* Depth of the actual sample */

        value = _dbin->getZVariable(iech, 0);
      }
      else
      {

        /* Depth of the collocated bottom sample */

        value = _dbout->getArray(igrid, _colrefB);
      }

      /* Centering to the reference Depth surface */

      if (_colrefD >= 0)
        value -= _dbout->getArray(igrid, _colrefD);

      /* Converting into velocities */

      if (_flagVel)
      {
        double dtime = _dbout->getFromLocator(_ptime, igrid, ilayer - 1);
        if (_colrefT >= 0)
          dtime -= _dbout->getArray(igrid, _colrefT);
        value /= dtime;
      }
      zval[iiech] = value;
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the Drift and subtract it from the Data
 **
 ** \param[in]  verbose   True for a  verbose option
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] zval      The data vector (Dimension: neq)
 **
 ** \remarks In the current version, the optimal coefficients of the Drift
 ** \remarks are output using the set_keypair mechanism using the keyword:
 ** \remarks "Optim_Drift_Coeffs" which returns 'ipos' values
 **
 *****************************************************************************/
Id MLayers::_subtractOptimalDrift(bool verbose,
                                  VectorInt& seltab,
                                  VectorDouble& zval)
{
  Id ipos;
  Id flag_subtract = 1;
  VectorDouble coor(2);

  /* Initializations */

  Id error = 1;
  Id neqd  = _nbfl * _nlayers;

  /* Core allocation */

  VectorDouble props(_nlayers);
  VectorDouble drift(neqd);
  VectorDouble coeff(neqd);
  MatrixSquare atab(neqd);
  VectorDouble btab(neqd);

  /* Find the vector of optimal mean values */

  for (Id iech = 0, iiech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech)) : _nlayers;

      /* Evaluate the proportion vector */

      if (_getPropsData(iech, ilayer, props)) goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = _dbin->getCoordinate(iech, 0);
      coor[1] = _dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (Id i = 0; i < _nlayers; i++)
      {
        double drext = _getDriftData(iech, i + 1);
        if (_fillDrift(coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Calculate the contribution to the different arrays */

      for (Id k1 = 0; k1 < ipos; k1++)
      {
        btab[k1] += zval[iiech] * drift[k1];
        for (Id k2 = 0; k2 < ipos; k2++)
          atab.setValue(k1, k2, atab.getValue(k1, k2) + drift[k1] * drift[k2]);
      }
    }
  }

  /* Find the optimal drift coefficients */

  if (matrix_invertFromMatrixSquare(&atab, neqd)) goto label_end;
  atab.prodMatVecInPlace(btab, coeff);

  /* Optional printout of the result */

  if (verbose)
    print_matrix("Estimated Drift", 0, 1, _nlayers, _nbfl, NULL, coeff.data());

  /* Subtract the optimal mean */

  for (Id iech = 0, iiech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech)) : _nlayers;

      /* Evaluate the proportion vector */

      if (_getPropsData(iech, ilayer, props)) goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = _dbin->getCoordinate(iech, 0);
      coor[1] = _dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (Id i = 0; i < _nlayers; i++)
      {
        double drext = _getDriftData(iech, i + 1);
        if (_fillDrift(coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Subtract the optimal estimation of the drift */

      if (flag_subtract)
        for (Id k1 = 0; k1 < ipos; k1++)
          zval[iiech] -= coeff[k1] * drift[k1];

      /* Print the residuals (optional) */

      if (OptDbg::query(EDbg::VARIOGRAM))
        message("Sample %d (Layer %d) - Coor = %lf %lf - Residual = %lf\n",
                iech + 1, ilayer, coor[0], coor[1], zval[iiech]);
    }
  }

  /* Set the error return code */
  error = 0;

label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Check if an intercept with the bottom layer is located close enough
 **  to the current sample
 **
 ** \return  1 if a duplicate must be generated; 0 otherwise
 **
 ** \param[in]  iech0     Rank of sample to be discarded (or -1)
 ** \param[in]  coor      Coordinates of the target
 ** \param[in]  eps       Tolerance
 **
 *****************************************************************************/
Id MLayers::_getCloseSample(Id iech0, const VectorDouble& coor, double eps)
{
  /* Check if a close sample has already been reviewed */

  for (Id iech = 0; iech < iech0; iech++)
  {
    double dx = _dbin->getCoordinate(iech, 0) - coor[0];
    if (ABS(dx) > eps) continue;
    double dy = _dbin->getCoordinate(iech, 1) - coor[1];
    if (ABS(dy) > eps) continue;
    return (0);
  }

  /* Check among the subsequent samples if a sample with matching coordinates */
  /* and belonging to the bottom surface exists */

  for (Id iech = iech0 + 1; iech < _dbin->getNSample(); iech++)
  {
    double dx = _dbin->getCoordinate(iech, 0) - coor[0];
    if (ABS(dx) > eps) continue;
    double dy = _dbin->getCoordinate(iech, 1) - coor[1];
    if (ABS(dy) > eps) continue;
    Id ilayer = static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech));
    if (ilayer == _nlayers) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Perform the per-calculation for estimation with collocated option
 **
 ** \param[in]  iechout   Rank of the target
 ** \param[in]  coor      Coordinates of the target
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  a         L.H.S. (square) inverted matrix
 ** \param[in]  zval      Data vector (extended)
 **
 ** \param[out] prop1     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b2        Working vector (Dimension = neq)
 ** \param[out] baux      Working vector (Dimension = neq)
 ** \param[out] ratio     Ratio value
 **
 *****************************************************************************/
Id MLayers::_prepareCollocation(Id iechout,
                                VectorDouble& coor,
                                VectorInt& seltab,
                                MatrixSquare* a,
                                VectorDouble& zval,
                                VectorDouble& prop1,
                                VectorDouble& prop2,
                                MatrixSquare& covtab,
                                VectorDouble& b2,
                                VectorDouble& baux,
                                double* ratio)
{
  (*ratio) = 0.;

  double botval = _dbout->getArray(iechout, _colrefB);
  if (_colrefD >= 0)
    botval -= _dbout->getArray(iechout, _colrefD);

  if (_getPropsResult(iechout, _nlayers, prop1)) return 1;
  double c0 = _getCIJ(_nlayers, prop1, _nlayers, prop1, VectorDouble(), covtab);
  if (FFFF(c0)) return 1;

  if (_computeLhsOne(seltab, iechout, _nlayers, coor,
                     prop1, prop2, covtab, baux)) return 1;
  a->prodMatVecInPlace(baux, b2);
  double coefz = VH::innerProductVD(b2, zval);
  double coefa = VH::innerProductVD(b2, baux);
  (*ratio)     = (ABS(c0 - coefa) > 1.e-6) ? (botval - coefz) / (c0 - coefa) : 0.;

  return 0;
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in regular case
 **
 ** \param[in]  c00       Variance for target
 ** \param[in]  a         L.H.S. (square) inverted matrix
 ** \param[in]  b         Working vector (Dimension = neq)
 ** \param[in]  dual      Dual vector
 ** \param[in]  wgt       Working array (Dimension = neq)
 **
 ** \param[out] estim     Estimated value
 ** \param[out] stdev     Standard deviation of estimation error
 **
 *****************************************************************************/
void MLayers::_estimateRegular(double c00,
                               MatrixSquare* a,
                               VectorDouble& b,
                               VectorDouble& dual,
                               VectorDouble& wgt,
                               double* estim,
                               double* stdev) const
{
  *estim = VH::innerProductVD(dual, b);

  /* Perform the variance of estimation error */

  if (_flagStd)
  {
    double c00val = c00;
    double stdv;
    if (FFFF(c00val))
      stdv = TEST;
    else
    {
      a->prodMatVecInPlace(b, wgt);
      stdv = c00val - VH::innerProductVD(b, wgt);
      stdv = (stdv > 0) ? sqrt(stdv) : 0.;
    }
    *stdev = stdv;
  }
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in the bayesian case
 **
 ** \param[in]  c00       Variance at target
 ** \param[in]  acov      L.H.S. (square) inverted matrix
 ** \param[in]  zval      Data vector
 ** \param[in]  b         Working vector (Dimension = neq)
 ** \param[in]  wgt       Working array (Dimension = neq)
 ** \param[out] post_mean Array of posterior mean
 ** \param[out] a0        Constant term
 ** \param[out] cc        Output value
 ** \param[out] ss        Output value
 ** \param[out] gs        Output value
 **
 ** \param[out] estim     Estimated value
 ** \param[out] stdev     Standard deviation of estimation error
 **
 *****************************************************************************/
void MLayers::_estimateBayes(double c00,
                             const MatrixSquare* acov,
                             VectorDouble& zval,
                             VectorDouble& b,
                             VectorDouble& wgt,
                             VectorDouble& post_mean,
                             MatrixDense& a0,
                             MatrixSquare& cc,
                             MatrixDense& ss,
                             const MatrixSquare& gs,
                             double* estim,
                             double* stdev) const
{
  VectorDouble ff0(b.begin() + _nech, b.end());

  /* Core allocation */

  VectorDouble temp(_npar);
  VectorDouble fsf0(_nech);
  VectorDouble c2(_nech);

  /* Perform the estimation */

  a0.prodMatVecInPlace(ff0, fsf0);
  for (Id iech = 0; iech < _nech; iech++)
    c2[iech] = b[iech] + fsf0[iech];
  cc.prodMatVecInPlace(c2, wgt);

  double estim1 = VH::innerProductVD(wgt, zval);
  double estim2 = VH::innerProductVD(ff0, post_mean);
  *estim        = estim1 + estim2;

  /* Calculate the standard deviation */

  if (_flagStd)
  {
    ss.prodVecMatInPlace(b, temp);
    for (Id ipar = 0; ipar < _npar; ipar++)
      temp[ipar] -= ff0[ipar];

    double stdv = c00;
    for (Id iech = 0; iech < _nech; iech++)
      for (Id jech = 0; jech < _nech; jech++)
        stdv -= b[iech] * acov->getValue(iech, jech) * b[jech];

    for (Id ipar = 0; ipar < _npar; ipar++)
      for (Id jpar = 0; jpar < _npar; jpar++)
        stdv += temp[ipar] * gs.getValue(ipar, jpar) * temp[jpar];

    stdv   = (stdv > 0) ? sqrt(stdv) : 0.;
    *stdev = stdv;
  }
}

/****************************************************************************/
/*!
**  Perform the estimation at the grid nodes
**
** \param[in]  seltab     Number of sample definition (0, 1 or 2)
** \param[in]  a          L.H.S. (square) inverted matrix
** \param[in]  zval       Data vector (extended)
** \param[in]  dual       Dual vector
** \param[in]  prior_mean Array of prior means
**
** \param[out] prop1      Working array (Dimension: nlayers)
** \param[out] prop2      Working array (Dimension: nlayers)
** \param[out] covtab     Working array (Dimension = nlayers * nlayers)
** \param[out] b          Working vector (Dimension = neq)
** \param[out] b2         Working vector (Dimension = neq)
** \param[out] baux       Working vector (Dimension = neq)
** \param[out] wgt        Working array (Dimension = neq)
** \param[out] c00        Working array (Dimension = nlayers)
** \param[out] a0         Constant term
** \param[out] cc         Output value
** \param[out] ss         Output value
** \param[out] gs         Output value
** \param[out] post_mean  Array of posterior means
**
*****************************************************************************/
void MLayers::_estimate(VectorInt& seltab,
                        MatrixSquare* a,
                        VectorDouble& zval,
                        VectorDouble& dual,
                        const VectorDouble& prior_mean,
                        VectorDouble& prop1,
                        VectorDouble& prop2,
                        MatrixSquare& covtab,
                        VectorDouble& b,
                        VectorDouble& b2,
                        VectorDouble& baux,
                        VectorDouble& wgt,
                        VectorDouble& c00,
                        MatrixDense& /*fftab*/,
                        MatrixDense& a0,
                        MatrixSquare& cc,
                        MatrixDense& ss,
                        MatrixSquare& gs,
                        VectorDouble& post_mean)
{
  DECLARE_UNUSED(prior_mean);
  double estim;
  VectorDouble coor(2);

  /* Loop on the grid nodes */

  double coefb = 0.;
  double ratio = 0.;
  double stdv  = 0.;
  if (_flagStd && !_flagCumul)
    _getC00(VectorDouble(), covtab, c00);

  for (Id iechout = 0; iechout < _dbout->getNSample(); iechout++)
  {
    OptDbg::setCurrentIndex(iechout + 1);
    if (!_dbout->isActive(iechout)) continue;
    coor[0] = _dbout->getCoordinate(iechout, 0);
    coor[1] = _dbout->getCoordinate(iechout, 1);
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(_dbout, iechout, 1, 0, 0, 0);
    }

    /* Correction in the case of collocation of the bottom surface */

    Id flag_correc = 0;
    if (_colrefB >= 0)
    {
      flag_correc = _getCloseSample(-1, coor);
      if (flag_correc)
      {
        if (_prepareCollocation(iechout, coor, seltab, a, zval, prop1, prop2, covtab, b2,
                                baux, &ratio)) continue;
      }
    }

    /* Loop on the layers */

    for (Id ilayer = 0; ilayer < _nlayers; ilayer++)
    {
      if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH))
        mestitle(2, "Layer #%d", ilayer + 1);

      /* Find the proportions for the target if flag_cumul */

      if (_flagCumul)
      {
        if (_getPropsResult(iechout, ilayer + 1, prop1))
          continue;
      }

      /* Find the C00 terms */

      if (_flagStd && _flagCumul)
        _getC00(prop1, covtab, c00);

      /* Establish the R.H.S. */

      if (_computeRHS(coor, seltab, iechout,
                      ilayer + 1, prop1, prop2, covtab, b)) continue;
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(1, _nech, _neq, _neq, NULL, b.data());

      /* Perform estimation */

      estim = stdv = TEST;
      if (_flagBayes)
        _estimateBayes(c00[ilayer], a, zval, b, wgt, post_mean, a0, cc, ss, gs,
                       &estim, &stdv);
      else
        _estimateRegular(c00[ilayer], a, b, dual, wgt, &estim, &stdv);

      /* Perform the correction (in case of collocated bottom) */

      if (flag_correc)
      {
        double cx = _getCI0(_nlayers, prop1, ilayer + 1, VectorDouble(), covtab);
        if (FFFF(cx)) continue;
        coefb = VH::innerProductVD(b2, b);
        estim += (cx - coefb) * ratio;
      }

      /* Store the result */

      _dbout->setLocVariable(ELoc::Z, iechout, ilayer, estim);
      if (_flagStd) _dbout->setLocVariable(ELoc::Z, iechout, _nlayers + ilayer, stdv);
      if (OptDbg::query(EDbg::RESULTS))
      {
        // TODO a supprimer dans la version finale ... apres DEBUG
        message("Traitement pour debug\n");
        message("flagbayes = %d\n", _flagBayes);
        VH::dump("dual", dual);
        VH::dump("b", b);
        message("Estimate = %lf", ilayer + 1, estim);
        if (_flagStd) message(" - Variance = %lf", stdv * stdv);
        message("\n");
      }
    }
  }
  OptDbg::setCurrentIndex(0);
}

/****************************************************************************/
/*!
 **  Check all the auxiliary variables
 **
 ** \return The number of active samples
 **
 ** \param[in,out]  seltab    Number of sample definition (0, 1 or 2)
 **
 *****************************************************************************/
Id MLayers::_checkAuxiliaryVariables(VectorInt& seltab)
{
  Id newval;
  VectorDouble coor(2);

  Id nechtot = 0;
  for (Id iech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0]   = _dbin->getCoordinate(iech, 0);
    coor[1]   = _dbin->getCoordinate(iech, 1);
    Id ilayer = static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech));
    Id igrid  = _locateSampleInDbout(iech);
    if (IFFFF(igrid)) goto label_suppress;

    // Case of an external drift
    if (_flagExt)
    {
      double drift = _getDriftData(iech, ilayer);
      if (FFFF(drift)) goto label_suppress;
    }

    // Case of a Depth Reference variable
    if (_colrefD >= 0)
    {
      double value = _dbout->getArray(igrid, _colrefD);
      if (FFFF(value)) goto label_suppress;
    }

    // Case of a Bottom Depth Reference variable
    newval = 1;
    if (_colrefB >= 0)
    {
      double value = _dbout->getArray(igrid, _colrefB);
      if (FFFF(value)) goto label_suppress;

      // Check if a duplicate sample must be added:       //
      // - the sample does not belong to the bottom layer //
      // - an analoguous sample does not already exist    //
      if (ilayer < _nlayers && _getCloseSample(iech, coor))
        newval = 2;
    }

    // The sample is valid
    seltab[iech] = newval;
    nechtot += newval;
    continue;

  label_suppress:
    seltab[iech] = 0;
  }

  return (nechtot);
}

/****************************************************************************/
/*!
 **  Convert the results in the Depth scale
 **
 ** \remarks The conversion is performed:
 ** \remarks - if the calculations have been performed in Velocity or Thickness
 ** \remarks - if the calculations have been performed in cumulative or not
 ** \remarks The standard deviation is also transformed
 **
 *****************************************************************************/
void MLayers::_convertResults()
{
  double time  = 0.;
  double stdv  = 0.;
  double depth = 0.;
  double delta = 0.;

  /* If Depth converion is not required, nothing to be done */

  if (!_flagZ) return;

  /* Loop on the target points */

  for (Id iechout = 0; iechout < _dbout->getNSample(); iechout++)
  {

    /* Identify the reference surface */

    double depth0 = (_colrefD >= 0) ? _dbout->getArray(iechout, _colrefD) : 0.;

    /* Identify the reference time (for velocity) */

    double time0 = (_colrefT >= 0) ? _dbout->getArray(iechout, _colrefT) : 0.;

    double depth_prev = depth0;
    double time_prev  = time0;

    /* Loop on the layers */

    for (Id ilayer = 0; ilayer < _nlayers; ilayer++)
    {

      /* Read the estimated value */

      double value = _dbout->getZVariable(iechout, ilayer);
      if (_flagStd) stdv = _dbout->getZVariable(iechout, _nlayers + ilayer);

      if (_flagCumul)
      {

        /* Case when calculations have been processed in cumulative way */

        if (_flagVel)
        {
          time  = _dbout->getFromLocator(_ptime, iechout, ilayer);
          delta = time - time0;
          depth = depth0 + value * delta;
          if (_flagStd) stdv *= delta;
        }
        else
        {
          depth = depth0 + value;
        }
      }
      else
      {

        /* Case when calculations have been processed in individual way */

        if (_flagVel)
        {
          time  = _dbout->getFromLocator(_ptime, iechout, ilayer);
          delta = time - time_prev;
          depth = depth_prev + value * delta;
          if (_flagStd) stdv *= delta;
        }
        else
        {
          depth = depth_prev + value;
        }
        time_prev  = time;
        depth_prev = depth;
      }

      /* Store the transformed results */

      _dbout->setLocVariable(ELoc::Z, iechout, ilayer, depth);
      if (_flagStd) _dbout->setLocVariable(ELoc::Z, iechout, _nlayers + ilayer, stdv);
    }
  }
}

/****************************************************************************/
/*!
 **  Fill the array of drift values at data points
 **
 ** \returns 1 Error return code
 **
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  prop1     Working array (Dimension: nlayers)
 ** \param[in]  drift     Working vector (Dimension: _npar)
 **
 ** \param[out] fftab     Drift Matrix (Dimension: npar * nech)
 **
 *****************************************************************************/
Id MLayers::_fillDriftData(VectorInt& seltab,
                           VectorDouble& prop1,
                           VectorDouble& drift,
                           MatrixDense& fftab)
{
  VectorDouble coor(2);
  fftab.fill(0.);

  for (Id iech = 0, iiech = 0; iech < _dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = _dbin->getCoordinate(iech, 0);
    coor[1] = _dbin->getCoordinate(iech, 1);
    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech)) : _nlayers;

      /* Evaluate the proportion vector */

      if (_getPropsData(iech, ilayer, prop1)) return (1);

      /* Loop on the second sample */

      Id ipos = 0;
      drift.fill(0.);
      for (Id i = 0; i < ilayer; i++)
      {
        double drext = _getDriftData(iech, i + 1);
        if (_fillDrift(coor, prop1[i], drext, &ipos, drift)) return (1);
      }
      for (Id ipar = 0; ipar < _npar; ipar++)
        fftab.setValue(ipar, iech, drift[ipar]);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the posterior mean and variances in Bayesian case
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id MLayers::_calculateDriftBayes(bool verbose,
                                 const VectorDouble& prior_mean,
                                 const VectorDouble& prior_vars,
                                 MatrixSquare* acov,
                                 VectorDouble& zval,
                                 MatrixDense& fftab,
                                 MatrixDense& a0,
                                 MatrixSquare& cc,
                                 MatrixDense& ss,
                                 MatrixSquare& gs,
                                 VectorDouble& post_mean,
                                 MatrixSquare& post_vars) const
{
  Id error = 1;
  MatrixDense ffc(_npar, _nech);
  VectorDouble fm1z(_npar);
  MatrixSquare gg(_npar);
  MatrixSquare invH(_npar);
  MatrixSquare invS(_npar);

  /* Constitute the prior Variance-Covariance matrix */

  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      invS.setValue(ipar, jpar, (ipar == jpar) ? prior_vars[ipar] : 0.);

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Prior Mean", 0, 1, _nlayers, _nbfl, NULL, prior_mean.data());
    print_matrix("Prior Variance", 0, 1, _npar, _npar, NULL, invS.getValues().data());
  }

  /* Invert the Data Variance-Covariance matrix */
  if (matrix_invertFromMatrixSquare(acov, _nech)) goto label_end;

  /* Invert the prior Variance-Covariance matrix */
  if (matrix_invertFromMatrixSquare(&invS, _npar)) goto label_end;

  /* Auxiliary calculations */
  ffc.prodMatMatInPlace(&fftab, acov);
  invH.prodMatMatInPlace(&ffc, &fftab, false, true);
  ffc.prodMatVecInPlace(zval, fm1z);

  /* Calculate the Posterior Variance-Covariance matrix */
  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      post_vars.setValue(ipar, jpar, invS.getValue(ipar, jpar) + invH.getValue(ipar, jpar));
  if (matrix_invertFromMatrixSquare(&post_vars, _npar)) goto label_end;

  /* Calculate the Posterior Mean vector */

  invS.prodMatVecInPlace(prior_mean, post_mean);
  for (Id i = 0; i < _npar; i++)
    fm1z[i] += post_mean[i];
  post_vars.prodMatVecInPlace(fm1z, post_mean);

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Posterior Mean", 0, 1,
                 _nlayers, _nbfl, NULL, post_mean.data());
    print_matrix("Posterior Variance", 0, 1,
                 _npar, _npar, NULL, post_vars.getValues().data());
  }

  /* Modify the Data vector */

  for (Id iech = 0; iech < _nech; iech++)
    for (Id ipar = 0; ipar < _npar; ipar++)
      zval[iech] -= fftab.getValue(ipar, iech) * post_mean[ipar];

  /* Auxiliary arrays prepared for estimation */

  a0.prodMatMatInPlace(&fftab, &post_vars, true);
  ss.prodMatMatInPlace(acov, &fftab, false, true);
  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      invS.setValue(ipar, jpar, post_vars.getValue(ipar, jpar));
  if (matrix_invertFromMatrixSquare(&invS, _npar)) goto label_end;
  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      gg.setValue(ipar, jpar, invH.getValue(ipar, jpar) + invS.getValue(ipar, jpar));
  if (matrix_invertFromMatrixSquare(&gg, _npar)) goto label_end;

  cc.prodNormMatMatInPlace(&ss, &gg, false);
  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      cc.setValue(ipar, jpar, acov->getValue(ipar, jpar) - cc.getValue(ipar, jpar));
  for (Id ipar = 0; ipar < _npar; ipar++)
    for (Id jpar = 0; jpar < _npar; jpar++)
      gs.setValue(ipar, jpar, invH.getValue(ipar, jpar) + invS.getValue(ipar, jpar));
  if (matrix_invertFromMatrixSquare(&gs, _npar)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the sill matrix for a given lag
 **
 ** \return  Error return code (proportions not calculatable)
 **
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  ifirst     Index of the first pair (included)
 ** \param[in]  ilast      Index of the last pair (excluded)
 ** \param[in]  zval       Array containing the sample values
 **
 ** \param[out] nval       Number of relevant pairs
 ** \param[out] distsum    Average distance
 ** \param[out] stat       Working array (Dimension: nlayers * nlayers)
 ** \param[out] phia       Working array for proportions (Dimension: nlayers)
 ** \param[out] phib       Working array for proportions (Dimension: nlayers)
 ** \param[out] atab       Working array (Dimension: nhalf * nhalf)
 ** \param[out] btab       Working array (Dimension: nhalf)
 **
 *****************************************************************************/
Id MLayers::_evaluateLag(Vario_Order* vorder,
                         Id ifirst,
                         Id ilast,
                         VectorDouble& zval,
                         Id* nval,
                         double* distsum,
                         VectorInt& stat,
                         VectorDouble& phia,
                         VectorDouble& phib,
                         MatrixSquare& atab,
                         VectorDouble& btab)
{
  Id iech, jech, iiech, jjech, ilayer, jlayer;
  double dist;

  /* Local initializations */

  (*nval)    = 0;
  (*distsum) = 0.;
  btab.fill(0.);
  atab.fill(0.);
  stat.fill(0.);

  /* Loop on the pairs contributing to the lag */

  for (Id ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iiech, &jjech, &dist);
    vario_order_get_auxiliary(vorder, ipair,
                              reinterpret_cast<char*>(&iech),
                              reinterpret_cast<char*>(&jech));
    double z1 = zval[iiech];
    double z2 = zval[jjech];
    (*distsum) += dist;

    ilayer = static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech));
    if (_getPropsData(iech, ilayer, phia)) return (1);

    jlayer = static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, jech));
    if (_getPropsData(jech, jlayer, phib)) return (1);

    Id ecr1 = 0;
    stat[(ilayer - 1) * _nlayers + (jlayer - 1)] += 1;

    for (Id il1 = 0; il1 < _nlayers; il1++)
      for (Id jl1 = 0; jl1 <= il1; jl1++, ecr1++)
      {
        double fact1 = phia[il1] * phib[jl1];
        if (il1 != jl1) fact1 += phia[jl1] * phib[il1];
        btab[ecr1] += fact1 * z1 * z2;

        Id ecr2 = 0;
        for (Id il2 = 0; il2 < _nlayers; il2++)
          for (Id jl2 = 0; jl2 <= il2; jl2++, ecr2++)
          {
            double fact2 = phia[il2] * phib[jl2];
            if (il2 != jl2) fact2 += phia[jl2] * phib[il2];
            atab.setValue(ecr1, ecr2, atab.getValue(ecr1, ecr2) + fact1 * fact2);
          }
      }

    (*nval)++;
  }
  (*distsum) /= (*nval);
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the variogram for each lag per direction
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose    True for a verbose option
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  zval       Data vector
 ** \param[in]  idir       Rank of the Direction
 **
 ** \param[out] vario      Vario structure
 **
 *****************************************************************************/
Id MLayers::_getVarioCHH(bool verbose,
                         Vario_Order* vorder,
                         VectorDouble& zval,
                         Id idir,
                         Vario* vario)
{
  double distsum;
  Id nval, ifirst, ilast;

  /* Initializations */

  Id error = 1;
  Id nhalf = _nlayers * (_nlayers + 1) / 2;

  /* Core allocation */

  VectorDouble phia(_nlayers);
  VectorDouble phib(_nlayers);
  VectorDouble btab(nhalf);
  MatrixSquare atab(nhalf);
  VectorDouble sill(nhalf);
  VectorInt stat(_nlayers * _nlayers);

  /* Loop on the lags */

  for (Id ilag = 0; ilag < vario->getNLag(idir); ilag++)
  {
    vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);
    Id number = ilast - ifirst;
    if (number <= 0) continue;

    /* Loop on the pairs contributing to the lag */

    if (_evaluateLag(vorder, ifirst, ilast,
                     zval, &nval, &distsum, stat, phia, phib, atab, btab))
      goto label_end;

    if (OptDbg::query(EDbg::VARIOGRAM))
    {
      message("Lag %d\n", ilag + 1);
      print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab.getValues().data());
      print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab.data());
    }

    if (matrix_invertFromMatrixSquare(&atab, nhalf))
    {
      messerr("--> Inversion problem for lag %d", ilag + 1);
      if (verbose)
      {
        /* Matrix must be evaluated (as it has been destroyed by inversion) */
        (void)_evaluateLag(vorder, ifirst,
                           ilast, zval, &nval, &distsum, stat, phia, phib,
                           atab, btab);
        messerr("Number of pairs  = %d", nval);
        messerr("Average distance = %lf", distsum);
        print_imatrix("Number of samples per layer", 0, 1, _nlayers, _nlayers,
                      NULL, stat.data());
        print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab.getValues().data());
        print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab.data());
      }
      continue;
    }
    atab.prodVecMatInPlace(btab, sill);

    /* Optional printout */

    if (OptDbg::query(EDbg::VARIOGRAM))
      print_trimat("C(h)", 2, _nlayers, sill.data());

    /* Store the covariance values */

    Id ijl    = 0;
    Id iadlag = 0;
    for (Id ilayer = 0; ilayer < _nlayers; ilayer++)
      for (Id jlayer = 0; jlayer <= ilayer; jlayer++, ijl++)
      {
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ilag, false, 1);
        vario->setGgByIndex(idir, iadlag, sill[ijl]);
        vario->setHhByIndex(idir, iadlag, distsum);
        vario->setSwByIndex(idir, iadlag, nval);
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ilag, false, -1);
        vario->setGgByIndex(idir, iadlag, sill[ijl]);
        vario->setHhByIndex(idir, iadlag, -distsum);
        vario->setSwByIndex(idir, iadlag, nval);
      }
  }

  /* Set the error return code */

  error = 0;

label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Determine the mean and variance of drift coefficients
 **
 ** \return  Error return code
 **
 ** \param[in]  nech       Number of (active) samples
 ** \param[in]  npar       Number of drift parameters
 ** \param[in]  zval       The data vector (Dimension: neq)
 ** \param[in]  fftab      Drift Matrix (Dimension: npar * nech)
 **
 ** \param[out] mean       Array of means
 ** \param[out] vars       Array of variances
 **
 *****************************************************************************/
Id MLayers::_getPrior(Id nech,
                      Id npar,
                      VectorDouble& zval,
                      MatrixDense& fftab,
                      VectorDouble& mean,
                      MatrixSquare& vars)
{
  MatrixSymmetric atab(npar);
  MatrixSymmetric atab0(npar);
  VectorDouble btab(npar);
  VectorDouble btab0(npar);
  VectorDouble result(npar);

  atab0.fill(0.);
  btab0.fill(0.);
  mean.fill(0.);
  vars.fill(0.);
  result.fill(0.);

  /* Loop on the data */

  for (Id iech = 0; iech < nech; iech++)
    for (Id ipar = 0; ipar < npar; ipar++)
    {
      btab0[ipar] += zval[iech] * fftab.getValue(ipar, iech);
      for (Id jpar = 0; jpar <= ipar; jpar++)
        atab0.updValue(ipar, jpar, EOperator::ADD,
                       fftab.getValue(ipar, iech) * fftab.getValue(jpar, iech));
    }

  /* Bootstrap for the variance-covariance */

  for (Id iech = 0; iech < nech; iech++)
  {
    btab = btab0;
    atab = atab0;

    /* Update the arrays by suppressing the current data */

    for (Id ipar = 0; ipar < npar; ipar++)
    {
      btab[ipar] -= zval[iech] * fftab.getValue(ipar, iech);
      for (Id jpar = 0; jpar <= ipar; jpar++)
        atab.setValue(ipar, jpar, atab.getValue(ipar, jpar) - fftab.getValue(ipar, iech) * fftab.getValue(jpar, iech));
    }

    /* Solve the system */
    if (atab.solve(btab, result)) return 1;

    /* Update the statistics */
    for (Id ipar = 0; ipar < npar; ipar++)
    {
      mean[ipar] += result[ipar];
      for (Id jpar = 0; jpar < npar; jpar++)
        vars.updValue(ipar, jpar, EOperator::ADD, result[ipar] * result[jpar]);
    }
  }

  /* Normalize the results */
  for (Id ipar = 0; ipar < npar; ipar++)
    mean[ipar] /= nech;
  for (Id ipar = 0; ipar < npar; ipar++)
    for (Id jpar = 0; jpar < npar; jpar++)
      vars.setValue(ipar, jpar,
                    vars.getValue(ipar, jpar) / nech - mean[ipar] * mean[jpar]);

  return 0;
}

/**
 * @brief Check if the environment is valid
 *
 * @return true
 * @return false
 */
bool MLayers::checkValid()
{
  if (_dbin->getNDim() != 2)
  {
    messerr("The input Db must be defined in 2-D");
    return false;
  }
  if (_dbout->getNDim() != 2)
  {
    messerr("The output Db must be defined in 2-D");
    return false;
  }
  if (!_dbin->isNVarComparedTo(1)) return false;
  if (!_flagSame && !_dbout->isGrid())
  {
    messerr("If Input and Output are different, Output should be a Grid Db");
    return false;
  }
  if (!_dbin->hasLocator(ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    return false;
  }
  if (_flagExt && _nlayers != _dbout->getNLoc(ELoc::F))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", _nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            _dbout->getNLoc(ELoc::F));
    return false;
  }
  if (_flagVel && _nlayers != get_LOCATOR_NITEM(_dbout, _ptime))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", _nlayers);
    messerr("- the number of time variables in the Output Db File (%d)",
            get_LOCATOR_NITEM(_dbout, _ptime));
    return false;
  }
  _flagBayes   = (!_priorMean.empty() && !_priorVars.empty());
  Id dim_prior = _getNDrift() * _nlayers;
  if (_flagBayes && dim_prior != _getNDrift() * _nlayers)
  {
    messerr("The dimension of the Prior information (%d)", dim_prior);
    messerr("must be equal to %d (nlayers) x %d (nbfl)", _nlayers, _getNDrift());
    return false;
  }
  if (_flagStd && _colrefB >= 0)
  {
    messerr("Calculation of the standard deviation of the estimation error");
    messerr("has not been programmed yet in collocation case");
    return false;
  }
  if (_flagBayes && _colrefB >= 0)
  {
    messerr("Use of Bayesian hypothesis has not been programmed yet");
    messerr("in collocation case");
    return false;
  }
  if (_flagCumul && _colrefB >= 0)
  {
    messerr("Collocation option is not coded when the results are expected");
    messerr("directly expressed in Depth (rather than Thickness)");
    return false;
  }
  return true;
}

VectorInt MLayers::_establishSelection(VectorDouble& props) const
{
  VectorInt seltab(_nechmax);
  for (Id iech = 0; iech < _nechmax; iech++)
  {
    seltab[iech] = 0;
    Id ilayer    = static_cast<Id>(_dbin->getFromLocator(ELoc::LAYER, iech));
    if (ilayer < 1 || ilayer > _nlayers) continue;
    if (_getPropsData(iech, ilayer, props)) continue;
    seltab[iech] = 1;
  }
  return seltab;
}

/****************************************************************************/
/*!
 **  Multi-layers architecture estimation
 **
 ** \return  Pointer to the newly created variables
 **
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
Id MLayers::kriging(bool verbose)
{
  Id nvar = _nlayers;
  if (_flagStd) nvar += _nlayers;
  Id iptr = _dbout->addColumnsByConstant(nvar, TEST, String(), ELoc::Z);

  /* Core allocation */

  VectorDouble prop1(_nlayers);
  VectorDouble prop2(_nlayers);

  /* Calculate the number of active samples */
  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */
  VectorInt seltab = _establishSelection(prop1);
  _nech            = _checkAuxiliaryVariables(seltab);

  /* Complementary allocation */

  VectorDouble b2(_neq);
  VectorDouble b(_neq);
  VectorDouble baux(_neq);
  VectorDouble zval(_neq);
  VectorDouble dual(_neq);
  VectorDouble wgt(_neq);
  VectorDouble c00(_nlayers);
  VectorDouble drift(_npar);
  MatrixSquare atot(_neq);
  MatrixSquare acov(_nech);
  MatrixSquare covtab(_nlayers);

  MatrixDense fftab;
  MatrixDense a0;
  MatrixSquare cc;
  MatrixDense ss;
  MatrixSquare gs;
  MatrixSquare postVars;
  VectorDouble postMean;

  if (_flagBayes)
  {
    fftab.reset(_npar, _nech);
    a0.reset(_nech, _npar);
    cc.reset(_nech, _nech);
    ss.reset(_nech, _npar);
    gs.reset(_npar, _npar);
    postVars.reset(_npar, _npar);
    postMean.resize(_npar);
  }

  /* Establish the kriging matrix */

  _computeLHS(seltab, prop1, prop2, covtab, atot, acov);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    krige_lhs_print(_nech, _neq, _neq, NULL, atot.getValues().data());

  /* Establish the data vector */

  _getDataVector(seltab, zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
  {
    mestitle(0, "Data Vector");
    message("Number of active samples  = %d\n", _nech);
    message("Total number of equations = %d\n", _neq);
    print_matrix("Data", 0, 1, 1, _nech, NULL, zval.data());
  }

  /* Assign the Variance-Covariance matrix */

  MatrixSquare* a = (_flagBayes) ? &acov : &atot;

  /* Calculate the Posterior in the Bayesian case */

  if (_flagBayes)
  {
    if (_fillDriftData(seltab, prop1, drift, fftab)) return -1;
    if (_calculateDriftBayes(verbose, _priorMean, _priorVars, a, zval,
                             fftab, a0, cc, ss, gs, postMean, postVars)) return -1;
  }
  else
  {
    if (matrix_invertFromMatrixSquare(a, _neq, -1)) return -1;
    a->prodMatVecInPlace(zval, dual);
  }

  /* Perform the estimation over the grid nodes */

  _estimate(seltab, a, zval, dual, _priorMean, prop1, prop2, covtab, b, b2,
            baux, wgt, c00, fftab, a0, cc, ss, gs, postMean);

  /* Reconstitute the surfaces (optional) */

  _convertResults();

  return iptr;
}

/****************************************************************************/
/*!
 **  Multi-layers architecture experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  vario      Vario structure
 ** \param[in]  verbose    True for a  verbose option
 **
 *****************************************************************************/
Id MLayers::vario(Vario* vario, bool verbose)
{
  Id error            = 1;
  Vario_Order* vorder = nullptr;

  /* Core allocation */

  VectorDouble prop1(_nlayers);
  VectorDouble zval(_nechmax);

  /* Calculate the number of active samples */
  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */
  VectorInt seltab = _establishSelection(prop1);
  _nech            = _checkAuxiliaryVariables(seltab);

  /* Establish the data vector */

  _getDataVector(seltab, zval);

  /* Subtract the optimal average or drift */

  if (_subtractOptimalDrift(verbose, seltab, zval)) goto label_end;

  /* Evaluate the Geometry */

  vorder = vario_order_manage(1, 1, sizeof(Id), NULL);
  if (vario->computeGeometryMLayers(_dbin, seltab, vorder)) goto label_end;

  /* Evaluate the variogram */

  for (Id idir = 0; idir < vario->getNDir(); idir++)
  {
    if (_getVarioCHH(verbose, vorder, zval, idir, vario)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

label_end:
  vario_order_manage(-1, 1, sizeof(Id), vorder);
  return (error);
}

/****************************************************************************/
/*!
 **  Multi-layers get the mean and prior matrices for Bayesian prior
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id MLayers::calculatePrior()
{
  Id error = 1;

  /* Core allocation */

  VectorDouble props(_nlayers);
  VectorDouble drift(_npar);

  /* Calculate the number of active samples */
  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */
  VectorInt seltab = _establishSelection(props);
  _nech            = _checkAuxiliaryVariables(seltab);

  /* Allocation */

  VectorDouble zval(_neq);
  MatrixDense fftab(_npar, _nech);
  VectorDouble mean(_npar);
  MatrixSquare vars(_npar);

  /* Establish the data vector */

  _getDataVector(seltab, zval);

  /* Establish the drift matrix */

  if (_fillDriftData(seltab, props, drift, fftab)) goto label_end;

  /* Estimate the optimal drift matrices */

  if (_getPrior(_nech, _npar, zval, fftab, mean, vars)) goto label_end;

  /* Print the resulting values */
  message("Number of parameters = %d\n", _npar);
  VH::dump("Means", mean);
  VH::dump("Variances", vars.getValues());

  /* Set the error return code */

  error = 0;

label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Multi-layers architecture estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  flag_same  True if input and output files coincide
 ** \param[in]  flag_Z     True if the output must be converted back into depth
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_cumul True if work is performed in Depth; False in Thickness
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  flag_std   True if the estimation error must be calculated
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  prior_mean Vector of prior means for drift coefficients
 ** \param[in]  prior_vars Vector of prior variances for drift coefficients
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  namerefb   Name of the Bottom Depth variable in Dbout (or -1)
 ** \param[in]  verbose    Verbose option
 ** \param[in]  namconv    Naming convention for output
 **
 *****************************************************************************/
Id multilayers_kriging(Db* dbin,
                       DbGrid* dbout,
                       Model* model,
                       bool flag_same,
                       bool flag_Z,
                       bool flag_vel,
                       bool flag_cumul,
                       bool flag_ext,
                       bool flag_std,
                       bool match_time,
                       Id irf_rank,
                       const VectorDouble& prior_mean,
                       const VectorDouble& prior_vars,
                       const String& namerefd,
                       const String& namereft,
                       const String& namerefb,
                       bool verbose,
                       const NamingConvention& namconv)
{
  MLayers lmlayers = MLayers();
  lmlayers.setGeneral(dbin, dbout, model,
                      flag_std, flag_same, flag_vel, flag_cumul, flag_ext, flag_Z,
                      match_time, irf_rank, prior_mean, prior_vars, namerefd, namereft, namerefb);

  if (!lmlayers.checkValid()) return 1;

  Id iptr = lmlayers.kriging(verbose);

  Id nvar = dbout->getNLoc(ELoc::Z);

  namconv.setNamesAndLocators(dbout, iptr, "estim", nvar);

  return 0;
}

/****************************************************************************/
/*!
 **  Multi-layers architecture experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  vario      Vario structure
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  verbose    True for a  verbose option
 **
 *****************************************************************************/
Id multilayers_vario(Db* dbin,
                     DbGrid* dbout,
                     Vario* vario,
                     bool flag_vel,
                     bool flag_ext,
                     Id irf_rank,
                     bool match_time,
                     const String& namerefd,
                     const String& namereft,
                     bool verbose)
{
  MLayers lmlayers = MLayers();
  lmlayers.setGeneral(dbin, dbout, nullptr,
                      false, false, flag_vel, false, flag_ext, false,
                      match_time, irf_rank, VectorDouble(), VectorDouble(), namerefd, namereft, String());

  if (!lmlayers.checkValid()) return 1;

  return lmlayers.vario(vario, verbose);
}

/****************************************************************************/
/*!
 **  Multi-layers get the mean and prior matrices for Bayesian prior
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  flag_same  True if input and output files coincide
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  namerefb   Name of the Bottom Depth variable in Dbout (or -1)
 **
 *****************************************************************************/
Id multilayers_getPrior(Db* dbin,
                        DbGrid* dbout,
                        Model* model,
                        bool flag_same,
                        bool flag_vel,
                        bool flag_ext,
                        Id irf_rank,
                        bool match_time,
                        const String& namerefd,
                        const String& namereft,
                        const String& namerefb)
{
  MLayers lmlayers = MLayers();
  lmlayers.setGeneral(dbin, dbout, model,
                      false, flag_same, flag_vel, false, flag_ext, false,
                      match_time, irf_rank, VectorDouble(), VectorDouble(), namerefd, namereft, namerefb);

  if (!lmlayers.checkValid()) return 1;

  return lmlayers.calculatePrior();
}

} // namespace gstlrn
