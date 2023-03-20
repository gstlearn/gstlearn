/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovCalcMode.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"

CovCalcMode::CovCalcMode(const ECalcMember& member,
                         bool asVario,
                         bool normalized,
                         bool filterNugget,
                         unsigned int keepOnlyCovIdx,
                         bool unitary,
                         int envelop,
                         int orderVario)
: AStringable(),
  _member(member),
  _asVario(asVario),
  _normalized(normalized),
  _filterNugget(filterNugget),
  _keepOnlyCovIdx(keepOnlyCovIdx),
  _unitary(unitary),
  _envelop(envelop),
  _orderVario(orderVario),
  _indexClass(0),
  _covFiltered()
{
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
    : AStringable(r),
      _member(r._member),
      _asVario(r._asVario),
      _normalized(r._normalized),
      _filterNugget(r._filterNugget),
      _keepOnlyCovIdx(r._keepOnlyCovIdx),
      _unitary(r._unitary),
      _envelop(r._envelop),
      _orderVario(r._orderVario),
      _indexClass(r._indexClass),
      _covFiltered(r._covFiltered)
{

}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _member = r._member;
    _asVario = r._asVario;
    _normalized = r._normalized;
    _filterNugget = r._filterNugget;
    _keepOnlyCovIdx = r._keepOnlyCovIdx;
    _unitary = r._unitary;
    _envelop = r._envelop;
    _orderVario = r._orderVario;
    _indexClass = r._indexClass;
    _covFiltered = r._covFiltered;
  }
  return *this;
}
CovCalcMode::~CovCalcMode()
{
}

bool CovCalcMode::isEqual(const CovCalcMode &r) const
{
  return (_member == r._member &&
          _asVario == r._asVario &&
          _normalized == r._normalized &&
          _filterNugget == r._filterNugget &&
          _keepOnlyCovIdx == r._keepOnlyCovIdx &&
          _unitary == r._unitary &&
          _envelop == r._envelop &&
          _orderVario == r._orderVario &&
          _indexClass == r._indexClass &&
          _covFiltered == r._covFiltered);
}

/**
 * Update the CovCalcMode structure according to the input arguments
 * @param member      Member of the Kriging System (ECalcMember)
 **                          (Default: ECalcMember::LHS)
 * @param nugget_opt  Option for nugget effect basic structure
 ** \li                      (Default: 0)
 ** \li                       0 : no particular option
 ** \li                       1 : discard the nugget effect
 ** \li                      -1 : only consider the nugget effect
 * @param nostd       0 standard; +-1 special; ITEST normalized
 **                          (Default: 0)
 * @param icov_r      rank of the target covariance or -1 for all
 **                          (Default: -1)
 * @param flag_norm   1 if the model is normalized
 **                          (Default: 0)
 * @param flag_cov    1 if the result must be given in covariance
 **                          (Default: 1)
 */
void CovCalcMode::update(const ECalcMember& member,
                         int nugget_opt,
                         int nostd,
                         int icov_r,
                         int flag_norm,
                         int flag_cov)
{
  if (nugget_opt == -1)
    my_throw("nugget_opt == -1 not yet implemented");
  else
    _filterNugget = (nugget_opt == 1);
  if (IFFFF(nostd))
    _unitary = true;
  else
  {
    if (nostd != 0) _envelop = nostd;
  }
  _member = member;
  _keepOnlyCovIdx = icov_r;
  _asVario = flag_cov == 0;
  _normalized = flag_norm == 1;
}

bool CovCalcMode::getCovFiltered(int i) const
{
  if (_covFiltered.empty()) return false;
  if (i < 0 || i >= (int) _covFiltered.size()) return false;
  return _covFiltered[i];
}

void CovCalcMode::setCovFiltered(int i, bool status)
{
  if (_covFiltered.empty()) return;
  if (i < 0 || i >= (int) _covFiltered.size()) return;
  _covFiltered[i] = status;
}

void CovCalcMode::setAllCovFiltered(int ncov, bool status)
{
  if (_covFiltered.empty())
    _covFiltered.resize(ncov, status);
  else
  {
    for (int i=0; i<(int) _covFiltered.size(); i++)
      _covFiltered[i] = status;
  }
}
