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
#include "Covariances/CovCalcMode.hpp"

#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

CovCalcMode::CovCalcMode(const ECalcMember& member,
                         bool asVario,
                         bool normalized,
                         bool filterNugget,
                         unsigned int keepOnlyCovIdx,
                         bool unitary,
                         int envelop,
                         int orderVario)
: _member(member),
  _asVario(asVario),
  _normalized(normalized),
  _filterNugget(filterNugget),
  _keepOnlyCovIdx(keepOnlyCovIdx),
  _unitary(unitary),
  _envelop(envelop),
  _orderVario(orderVario)
{
}

CovCalcMode::CovCalcMode(const CovCalcMode &r)
: _member(r._member),
  _asVario(r._asVario),
  _normalized(r._normalized),
  _filterNugget(r._filterNugget),
  _keepOnlyCovIdx(r._keepOnlyCovIdx),
  _unitary(r._unitary),
  _envelop(r._envelop),
  _orderVario(r._orderVario)
{

}

CovCalcMode& CovCalcMode::operator=(const CovCalcMode &r)
{
  if (this != &r)
  {
    _member = r._member;
    _asVario = r._asVario;
    _normalized = r._normalized;
    _filterNugget = r._filterNugget;
    _keepOnlyCovIdx = r._keepOnlyCovIdx;
    _unitary = r._unitary;
    _envelop = r._envelop;
    _orderVario = r._orderVario;
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
          _orderVario == r._orderVario);
}

/**
 * Update the CovCalcMode structure according to the input arguments
 * @param nugget_opt  Option for nugget effect basic structure
 ** \li                      (Default: 0)
 ** \li                       0 : no particular option
 ** \li                       1 : discard the nugget effect
 ** \li                      -1 : only consider the nugget effect
 * @param nostd       0 standard; +-1 special; ITEST normalized
 **                          (Default: 0)
 * @param member      Member of the Kriging System (ECalcMember)
 **                          (Default: ECalcMember::LHS)
 * @param icov_r      rank of the target covariance or -1 for all
 **                          (Default: -1)
 * @param flag_norm   1 if the model is normalized
 **                          (Default: 0)
 * @param flag_cov    1 if the result must be given in covariance
 **                          (Default: 1)
 */
void CovCalcMode::update(int nugget_opt,
                         int nostd,
                         const ECalcMember& member,
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
