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
#pragma once

#include "Covariances/TabNoStat.hpp"
#include "Enum/EConsElem.hpp"

class NoStatElem;

namespace gstlrn
{
 
class GSTLEARN_EXPORT TabNoStatCovAniso : public TabNoStat
{
  IMPLEMENT_CLONING(TabNoStatCovAniso)
  public:
  TabNoStatCovAniso();
  TabNoStatCovAniso(const TabNoStatCovAniso &m);
  TabNoStatCovAniso& operator= (const TabNoStatCovAniso &m);
  virtual ~TabNoStatCovAniso();
  Id getNAngles() const {return _nAngles;}
  Id getNRanges() const {return _nRanges;}
  Id getNScales() const {return _nScales;}
  
  bool isParam()   const {return _param;}
  bool isDefinedForTensor()  const {return _definedForTensor;}
  bool isDefinedForAnisotropy() const;
  bool isDefinedForRotation() const;
  Id addElem(std::shared_ptr<ANoStat> &nostat, const EConsElem &econs, Id iv1=0, Id iv2 = 0) override;
  Id removeElem(const EConsElem &econs, Id iv1=0, Id iv2 = 0) override;

  private:
    void _clear() override;
    void _updateDescription() override;
    bool _isValid(const EConsElem &econs) const override;

  private:
    Id  _nAngles;
    Id  _nRanges;
    Id  _nScales;
    Id  _nTensor;
    bool _param;
    bool _definedForAnisotropy;
    bool _definedByAnglesAndScales;
    bool _definedForRotation;
    bool _definedForTensor;

};
}