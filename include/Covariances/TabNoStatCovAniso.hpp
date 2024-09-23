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

class GSTLEARN_EXPORT TabNoStatCovAniso : public TabNoStat
{
  public:
  TabNoStatCovAniso();
  TabNoStatCovAniso(const TabNoStatCovAniso &m);
  TabNoStatCovAniso& operator= (const TabNoStatCovAniso &m);
  virtual ~TabNoStatCovAniso();
  int getNAngles() const {return _nAngles;}
  int getNRanges() const {return _nRanges;}
  int getNScales() const {return _nScales;}
  
  bool isParam()   const {return _param;}
  bool isDefinedForTensor()  const {return _definedForTensor;}
  bool isDefinedForAnisotropy() const;
  bool isDefinedForRotation() const;
  int addElem(std::shared_ptr<ANoStat> &nostat, const EConsElem &econs, int iv1=0, int iv2 = 0) override;
  int removeElem(const EConsElem &econs, int iv1=0, int iv2 = 0) override;

  private:
    void _updateDescription() override;
    bool _isValid(const EConsElem &econs) const override;

  private:
    int  _nAngles;
    int  _nRanges;
    int  _nScales;
    int  _nTensor;
    bool _param;
    bool _definedForAnisotropy;
    bool _definedByAnglesAndScales;
    bool _definedForRotation;
    bool _definedForTensor;

};
