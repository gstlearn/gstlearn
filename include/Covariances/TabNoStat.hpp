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

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ANoStat.hpp"
#include "Covariances/ParamId.hpp"
#include "Enum/EConsElem.hpp"
#include <memory>
#include <unordered_map>

class Db;
typedef std::unordered_map<ParamId,std::shared_ptr<ANoStat>,ParamIdHash,ParamIdEqual> mapNoStat;

class GSTLEARN_EXPORT TabNoStat : public AStringable, public ICloneable
{
  IMPLEMENT_CLONING(TabNoStat)
  public:
  TabNoStat();
  TabNoStat(const TabNoStat &m);
  TabNoStat& operator= (const TabNoStat &m);
  bool isNoStat() const { return !_items.empty();}
  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;
  void informMeshByMesh(const AMesh* amesh, const EConsElem & econs) const;
  void informMeshByApex(const AMesh* amesh, const EConsElem & econs) const;
  void informDbIn(const Db* dbin, const EConsElem & econs) const;
  void informDbOut(const Db* dbout, const EConsElem & econs) const;
  bool isDefinedForVariance() const {return _definedForVariance;}
  int getNSills()  const {return _nSills;}
  void updateDescription();
  const mapNoStat&  getTable() const { return _items;}
  bool isValid(const EConsElem &econs) const;
  virtual int addElem(std::shared_ptr<ANoStat> &nostat, const EConsElem &econs, int iv1=0, int iv2 = 0);
  virtual int removeElem(const EConsElem &econs, int iv1=0, int iv2 = 0);
  virtual ~TabNoStat();
  void clear();
  void setDbNoStatRef(const Db* dbref){ _dbNoStatRef = dbref;}
  const Db* getDbNoStatRef() const {return _dbNoStatRef;}
  void informCoords(const VectorVectorDouble &coords,
                    const EConsElem &econs,
                    int iv1, 
                    int iv2, 
                    VectorDouble& result) const;
  String toString(const AStringFormat* strfmt = nullptr) const override;
  String toStringInside(const AStringFormat* strfmt = nullptr,int i = 0) const;
  bool isElemDefined(const EConsElem &econs, int iv1 = 0, int iv2 = 0) const;
  std::shared_ptr<ANoStat> getElem(const EConsElem &econs, int iv1 = 0, int iv2 =0);
protected:
private:
  virtual void _clear(){};
  virtual void _updateDescription() {};
  virtual bool _isValid(const EConsElem &econs) const;
private :
  mapNoStat _items;
  const Db* _dbNoStatRef;
  bool _definedForVariance;
  int  _nSills;
};
