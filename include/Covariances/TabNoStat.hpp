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

#include "Basic/VectorNumT.hpp"
#include "Covariances/ANoStat.hpp"
#include "Covariances/ParamId.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/EConsElem.hpp"
#include <memory>
#include <unordered_map>

typedef std::unordered_map<ParamId,std::shared_ptr<ANoStat>,ParamIdHash,ParamIdEqual> mapNoStat;

class GSTLEARN_EXPORT TabNoStat 
{
  public:
  TabNoStat();
  TabNoStat(const TabNoStat &m);
  TabNoStat& operator= (const TabNoStat &m);
  bool isNoStat() const { return !_items.empty();}
  void informMeshByMesh(const AMesh* amesh);
  void informMeshByApex(const AMesh* amesh);
  void informDbIn(const Db* dbin);
  void informDbOut(const Db* dbout);
  void informMeshByMesh(const AMesh* amesh, const EConsElem & econs);
  void informMeshByApex(const AMesh* amesh, const EConsElem & econs);
  void informDbIn(const Db* dbin, const EConsElem & econs);
  void informDbOut(const Db* dbout, const EConsElem & econs);
  bool isDefinedForVariance() const {return _definedForVariance;}
  int getNSills()  const {return _nSills;}
  void updateDescription();
  const mapNoStat&  getTable() const { return _items;}
  bool isValid(const EConsElem &econs) const;
  virtual int addElem(std::shared_ptr<ANoStat> &nostat, const EConsElem &econs, int iv1=0, int iv2 = 0);
  virtual int removeElem(const EConsElem &econs, int iv1=0, int iv2 = 0);
  virtual ~TabNoStat();
  void setDbNoStatRef(const DbGrid* dbref){ _dbNoStatRef = dbref;}
  const DbGrid* getDbNoStatRef() const {return _dbNoStatRef;}
  void informCoords(const VectorVectorDouble &coords,
                    const EConsElem &econs,
                    int iv1, 
                    int iv2, 
                    VectorDouble& result) const;

  bool isElemDefined(const EConsElem &econs, int iv1 = 0, int iv2 = 0) const;
  std::shared_ptr<ANoStat> getElem(const EConsElem &econs, int iv1 = 0, int iv2 =0);
protected:
private:
  virtual void _updateDescription() {};
  virtual bool _isValid(const EConsElem &econs) const;
private :
  mapNoStat _items;
  const DbGrid* _dbNoStatRef;
  bool _definedForVariance;
  int  _nSills;
};
