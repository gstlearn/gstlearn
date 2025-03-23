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
#include "Db/Db.hpp"
#include "Covariances/ANoStat.hpp"
#include "Covariances/ParamId.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include <memory>
#include <unordered_map>

typedef std::unordered_map<ParamId, std::shared_ptr<ANoStat>, ParamIdHash, ParamIdEqual> mapNoStat;

class GSTLEARN_EXPORT TabNoStat: public AStringable, public ICloneable
{
public:
  TabNoStat();
  TabNoStat(const TabNoStat& m);
  TabNoStat& operator=(const TabNoStat& m);
  virtual ~TabNoStat();

  IMPLEMENT_CLONING(TabNoStat)
  bool isNoStat() const { return !_items.empty(); }
  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;
  void informMeshByMesh(const AMesh* amesh, const EConsElem& econs) const;
  void informMeshByApex(const AMesh* amesh, const EConsElem& econs) const;
  void informDbIn(const Db* dbin, const EConsElem& econs) const;
  void informDbOut(const Db* dbout, const EConsElem& econs) const;
  int size() const { return _items.size(); }
  bool empty() const { return _items.empty(); }
  void updateDescription();
  const mapNoStat& getTable() const { return _items; }
  bool isValid(const EConsElem& econs) const;
  virtual int addElem(std::shared_ptr<ANoStat>& nostat, const EConsElem& econs, int iv1 = 0, int iv2 = 0);
  virtual int removeElem(const EConsElem& econs, int iv1 = 0, int iv2 = 0);
  void clear();
  void setDbNoStatRef(const Db* dbref);
  void setDbNoStatRef(std::shared_ptr<const Db>& dbref);
  std::shared_ptr<const Db> getDbNoStatRef() const;
  const Db* getDbNoStatRefRaw() const;
  void informCoords(const VectorVectorDouble& coords,
                    const EConsElem& econs,
                    int iv1,
                    int iv2,
                    VectorDouble& result) const;
  String toString(const AStringFormat* strfmt = nullptr) const override;
  String toStringInside(const AStringFormat* strfmt = nullptr, int i = 0) const;
  bool isElemDefined(const EConsElem& econs, int iv1 = 0, int iv2 = 0) const;
  std::shared_ptr<ANoStat> getElem(const EConsElem& econs, int iv1 = 0, int iv2 = 0);

protected:

private:
  virtual void _clear() {};
  virtual void _updateDescription() {};
  virtual bool _isValid(const EConsElem& econs) const
  {
    DECLARE_UNUSED(econs)
    return false;
  };

private:
  mapNoStat _items;
  std::shared_ptr<const Db> _dbNoStatRef;
};
