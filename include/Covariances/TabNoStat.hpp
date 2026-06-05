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

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ANoStat.hpp"
#include "Covariances/ParamId.hpp"
#include "Db/Db.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include <memory>

#ifdef USE_BOOST_SPAN
// Hash containers use user-specified hash operators (e.g., ParamIdHash) that
// return a 64 bits hash. This is internally folded into a bucket index which
// depends on the current container size. With VS2017 (and maybe other
// versions), this is done internally by std::_Hash::_Hashval() in xhash by
// masking (&) with the size. This can cause some hash collisions, because in
// effect this amounts to only looking at the lowest-order bits of the full
// hash (those that match the size).
//
// boost::unordered_map implements a smarter folding strategy that uses all
// bits, to avoid this.
//
// gcc does not seem to have this problem, even though the folding of the hash
// is apparently also done by default with a modulo (which is the same), see
// the _Mod_range_hashing struct in hashtable_policy.h?
#include <boost/unordered_map.hpp>
#else
#include <unordered_map>
#endif

namespace gstlrn
{
#ifdef USE_BOOST_SPAN
  typedef boost::
    unordered_map<ParamId, std::shared_ptr<ANoStat>, ParamIdHash, ParamIdEqual>
      mapNoStat;
#else
  typedef std::
    unordered_map<ParamId, std::shared_ptr<ANoStat>, ParamIdHash, ParamIdEqual>
      mapNoStat;
#endif

  class GSTLEARN_EXPORT TabNoStat: public AStringable,
                                   public ICloneable,
                                   public ASerializable
  {
  public:
    TabNoStat();
    TabNoStat(const TabNoStat& m);
    TabNoStat& operator=(const TabNoStat& m);
    virtual ~TabNoStat();

    /// ICloneable Interface
    IMPLEMENT_CLONING(TabNoStat)

    /// ASerializable Interface
    String getNFName() const override { return "TabNoStat"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    bool isNoStat() const { return !_items.empty(); }

    void informMeshByMesh(const AMesh* amesh) const;
    void informMeshByApex(const AMesh* amesh) const;
    void informDbIn(const Db* dbin) const;
    void informDbOut(const Db* dbout) const;
    void informMeshByMesh(const AMesh* amesh, const EConsElem& econs) const;
    void informMeshByApex(const AMesh* amesh, const EConsElem& econs) const;
    void informDbIn(const Db* dbin, const EConsElem& econs) const;
    void informDbOut(const Db* dbout, const EConsElem& econs) const;

    Id size() const { return static_cast<Id>(_items.size()); }

    bool empty() const { return _items.empty(); }

    void updateDescription();

    const mapNoStat& getTable() const { return _items; }

    bool isValid(const EConsElem& econs) const;
    virtual Id addElem(
      std::shared_ptr<ANoStat>& nostat,
      const EConsElem& econs,
      Id iv1 = 0,
      Id iv2 = 0);
    virtual Id removeElem(const EConsElem& econs, Id iv1 = 0, Id iv2 = 0);
    void clear();
    void setDbNoStatRef(const Db* dbref);
    void setDbNoStatRef(std::shared_ptr<const Db>& dbref);
    std::shared_ptr<const Db> getDbNoStatRef() const;
    const Db* getDbNoStatRefRaw() const;
    void informCoords(
      const VectorVectorDouble& coords,
      const EConsElem& econs,
      Id iv1,
      Id iv2,
      VectorDouble& result) const;
    String toString(const AStringFormat* strfmt = nullptr) const override;
    bool isElemDefined(const EConsElem& econs, Id iv1 = 0, Id iv2 = 0) const;
    std::shared_ptr<ANoStat>
      getElem(const EConsElem& econs, Id iv1 = 0, Id iv2 = 0);
    bool variableExistsInDb(const String& namecol) const;
    String
      toStringInside(const AStringFormat* strfmt = nullptr, Id i = 0) const;

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
} // namespace gstlrn
