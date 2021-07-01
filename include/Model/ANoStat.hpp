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
#pragma once

#include "Model/ConsItem.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_enum.h"

class Model;

class ANoStat : public AStringable
{
public:
  ANoStat();
  ANoStat(const VectorString& codes);
  ANoStat(const ANoStat &m);
  ANoStat& operator= (const ANoStat &m);
  virtual ~ANoStat();

  virtual String toString(int level = 0) const override;

  bool isDefined() const { return ! _items.empty(); }
  bool isDefinedByCov(int igrf, int icov);
  bool isDefinedByType(int igrf, ENUM_CONS type);
  bool isDefinedByCovType(int igrf, int icov, ENUM_CONS type);
  bool isDefined(int igrf, int icov, ENUM_CONS type, int iv1=0, int iv2=0);
  bool isDefinedforAnisotropy(int igrf, int icov);

  virtual double getValue(int igrf,
                          int icov,
                          ENUM_CONS type,
                          int iv1,
                          int iv2,
                          int icas,
                          int rank) const = 0;
  virtual double getValue(int ipar, int icas, int iech) const = 0;

  void addNoStatElem(int igrf, int icov, ENUM_CONS type, int iv1, int iv2);
  void addNoStatElems(const VectorString& codes);
  void addNoStatElem(const ConsItem& item);

  int getRank   (int igrf, int icov, ENUM_CONS type, int iv1, int iv2) const;
  int getIGrf(int ipar) const { return _items[ipar].getIGrf(); }
  int getICov(int ipar) const { return _items[ipar].getICov(); }
  ENUM_CONS getType(int ipar) const { return _items[ipar].getType(); }
  int getIV1 (int ipar) const { return _items[ipar].getIV1(); }
  int getIV2 (int ipar) const { return _items[ipar].getIV2(); }
  int getNoStatElemNumber() const { return _items.size(); }
  const std::vector<ConsItem>& getNoStat() const { return _items; }
  const ConsItem getNoStat(int ipar) const { return _items[ipar]; }

  const int attachModel(const Model* model);

  bool matchIGrf(int ipar, int igrf0) const { return _items[ipar].matchIGrf(igrf0); }
  bool matchICov(int ipar, int icov0) const { return _items[ipar].matchICov(icov0); }
  bool matchType(int ipar, ENUM_CONS type0) const { return _items[ipar].matchType(type0); }
  bool matchIV1(int ipar, int iv10) const { return _items[ipar].matchIV1(iv10); }
  bool matchIV2(int ipar, int iv20) const { return _items[ipar].matchIV2(iv20); }

  const std::vector<ConsItem>& getItems() const { return _items; }
  const ConsItem getItems(int ipar) const { return _items[ipar]; }

private:
  int _understandCode(const String& code,
                      int *igrf,
                      int *icov,
                      ENUM_CONS *type,
                      int *iv1,
                      int *iv2);
  void _updateFromModel(const Model* model);

private:
  std::vector<ConsItem> _items;
};
