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

#include "gstlearn_export.hpp"

#include "ACalcDbToDb.hpp"

#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT CalcGridToGrid: public ACalcDbToDb
{
public:
  CalcGridToGrid();
  CalcGridToGrid(const CalcGridToGrid& r)            = delete;
  CalcGridToGrid& operator=(const CalcGridToGrid& r) = delete;
  virtual ~CalcGridToGrid();

  void setFlagCopy(bool flagCopy) { _flagCopy = flagCopy; }
  void setFlagExpand(bool flagExpand) { _flagExpand = flagExpand; }
  void setFlagShrink(bool flagShrink) { _flagShrink = flagShrink; }
  void setFlagInter(bool flagInter) { _flagInter = flagInter; }
  void setNameBots(const VectorString& name_bots) { _nameBots = name_bots; }
  void setNameTops(const VectorString& name_tops) { _nameTops = name_tops; }

protected:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

private:
  bool _g2gCopy();
  bool _g2gExpand();
  bool _g2gShrink();
  bool _g2gInter();
  Id _compareInMinusOut() const;
  static void _reduceIndices(const VectorInt& indgIn, VectorInt& indgOut);
  bool _loadExtrema(Id nvar, Id iech, const VectorInt& iuids, VectorDouble& coor);
  double _interpolate(Id nvar,
                      double valTop,
                      double valBot,
                      const VectorDouble& coorTop,
                      const VectorDouble& coorBot,
                      const VectorDouble& coorOut);

private:
  Id _iattOut;
  bool _flagCopy;
  bool _flagExpand;
  bool _flagShrink;
  Id _iattAux;
  bool _flagInter;
  VectorString _nameTops;
  VectorString _nameBots;
};

GSTLEARN_EXPORT Id dbg2gCopy(DbGrid* dbin,
                              DbGrid* dbout,
                              const NamingConvention& namconv = NamingConvention(
                                "Copy"));
GSTLEARN_EXPORT Id dbg2gExpand(DbGrid* dbin,
                                DbGrid* dbout,
                                const NamingConvention& namconv = NamingConvention(
                                  "Expand"));
GSTLEARN_EXPORT Id dbg2gShrink(DbGrid* dbin,
                                DbGrid* dbout,
                                const NamingConvention& namconv = NamingConvention(
                                  "Shrink"));
GSTLEARN_EXPORT Id dbg2gInterpolate(DbGrid* dbin,
                                     DbGrid* dbout,
                                     const VectorString& tops        = VectorString(),
                                     const VectorString& bots        = VectorString(),
                                     const NamingConvention& namconv = NamingConvention(
                                       "Interpolation",
                                       false));
} // namespace gstlrn
