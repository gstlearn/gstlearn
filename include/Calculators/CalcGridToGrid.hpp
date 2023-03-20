/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "ACalcDbToDb.hpp"

#include "Db/DbGrid.hpp"
#include "Basic/NamingConvention.hpp"

class GSTLEARN_EXPORT CalcGridToGrid: public ACalcDbToDb
{
public:
  CalcGridToGrid();
  CalcGridToGrid(const CalcGridToGrid &r) = delete;
  CalcGridToGrid& operator=(const CalcGridToGrid &r) = delete;
  virtual ~CalcGridToGrid();

  void setFlagCopy(bool flagCopy)     { _flagCopy = flagCopy; }
  void setFlagExpand(bool flagExpand) { _flagExpand = flagExpand; }
  void setFlagShrink(bool flagShrink) { _flagShrink = flagShrink; }
  void setFlagInter(bool flagInter)   { _flagInter = flagInter;   }
  void setNameBots(const VectorString name_bots) { _nameBots = name_bots; }
  void setNameTops(const VectorString name_tops) { _nameTops = name_tops; }

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  virtual int  _getNVar() const;

private:
  bool _g2gCopy();
  bool _g2gExpand();
  bool _g2gShrink();
  bool _g2gInter();
  int _compareInMinusOut() const;
  void _reduceIndices(const VectorInt& indgIn, VectorInt& indgOut);
  bool _loadExtrema(int nvar, int iech, const VectorInt& iuids, VectorDouble& coor);
  double _interpolate(int nvar,
                      double valTop,
                      double valBot,
                      const VectorDouble &coorTop,
                      const VectorDouble &coorBot,
                      const VectorDouble &coorOut);

private:
  int  _iattOut;
  bool _flagCopy;
  bool _flagExpand;
  bool _flagShrink;
  int  _iattAux;
  bool _flagInter;
  VectorString _nameTops;
  VectorString _nameBots;
};

GSTLEARN_EXPORT int dbg2gCopy(DbGrid *dbin,
                              DbGrid *dbout,
                              const NamingConvention &namconv = NamingConvention(
                                  "Copy"));
GSTLEARN_EXPORT int dbg2gExpand(DbGrid *dbin,
                                DbGrid *dbout,
                                const NamingConvention &namconv = NamingConvention(
                                    "Expand"));
GSTLEARN_EXPORT int dbg2gShrink(DbGrid *dbin,
                                DbGrid *dbout,
                                const NamingConvention &namconv = NamingConvention(
                                    "Shrink"));
GSTLEARN_EXPORT int dbg2gInterpolate(DbGrid *dbin,
                                     DbGrid *dbout,
                                     const VectorString &tops = VectorString(),
                                     const VectorString &bots = VectorString(),
                                     const NamingConvention &namconv = NamingConvention(
                                         "Interpolation", false));
