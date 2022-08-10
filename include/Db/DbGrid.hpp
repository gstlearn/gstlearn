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

#include "gstlearn_export.hpp"
#include "geoslib_d.h"

#include "Db/ELoadBy.hpp"
#include "Db/PtrGeos.hpp"
#include "Db/Db.hpp"

#include "Basic/Grid.hpp"
#include "Basic/Limits.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/CSVformat.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"

class Polygons;

/**
 * Class containing a Data Set organized as a regular Grid
 */
class GSTLEARN_EXPORT DbGrid: public Db
{
public:
  DbGrid();
  DbGrid(const DbGrid& r);
  DbGrid& operator=(const DbGrid& r);
  virtual ~DbGrid();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbGrid)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  inline bool isGrid() const override { return true; }
  double getCoordinate(int iech, int idim, bool flag_rotate=true) const override;
  double getUnit(int idim = 0) const override;
  int getNDim() const override;
  bool mayChangeSampleNumber() const override { return false; }

  static DbGrid* createFromNF(const String& neutralFilename,
                              bool verbose = true);
  int reset(const VectorInt& nx,
            const VectorDouble& dx = VectorDouble(),
            const VectorDouble& x0 = VectorDouble(),
            const VectorDouble& angles = VectorDouble(),
            const ELoadBy& order = ELoadBy::SAMPLE,
            const VectorDouble& tab = VectorDouble(),
            const VectorString& names = VectorString(),
            const VectorString& locatorNames = VectorString(),
            int flag_add_rank = 1);
  int resetCoveringDb(Db* db,
                      const VectorInt& nodes,
                      const VectorDouble& dcell = VectorDouble(),
                      const VectorDouble& origin = VectorDouble(),
                      const VectorDouble& margin = VectorDouble());
  int resetFromPolygon(Polygons* polygon,
                       const VectorInt& nodes,
                       const VectorDouble& dcell,
                       int flag_add_rank);
  static DbGrid* create(const VectorInt& nx,
                        const VectorDouble& dx = VectorDouble(),
                        const VectorDouble& x0 = VectorDouble(),
                        const VectorDouble& angles = VectorDouble(),
                        const ELoadBy& order = ELoadBy::SAMPLE,
                        const VectorDouble& tab = VectorDouble(),
                        const VectorString& names = VectorString(),
                        const VectorString& locatorNames = VectorString(),
                        int flag_add_rank = 1);
  static DbGrid* createCoveringDb(Db* dbin,
                                  const VectorInt& nodes = VectorInt(),
                                  const VectorDouble& dcell = VectorDouble(),
                                  const VectorDouble& origin = VectorDouble(),
                                  const VectorDouble& margin = VectorDouble());
  static DbGrid* createFromPolygon(Polygons* polygon,
                                   const VectorInt& nodes,
                                   const VectorDouble& dcell,
                                   int flag_add_rank = 1);
  static DbGrid* createCoarse(DbGrid *dbin,
                              const VectorInt &nmult,
                              int flag_add_rank);
  static DbGrid* createRefine(DbGrid *dbin,
                              const VectorInt &nmult,
                              int flag_add_rank);
  static bool migrateAllVariables(Db *dbin, Db *dbout, int flag_add_rank);

  inline const Grid& getGrid() const { return _grid; }
  void generateCoordinates(const String& radix = "x");

  VectorDouble getColumnSubGrid(const String& name,
                                int idim0,
                                int rank,
                                bool useSel = false);

  int gridDefine(const VectorInt& nx,
                 const VectorDouble& dx = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 const VectorDouble& angles = VectorDouble());
  void gridCopyParams(int mode, const Grid& gridaux);
  bool isSameGrid(const Grid& grid) const;
  bool isSameGridMesh(const DbGrid& dbaux) const;
  bool isSameGridMeshOldStyle(const DbGrid* dbaux) const;
  bool isSameGridRotation(const DbGrid& dbaux) const;
  bool isSameGridRotationOldStyle(const DbGrid* dbaux) const;
  bool isGridRotated() const;

  int  getNTotal() const { return _grid.getNTotal(); }
  double getCellSize() const { return _grid.getCellSize(); }
  double getVolume() const { return _grid.getVolume(); }
  VectorDouble getExtends() const { return _grid.getExtends(); }
  double getExtend(int idim) const { return _grid.getExtend(idim); }

  int  getNX(int idim) const { return _grid.getNX(idim); }
  VectorInt getNXs() const { return _grid.getNXs(); }
  double getDX(int idim) const { return _grid.getDX(idim); }
  VectorDouble getDXs() const { return _grid.getDXs(); }
  double getX0(int idim) const { return _grid.getX0(idim); }
  VectorDouble getX0s() const { return _grid.getX0s(); }
  double getAngle(int idim) const { return _grid.getRotAngle(idim); }
  VectorDouble getAngles() const { return _grid.getRotAngles(); }
  VectorDouble getRotMat() const { return _grid.getRotMat(); }
  void setNX(int idim, int value) { _grid.setNX(idim, value); }
  void setX0(int idim, double value) { _grid.setX0(idim, value); }
  void setDX(int idim, double value) { _grid.setDX(idim, value); }
  VectorDouble getGridAxis(int idim) const { return _grid.getAxis(idim); }
  VectorDouble getCoordinateFromCorner(const VectorInt& icorner) const
  {
    return _grid.getCoordinatesByCorner(icorner);
  }
  int coordinateToRank(const VectorDouble& coor, double eps = EPSILON6) const
  {
    return _grid.coordinateToRank(coor,eps);
  }
  int indiceToRank(const VectorInt& indice) const
  {
    return _grid.indiceToRank(indice);
  }
  void rankToIndice(int node, VectorInt& indice, bool minusOne = false) const
  {
    _grid.rankToIndice(node,indice, minusOne);
  }
  void rankToCoordinate(int rank,
                        VectorDouble& coor,
                        const VectorDouble& percent = VectorDouble()) const
  {
    _grid.rankToCoordinatesInPlace(rank, coor, percent);
  }
  int getMirrorIndex(int idim, int ix) const
  {
    return _grid.getMirrorIndex(idim, ix);
  }

  VectorVectorDouble getSlice(const String& name,
                              int pos,
                              int indice,
                              bool useSel = false) const;
  VectorDouble getOneSlice(const String& name,
                           int posx,
                           int posy,
                           const VectorInt& corner,
                           bool useSel = false) const;
  int assignGridColumn(const String& name,
                       int idim,
                       int rank,
                       double value,
                       bool useSel = false);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "DbGrid"; }

private:
  void _createCoordinatesGrid(int icol0);


private:
  Grid _grid;                //!< Grid characteristics
};

