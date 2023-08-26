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
#include "geoslib_d.h"

#include "Enum/ELoadBy.hpp"

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
class EMorpho;
class ANeigh;
class SpaceTarget;

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
  void getCoordinatesPerSampleInPlace(int iech, VectorDouble& coor, bool flag_rotate=true) const override;
  double getUnit(int idim = 0) const override;
  int getNDim() const override;
  bool mayChangeSampleNumber() const override { return false; }
  void resetDims(int ncol, int nech) override;

  static DbGrid* createFromNF(const String& neutralFilename,
                              bool verbose = true);
  int reset(const VectorInt& nx,
            const VectorDouble& dx = VectorDouble(),
            const VectorDouble& x0 = VectorDouble(),
            const VectorDouble& angles = VectorDouble(),
            const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
            const VectorDouble& tab = VectorDouble(),
            const VectorString& names = VectorString(),
            const VectorString& locatorNames = VectorString(),
            int flag_add_rank = 1,
            bool flag_add_coordinates = true);
  int resetCoveringDb(const Db* db,
                      const VectorInt& nx = VectorInt(),
                      const VectorDouble& dx = VectorDouble(),
                      const VectorDouble& x0 = VectorDouble(),
                      const VectorDouble& margin = VectorDouble());
  int resetFromPolygon(Polygons* polygon,
                       const VectorInt& nodes,
                       const VectorDouble& dcell,
                       int flag_add_rank = 1);

  static DbGrid* create(const VectorInt& nx,
                        const VectorDouble& dx = VectorDouble(),
                        const VectorDouble& x0 = VectorDouble(),
                        const VectorDouble& angles = VectorDouble(),
                        const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
                        const VectorDouble& tab = VectorDouble(),
                        const VectorString& names = VectorString(),
                        const VectorString& locatorNames = VectorString(),
                        int flag_add_rank = 1,
                        bool flag_add_coordinates = true);
  static DbGrid* createCoveringDb(const Db* dbin,
                                  const VectorInt& nx = VectorInt(),
                                  const VectorDouble& dx = VectorDouble(),
                                  const VectorDouble& x0 = VectorDouble(),
                                  const VectorDouble& margin = VectorDouble());
  static DbGrid* createFromPolygon(Polygons* polygon,
                                   const VectorInt& nodes,
                                   const VectorDouble& dcell,
                                   int flag_add_rank = 1);
  static DbGrid* createCoarse(DbGrid *dbin,
                              const VectorInt &nmult,
                              int flag_cell = 1,
                              int flag_add_rank = 1);
  static DbGrid* createRefine(DbGrid *dbin,
                              const VectorInt &nmult,
                              int flag_cell = 1,
                              int flag_add_rank = 1);
  static DbGrid* createFromGridExtend(const DbGrid& gridIn,
                                      const VectorString &tops,
                                      const VectorString &bots,
                                      const VectorInt &nxnew,
                                      bool verbose = false,
                                      double eps = EPSILON3);
  static DbGrid* createFromGridShrink(const DbGrid& gridIn,
                                      const VectorInt& deletedRanks);
  static DbGrid* createGrid2D(const ELoadBy &order,
                              int nx,
                              int ny,
                              double x0 = 0.,
                              double y0 = 0.,
                              double dx = 1.,
                              double dy = 1.,
                              double angle = 0.,
                              int flag_add_rank = 1,
                              const VectorDouble &tab = VectorDouble());

  DbGrid* coarsify(const VectorInt &nmult);
  DbGrid* refine(const VectorInt &nmult);
  static bool migrateAllVariables(Db *dbin, Db *dbout, int flag_add_rank = 1);

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
  bool isSameGridRotation(const DbGrid& dbaux) const;
  bool isGridRotated() const;

  int  getNTotal() const { return _grid.getNTotal(); }
  double getCellSize() const { return _grid.getCellSize(); }
  double getVolume() const { return _grid.getVolume(); }
  VectorDouble getExtends() const { return _grid.getExtends(); }
  double getExtend(int idim) const { return _grid.getExtend(idim); }

  int  getNX(int idim) const { return _grid.getNX(idim); }
  VectorInt getNXs() const { return _grid.getNXs(); }
  VectorInt getNXsExt(int ndimMax) const;
  bool hasSingleBlock() const;
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
  VectorDouble getCoordinatesByIndice(const VectorInt &indice,
                                      bool flag_rotate = true,
                                      const VectorInt& shift = VectorInt(),
                                      const VectorDouble& dxsPerCell = VectorDouble()) const
  {
    return _grid.getCoordinatesByIndice(indice, flag_rotate, shift, dxsPerCell);
  }
  int coordinateToRank(const VectorDouble &coor,
                       bool centered = false,
                       double eps = EPSILON6) const;
  VectorInt coordinateToIndices(const VectorDouble &coor,
                                bool centered = false,
                                double eps = EPSILON6) const;
  int coordinateToIndicesInPlace(const VectorDouble &coor,
                                 VectorInt &indices,
                                 bool centered = false,
                                 double eps = EPSILON6) const;
  VectorInt getCenterIndices() const
  {
    return _grid.getCenterIndices();
  }
  int indiceToRank(const VectorInt& indice) const
  {
    return _grid.indiceToRank(indice);
  }
  void rankToIndice(int node, VectorInt& indices, bool minusOne = false) const
  {
    _grid.rankToIndice(node,indices, minusOne);
  }
  void rankToCoordinatesInPlace(int rank,
                                VectorDouble &coor,
                                const VectorDouble &percent = VectorDouble()) const
  {
    _grid.rankToCoordinatesInPlace(rank, coor, percent);
  }
  VectorDouble rankToCoordinates(int rank,
                                 const VectorDouble& percent = VectorDouble()) const
  {
    return _grid.rankToCoordinates(rank,percent);
  }

  void indicesToCoordinateInPlace(const VectorInt& indice,
                                  VectorDouble& coor,
                                  const VectorDouble& percent = VectorDouble()) const
  {
    _grid.indicesToCoordinateInPlace(indice, coor, percent);
  }
  VectorDouble indicesToCoordinate(const VectorInt& indice,
                                   const VectorDouble& percent = VectorDouble()) const
  {
    return _grid.indicesToCoordinate(indice, percent);
  }
  bool sampleBelongsToCell(const VectorDouble& coor,
                           int rank,
                           const VectorDouble &dxsPerCell = VectorDouble()) const
  {
    return _grid.sampleBelongsToCell(coor, rank, dxsPerCell);
  }

  int centerCoordinateInPlace(VectorDouble &coor,
                              bool centered = false,
                              bool stopIfOut = false,
                              double eps = EPSILON6) const;

  VectorInt locateDataInGrid(const Db *data,
                             const VectorInt &rankIn = VectorInt(),
                             bool centered = false,
                             bool useSel = false) const;

  int getMirrorIndex(int idim, int ix) const
  {
    return _grid.getMirrorIndex(idim, ix);
  }

  VectorVectorDouble getSlice(const String& name,
                              int pos = 0,
                              int indice = 0,
                              bool useSel = false) const;
  VectorDouble getOneSlice(const String& name,
                           int posx = 0,
                           int posy = 1,
                           const VectorInt& corner = VectorInt(),
                           bool useSel = false) const;
  int assignGridColumn(const String& name,
                       int idim,
                       int rank,
                       double value,
                       bool useSel = false);
  VectorDouble getBlockExtensions(int node) const;
  VectorVectorDouble getCellEdges(int node = 0, bool forceGridMesh = false) const;
  VectorVectorDouble getAllCellsEdges(bool forceGridMesh = false) const;
  VectorVectorDouble getGridEdges() const;
  VectorDouble getCodir(const VectorInt& grincr) const;

  int morpho(const EMorpho &oper,
             double vmin = 0.5,
             double vmax = 1.5,
             int option = 0,
             const VectorInt &radius = VectorInt(),
             bool flagDistErode = false,
             bool verbose = false,
             const NamingConvention &namconv = NamingConvention("Morpho"));
  int smooth(ANeigh *neigh,
             int type = 1,
             double range = 1.,
             const NamingConvention &namconv = NamingConvention("Smooth"));

  int addSelectionFromDbByMorpho(Db *db,
                                 int nmin = 0,
                                 int radius = 0,
                                 int option = 0,
                                 const VectorInt &dilation = VectorInt(),
                                 bool verbose = false,
                                 const NamingConvention &namconv = NamingConvention("Morpho", false, false, true,
                                     ELoc::fromKey("SEL")));

  void getSampleAsST(int iech, SpaceTarget& P) const;

  VectorVectorDouble getDiscretizedBlock(const VectorInt &ndiscs,
                                         int iech = 0,
                                         bool flagPerCell = false,
                                         bool flagRandom = false,
                                         int seed = 132433) const;

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

