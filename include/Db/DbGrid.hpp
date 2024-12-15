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

#include "Enum/ELoadBy.hpp"

#include "Db/Db.hpp"
#include "Basic/Grid.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/ICloneable.hpp"

class Polygons;
class EMorpho;
class ANeigh;
class SpaceTarget;

/**
 * \brief
 * Class containing the Data Information organized as a Regular Grid
 *
 * This class is derived from the Db class, with a specific decoration: its samples correspond to the nodes
 * of a regular grid defined in the current space.
 *
 * A regular grid is a tessellation of n-dimensional Euclidean space by congruent parallelotopes (e.g. bricks).
 * The number of meshes along each space direction is defined by NX and the total number of samples is equal to:
 * NX[1] x NX[2] x ... x NX[NDIM] where NDIM is he space dimension.
 *
 * Note that this particular Db does not allow the modification of the sample number by addition or deletion.
 *
 * The grid decoration contains the following information:
 * - the number of meshes per direction (**NX**)
 * - the coordinates of the grid origin (**X0**) defined as the node located in first position (smallest index)
 * along each space direction
 * - the mesh value along each space direction (**DX**)
 * - the rotation angles around the grid origin (the origin is invariant through this rotation).
 *
 * Note that rotation is defined by as many angles as the space dimension. This general property is not valid:
 * - for 1D where the rotation concept does not make sense
 * - for 2D where a single angle would be sufficient: however a second dummy angle (always set to 0) is used for generality sake.
 *
 * Note that this grid decoration essentially describes the organization of the grid nodes. Then the same grid can be considered
 * as a grid of points (**punctual** grid) or a grid of cells (**block** grid). In the latter case:
 * - the origin is (currently) the lower-left corner (for 2D) or the origin cell.
 * - the extensions of the cells are considered as constant by default and equal to DX. However, they can be replaced
 * by a set of variables allowing the user to assign a variable cell extension for each cell.
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
  double getCoordinate(int iech, int idim, bool flag_rotate = true) const override;
  void getCoordinatesPerSampleInPlace(int iech, VectorDouble& coor, bool flag_rotate=true) const override;
  double getUnit(int idim = 0) const override;
  int getNDim() const override;
  bool mayChangeSampleNumber() const override { return false; }
  void resetDims(int ncol, int nech) override;
  bool isConsistent() const override;

  static DbGrid* createFromNF(const String& neutralFilename,
                              bool verbose = true);
  static DbGrid* create(const VectorInt& nx,
                        const VectorDouble& dx = VectorDouble(),
                        const VectorDouble& x0 = VectorDouble(),
                        const VectorDouble& angles = VectorDouble(),
                        const ELoadBy& order = ELoadBy::fromKey("SAMPLE"),
                        const VectorDouble& tab = VectorDouble(),
                        const VectorString& names = VectorString(),
                        const VectorString& locatorNames = VectorString(),
                        bool flagAddSampleRank = true,
                        bool flagAddCoordinates = true);
  static DbGrid* createCoveringDb(const Db* dbin,
                                  const VectorInt& nx = VectorInt(),
                                  const VectorDouble& dx = VectorDouble(),
                                  const VectorDouble& x0 = VectorDouble(),
                                  const VectorDouble& margin = VectorDouble());
  static DbGrid* createFromPolygon(Polygons* polygon,
                                   const VectorInt& nodes,
                                   const VectorDouble& dcell,
                                   bool flagAddSampleRank = true);
  static DbGrid* createCoarse(DbGrid *dbin,
                              const VectorInt &nmult,
                              bool flagCell = true,
                              bool flagAddSampleRank = true);
  static DbGrid* createRefine(DbGrid *dbin,
                              const VectorInt &nmult,
                              bool flagCell = true,
                              bool flagAddSampleRank = true);
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
                              bool flagAddSampleRank = true,
                              const VectorDouble &tab = VectorDouble());
  static DbGrid* createSqueezeAndStretchForward(const DbGrid* grid3Din,
                                                const DbGrid *surf2D,
                                                const String &nameTop,
                                                const String &nameBot,
                                                const VectorString &names,
                                                int nzout,
                                                double thickmin = 0.);
  static DbGrid* createSqueezeAndStretchBackward(const DbGrid *grid3Din,
                                                 const DbGrid *surf2D,
                                                 const String &nameTop,
                                                 const String &nameBot,
                                                 const VectorString &names,
                                                 int nzout,
                                                 double z0out,
                                                 double dzout);
  static DbGrid* createSubGrid(const DbGrid* gridin,
                               VectorVectorInt limits,
                               bool flagAddCoordinates = false);
  static DbGrid* createMultiple(DbGrid* dbin, const VectorInt& nmult, bool flagAddSampleRank);
  static DbGrid* createDivider(DbGrid* dbin, const VectorInt& nmult, bool flagAddSampleRank);

  int reset(const VectorInt& nx,
            const VectorDouble& dx           = VectorDouble(),
            const VectorDouble& x0           = VectorDouble(),
            const VectorDouble& angles       = VectorDouble(),
            const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
            const VectorDouble& tab          = VectorDouble(),
            const VectorString& names        = VectorString(),
            const VectorString& locatorNames = VectorString(),
            bool flagAddSampleRank           = true,
            bool flagAddCoordinates          = true);
  int resetCoveringDb(const Db* db,
                      const VectorInt& nx = VectorInt(),
                      const VectorDouble& dx = VectorDouble(),
                      const VectorDouble& x0 = VectorDouble(),
                      const VectorDouble& margin = VectorDouble());
  int resetFromPolygon(Polygons* polygon,
                       const VectorInt& nodes,
                       const VectorDouble& dcell,
                       bool flagAddSampleRank = true);

  DbGrid* coarsify(const VectorInt &nmult);
  DbGrid* refine(const VectorInt &nmult);
  bool migrateAllVariables(Db *dbin, bool flag_fill = true,
                           bool flag_inter = true, bool flag_ball = false,
                           bool flagAddSampleRank = true);
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
  VectorDouble getCoordinatesPerSample(int iech, bool flag_rotate = true) const
  {
    return _grid.getCoordinatesByRank(iech, flag_rotate);
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
  VectorVectorInt getLimitsFromVariableExtend(const String &nameTop,
                                              const String &nameBot,
                                              const VectorInt &dimExclude = VectorInt()) const;
  int setSelectionFromVariableExtend(const String &nameTop, const String &nameBot);
  void clean3DFromSurfaces(const VectorString& names,
                           const DbGrid* surf2D,
                           const String& nameTop = String(),
                           const String& nameBot = String(),
                           bool verbose = false);

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

  void getSampleAsSTInPlace(int iech, SpaceTarget& P) const override;

  VectorVectorDouble getDiscretizedBlock(const VectorInt &ndiscs,
                                         int iech = 0,
                                         bool flagPerCell = false,
                                         bool flagRandom = false,
                                         int seed = 132433) const;

  void getGridPileInPlace(int iuid,
                            const VectorInt &indg,
                            int idim0,
                            VectorDouble &vec) const;
  void setGridPileInPlace(int iuid,
                            const VectorInt &indg,
                            int idim0,
                            const VectorDouble &vec);
  VectorDouble getDistanceToOrigin(const VectorInt& origin,
                                   const VectorDouble& radius = VectorDouble());
  
protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "DbGrid"; }

private:
  void _createGridCoordinates(int icol0);
  void _interpolate(const DbGrid *grid3D,
                    int idim0,
                    double top,
                    double bot,
                    const VectorDouble &vecin,
                    VectorDouble &vecout) const;

private:
  Grid _grid;                //!< Grid characteristics
};
