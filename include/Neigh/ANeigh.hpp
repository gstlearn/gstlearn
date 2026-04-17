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

#include "Basic/ICloneable.hpp"
#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Enum/ENeigh.hpp"

#include "Basic/ASerializable.hpp"
#include "Space/ASpaceObject.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  class Db;
  class DbGrid;
  class Ball;

  /**
   * \brief
   * Class to define the a subset of an input Data Base ('Db') called a Neighborhood.
   * This Neighborhood feature is invoked when the geostatistical processing cannot
   * handle the whole data set (usually due to core limitations) and requires a fine
   * selection of the most suitable part of the data set.
   *
   * Several implementations can be defined, such as:
   * - Unique Neighborhood: all active samples are selected (see NeighUnique)
   * - Moving Neighborhood: the sub-population essentially selects the samples close
   * enough to the target. This sub-population evolves with the location of the target,
   * hence the name of this Neighborhood feature (see NeighMoving).
   * - Bench Neighborhood: the sub-population selects all samples located in the same
   * 'bench' as the target. A bench is a portion of the space characterized sliced
   * along the highest dimension of the space (e.g. horizontal bench for 3D space) (see NeighBench).
   * - Cell Neighborhood: the sub-population selects all samples belonging to the same
   * 'cell' as the target. Obviously this feature is only valid if the target data base
   * is defined as a Grid (see NeighCell).
   * - Image Neighborhood: the sub-population selects a pattern of constant dimensions
   * centered on the target. This Neighborhood is only valid when the target and
   * the input data base are matching grids (usually they are the same file) (see NeighImage).
   *
   * Several other topics are considered as belonging to the Neighborhood selection procedure,
   * such as:
   * - Possibility to add some information, collected from the Target File, in the
   * Neighborhood calculated from the input data file: this is the Colocation option
   * - Possibility to exclude the target (or samples sharing some characteristics with
   * the Target). This is the cross-validation option.
   */
  class GSTLEARN_EXPORT ANeigh: public ASpaceObject,
                                public ASerializable,
                                public ICloneable
  {
  public:
    ANeigh(const ASpaceSharedPtr& space = ASpaceSharedPtr());
    ANeigh(const ANeigh& r);
    ANeigh& operator=(const ANeigh& r);
    virtual ~ANeigh();

    /// ASpaceObject Interface
    bool isConsistent(const ASpace* space) const override
    {
      DECLARE_UNUSED(space);
      return true;
    }

    /// ASerializable Interface
    String getNFName() const override { return "ANeigh"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// Interface for ANeigh
    virtual Id attach(const Db* dbin, const Db* dbout);
    virtual void getNeigh(Id iech_out, VectorInt& ranks) = 0;
    virtual Id getNSampleMax(const Db* db) const = 0;

    virtual bool hasChanged(Id iech_out) const
    {
      DECLARE_UNUSED(iech_out);
      return true;
    }

    virtual VectorDouble summary(Id iech_out)
    {
      DECLARE_UNUSED(iech_out);
      return VectorDouble();
    }

    virtual ENeigh getType() const = 0;

    virtual bool getFlagContinuous() const { return false; }

    static ANeigh*
      createDefaultNeighborhood(ANeigh* neigh, Db* dbin, Id maxNumber = 500);

    void displayDebug(VectorInt& ranks) const;
    void select(Id iech_out, VectorInt& ranks);

    bool isUnchanged() const { return _flagIsUnchanged; }

    void setIsChanged(bool status = false);
    void reset();

    bool getFlagXvalid() const { return _flagXvalid; }

    bool getFlagKFold() const { return _flagKFold; }

    void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }

    void setFlagKFold(bool flagKFold) { _flagKFold = flagKFold; }

    void setFlagSimu(bool flagSimu) { _flagSimu = flagSimu; }

    void setBallSearch(bool status, Id leaf_size = 10);
    void attachBall();

  protected:
    bool _isNbghMemoEmpty() const { return _nbghMemo.empty(); }

    static void _neighCompress(VectorInt& ranks);
    void _display(const VectorInt& ranks, bool flagCompress = false) const;
    bool _discardUndefined(Id iech);
    Id _xvalid(Id iech_in, Id iech_out, double eps = EPSILON9);
    bool _isDimensionValid(Id idim) const;

    Ball& getBall() { return _ball; }

  protected:
    bool _deserializeAscii(std::istream& is) override;
    bool _serializeAscii(std::ostream& os) const override;

  private:
    bool _isSameTarget(Id iech_out);
    void _checkUnchanged(Id iech_out, const VectorInt& ranks);

  protected:
    const Db* _dbin;
    const Db* _dbout;
    const DbGrid* _dbgrid; // Equivalent to dbout, defined only for grid

    Id _iechMemo;
    bool _flagSimu;
    bool _flagXvalid; /* True to suppress the target */
    bool _flagKFold; /* True to perform a KFold Cross-validation */
    bool
      _useBallSearch; /* If Neighborhood search favors Ball Tree algorithms */
    Id _ballLeafSize; /* Dimension of ultimate Leaf for Ball-Tree algorithm */

  private:
    bool _flagIsUnchanged;
    mutable VectorInt _nbghMemo;
    mutable Ball _ball;
  };
} // namespace gstlrn
