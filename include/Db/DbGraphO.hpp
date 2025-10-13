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

#include "Matrix/MatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "gstlearn_export.hpp"

#include "Db/Db.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/ICloneable.hpp"

namespace gstlrn 
{

/**
 * \brief
 * Class containing the Data Information organized as a Oriented Graph
 *
 * This class is derived from the Db class, with a specific decoration: samples
 * or 'nodes' are connected via oriented 'arcs' within a Graph.
 *
 * Note that this particular Db does not allow the modification of the sample
 * number by addition or deletion.
 *
 * The graph is stored as a non-symmetric sparse square matrix.
 * Its number of rows and columns (regardless of sparcity) must be equal to the
 * number of samples contained in the Db)
 * By convention, a dowstream arc starts from a sample (whose rank is given by
 * row index) and ends to another sample (whose rank is given by column index).
 */
class GSTLEARN_EXPORT DbGraphO: public Db
{
public:
  DbGraphO();
  DbGraphO(const DbGraphO& r);
  DbGraphO& operator=(const DbGraphO& r);
  virtual ~DbGraphO();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbGraphO)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isLine() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  Id resetFromSamples(Id nech,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       NF_Triplet& NF_arcs,
                       const VectorString& names        = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       bool flagAddSampleRank           = true);
  Id resetFromMatrix(Id nech,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const MatrixSparse& MatArcs,
                       const VectorString& names        = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       bool flagAddSampleRank           = true);
  static DbGraphO* createFromSamples(Id nech,
                                     const ELoadBy& order,
                                     const VectorDouble& tab,
                                     NF_Triplet& NF_arcs,
                                     const VectorString& names = VectorString(),
                                     const VectorString& locatorNames = VectorString(),
                                     bool flagAddSampleRank           = true);

  static DbGraphO* createFromMatrix(Id nech,
                                    const ELoadBy& order,
                                    const VectorDouble& tab,
                                    const MatrixSparse& MatArcs,
                                    const VectorString& names        = VectorString(),
                                    const VectorString& locatorNames = VectorString(),
                                    bool flagAddSampleRank           = true);

  static DbGraphO* createFromNF(const String& NFFilename, bool verbose = true);

  Id getNArc() const;
  Id getNNode() const;
  VectorDouble getArc(Id iarc, Id idim) const;
  double getArcValue(Id iarc) const;
  VectorInt getOrderDown(Id node = 0) const;
  VectorDouble getCumulDown(Id node) const;
  VectorInt getIndicesNextDown(Id node = 0) const;
  VectorInt getIndicesNextUp(Id node = 0) const;
  bool isEndDown(Id node = 0) const;
  bool isEndUp(Id node = 0) const;
  bool areConnected(Id node1, Id node2) const;
  VectorInt getEndsDown() const;
  VectorInt getEndsUp() const;
  VectorInt getOrphans() const;

  const MatrixSparse& getMatArcs() const { return _downArcs; }

  void setArcLine(const VectorInt& nodes, double value = 1.);

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "DbGraphO"; }

private:
  Id  _arcLinkage(NF_Triplet& NF_arcs, Id nech);
  bool _isValidArcRank(Id iarc) const;
  bool _isValidNode(Id node) const;
  void _checkForceDimension(Id nech);
  static VectorInt _getNoneZeroIndices(const VectorDouble& v);
  static VectorInt _getRanks(const VectorDouble& v);
  static void _updateOrder(Id rank, const VectorDouble& v, VectorInt& order);
  static void _updateCumul(Id rank, const VectorDouble& v, VectorDouble& cumul);
  void _iterateCumul(const VectorInt& inds,
                     VectorDouble& cumul,
                     VectorDouble& v1,
                     VectorDouble& v2) const;

private:
  // Information on arcs connecting nodes of the Db:
  // a sample 'irow' is connected to any other sample 'icol' if the
  // corresponding matrix value is non zero.
  MatrixSparse _downArcs;
};
}
