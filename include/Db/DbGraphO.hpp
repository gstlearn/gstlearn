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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isLine() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  int resetFromSamples(int nech,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const NF_Triplet& NF_arcs,
                       const VectorString& names        = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       bool flagAddSampleRank           = true);
  int resetFromMatrix(int nech,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const MatrixSparse& MatArcs,
                       const VectorString& names        = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       bool flagAddSampleRank           = true);
  static DbGraphO* createFromSamples(int nech,
                                     const ELoadBy& order,
                                     const VectorDouble& tab,
                                     const NF_Triplet& NF_arcs,
                                     const VectorString& names = VectorString(),
                                     const VectorString& locatorNames = VectorString(),
                                     bool flagAddSampleRank           = true);

  static DbGraphO* createFromMatrix(int nech,
                                    const ELoadBy& order,
                                    const VectorDouble& tab,
                                    const MatrixSparse& MatArcs,
                                    const VectorString& names        = VectorString(),
                                    const VectorString& locatorNames = VectorString(),
                                    bool flagAddSampleRank           = true);

  static DbGraphO* createFromNF(const String& neutralFilename,
                                bool verbose = true);

  int getArcNumber() const;
  int getNodeNumber() const;
  VectorDouble getArc(int iarc, int idim) const;
  double getArcValue(int iarc) const;
  VectorInt getOrderDown(int node = 0) const;
  VectorDouble getCumulDown(int node) const;
  VectorInt getIndicesNextDown(int node = 0) const;
  VectorInt getIndicesNextUp(int node = 0) const;
  bool isEndDown(int node = 0) const;
  bool isEndUp(int node = 0) const;
  bool areConnected(int node1, int node2) const;
  VectorInt getEndsDown() const;
  VectorInt getEndsUp() const;
  VectorInt getOrphans() const;

  const MatrixSparse& getMatArcs() const { return _downArcs; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,
                          bool verbose = false) const override;
  String _getNFName() const override { return "DbGraphO"; }

private:
  int _arcLinkage(const NF_Triplet& NF_arcs);
  bool _isValidArcRank(int iarc) const;
  bool _isValidNode(int node) const;
  void _checkForceDimension(int nech);
  static VectorInt _getNoneZeroIndices(const VectorDouble& v);
  static VectorInt _getRanks(const VectorDouble& v);
  static void _updateOrder(int rank, const VectorDouble& v, VectorInt& order);
  static void _updateCumul(int rank, const VectorDouble& v, VectorDouble& cumul);
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
