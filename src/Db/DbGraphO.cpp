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
#include "Db/DbGraphO.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Polygon/Polygons.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
DbGraphO::DbGraphO()
  : Db()
  , _downArcs()
{
  _clear();
}

DbGraphO::DbGraphO(const DbGraphO& r)
  : Db(r)
  , _downArcs(r._downArcs)
{
}

DbGraphO& DbGraphO::operator=(const DbGraphO& r)
{
  if (this != &r)
  {
    Db::operator=(r);
    _downArcs = r._downArcs;
  }
  return *this;
}

DbGraphO::~DbGraphO()
{
}

String DbGraphO::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const auto* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Oriented Graph Characteristics");

  sstr << _toStringCommon(&dsf);

  if (dsf.matchResume())
  {
    sstr << _summaryString();
  }

  sstr << _downArcs.toString(strfmt);

  return sstr.str();
}

DbGraphO* DbGraphO::createFromSamples(Id nech,
                                      const ELoadBy& order,
                                      const VectorDouble& tab,
                                      NF_Triplet& NF_arcs,
                                      const VectorString& names,
                                      const VectorString& locatorNames,
                                      bool flagAddSampleRank)
{
  DbGraphO* dbgraphO = new DbGraphO;
  if (dbgraphO->resetFromSamples(nech, order, tab, NF_arcs, names, locatorNames,
                                 flagAddSampleRank))
  {
    messerr("Error when creating DbGraphO from Samples");
    delete dbgraphO;
    return nullptr;
  }
  return dbgraphO;
}

DbGraphO* DbGraphO::createFromMatrix(Id nech,
                                     const ELoadBy& order,
                                     const VectorDouble& tab,
                                     const MatrixSparse& MatArcs,
                                     const VectorString& names,
                                     const VectorString& locatorNames,
                                     bool flagAddSampleRank)
{
  DbGraphO* dbgraphO = new DbGraphO;
  if (dbgraphO->resetFromMatrix(nech, order, tab, MatArcs, names, locatorNames,
                                flagAddSampleRank))
  {
    messerr("Error when creating DbGraphO from Samples");
    delete dbgraphO;
    return nullptr;
  }
  // message("arcs\n");
  // message("nrows=%d ncols=%d\n", MatArcs.getNRows(), MatArcs.getNCols());
  return dbgraphO;
}

Id DbGraphO::_arcLinkage(NF_Triplet& NF_arcs, Id nech)
{
  NF_arcs.force(nech, nech);
  _downArcs.resetFromTriplet(NF_arcs);
  return 0;
}

void DbGraphO::_checkForceDimension(Id nech)
{
  if (_downArcs.getValue(nech - 1, nech - 1) > 0) return;

  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  nft.force(nech, nech);
  _downArcs.resetFromTriplet(nft);
}

/**
 * @brief Reset the contents of a DbGraphO from arguments (previous contents is
 * cleared beforehand). The line contents is provided in 'lineCounts'.
 *
 * @param nech Number of samples to be loaded
 * @param order Ordering mode used for storing in 'tab' (by column or by sample)
 * @param tab Vector containing the values to be imported
 * @param NF_arcs List of connected arcs provided as a triplet structure
 * @param names Names given to the output variables
 * @param locatorNames Name of the locators given to the output variables
 * @param flagAddSampleRank When TRUE, the 'rank' variable is added
 * @return Id Error returned code
 */
Id DbGraphO::resetFromSamples(Id nech,
                               const ELoadBy& order,
                               const VectorDouble& tab,
                               NF_Triplet& NF_arcs,
                               const VectorString& names,
                               const VectorString& locatorNames,
                               bool flagAddSampleRank)
{
  if (Db::resetFromSamples(nech, order, tab, names, locatorNames,
                           flagAddSampleRank) != 0)
    return 1;

  // Create the Arcs Linkage (forcing the dimension if needed)

  if (_arcLinkage(NF_arcs, nech) != 0) return 1;

  return (!isConsistent());
}

/**
 * @brief Reset the contents of a DbGraphO from arguments (previous contents is
 * cleared beforehand). The line contents is provided in 'lineCounts'.
 *
 * @param nech Number of samples to be loaded
 * @param order Ordering mode used for storing in 'tab' (by column or by sample)
 * @param tab Vector containing the values to be imported
 * @param MatArcs Sparse Matrix giving the arcs
 * @param names Names given to the output variables
 * @param locatorNames Name of the locators given to the output variables
 * @param flagAddSampleRank When TRUE, the 'rank' variable is added
 * @return Id Error returned code
 */
Id DbGraphO::resetFromMatrix(Id nech,
                              const ELoadBy& order,
                              const VectorDouble& tab,
                              const MatrixSparse& MatArcs,
                              const VectorString& names,
                              const VectorString& locatorNames,
                              bool flagAddSampleRank)
{
  if (Db::resetFromSamples(nech, order, tab, names, locatorNames,
                           flagAddSampleRank) != 0)
    return 1;

  // Create the Arcs Linkage

  _downArcs = MatArcs;

  _checkForceDimension(nech);

  return (!isConsistent());
}

bool DbGraphO::_deserializeAscii(std::istream& is, bool verbose)
{
  Id ndim  = 0;
  Id narcs = 0;
  VectorString locators;
  VectorString names;
  VectorDouble values;
  VectorDouble allvalues;

  /* Initializations */

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Space Dimension", ndim);

  // Reading the set of arcs for the Oriented Graph organization

  NF_Triplet nft;
  ret = ret && _recordRead<Id>(is, "Number of arcs", narcs);
  VectorDouble tab(3);
  for (Id i = 0; i < narcs; i++)
  {
    ret = ret && _recordReadVec<double>(is, "", tab, 3);
    nft.add(static_cast<Id>(tab[0]), static_cast<Id>(tab[1]), tab[2]);
  }
  _downArcs.resetFromTriplet(nft);

  // Writing the set of addresses for Line organization

  ret = ret && Db::_deserializeAscii(is, verbose);

  return ret;
}

bool DbGraphO::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());

  // Writing the set of arcs for the Oriented Graph organization

  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  ret            = ret && _recordWrite<Id>(os, "Number of arcs", getNArc());

  VectorDouble tab(3);
  for (Id i = 0, n = getNArc(); i < n; i++)
  {
    tab[0] = nft.getRow(i);
    tab[1] = nft.getCol(i);
    tab[2] = nft.getValue(i);
    ret    = ret && _recordWriteVec<double>(os, "", tab);
  }

  /* Writing the tail of the file */

  ret = ret && Db::_serializeAscii(os, verbose);

  return ret;
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (Db format)
 * @param verbose    Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbGraphO* DbGraphO::createFromNF(const String& NFFilename, bool verbose)
{
  DbGraphO* dbgraphO = new DbGraphO;
  if (dbgraphO->_fileOpenAndDeserialize(NFFilename, verbose)) return dbgraphO;
  delete dbgraphO;
  return nullptr;
}

/**
 * @brief Check if the contents of private member of this class is compatible
 * with the number of samples stored in the Db
 * @return true if everything is OK; false if a problem occurs
 */
bool DbGraphO::isConsistent() const
{
  // Check on the count of addresses
  auto nech = getNNode();
  if (_downArcs.getNRows() > nech)
  {
    messerr("Number of rows of '_connectedArcs' (%d)", _downArcs.getNRows());
    messerr("must not be larger than Sample Number (%d)", nech);
    return false;
  }
  if (_downArcs.getNCols() > nech)
  {
    messerr("Number of columns of '_connectedArcs' (%d)",
            _downArcs.getNCols());
    messerr("must not be larger than Sample Number (%d)", nech);
    return false;
  }

  // Check that all arcs valuation are positive
  for (Id irow = 0, nrows = _downArcs.getNRows(); irow < nrows; irow++)
    for (Id icol = 0, ncols = _downArcs.getNCols(); icol < ncols; icol++)
    {
      double value = _downArcs.getValue(irow, icol);
      if (value < 0)
      {
        messerr("The value for Arc(%d; %d) may not be negative (%lf)", irow,
                icol, value);
        return false;
      }
    }

  return true;
}

Id DbGraphO::getNArc() const
{
  return _downArcs.getNonZeros();
}

Id DbGraphO::getNNode() const
{
  return getNSample();
}

bool DbGraphO::_isValidArcRank(Id iarc) const
{
  if (iarc < 0)
  {
    messerr("Argument 'iarc' (%d) should not be negative", iarc);
    return false;
  }
  auto narcs = getNArc();
  if (iarc >= narcs)
  {
    messerr("Argument 'iarc' (%d) should be smaller than Number of arcs (%d)", iarc, narcs);
    return false;
  }
  return true;
}

bool DbGraphO::_isValidNode(Id node) const
{
  if (node < 0)
  {
    messerr("Argument 'node' (%d) should not be negative", node);
    return false;
  }
  auto nodeNumber = getNNode();
  if (node >= nodeNumber)
  {
    messerr("Argument 'node' (%d) should be smaller than Number of Samples (%d)",
            node, nodeNumber);
    return false;
  }
  return true;
}

VectorDouble DbGraphO::getArc(Id iarc, Id idim) const
{
  if (!_isValidArcRank(iarc)) return VectorDouble();
  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  VectorDouble vec(2);
  vec[0] = getCoordinate(nft.getRow(iarc), idim);
  vec[1] = getCoordinate(nft.getCol(iarc), idim);
  return vec;
}

double DbGraphO::getArcValue(Id iarc) const
{
  if (!_isValidArcRank(iarc)) return TEST;
  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  return nft.getValue(iarc);
}

void DbGraphO::_updateOrder(Id rank, const VectorDouble& v, VectorInt& order)
{
  Id nech = static_cast<Id>(v.size());

  for (Id iech = 0; iech < nech; iech++)
  {
    if (v[iech] > 0.) order[iech] = MAX(order[iech], rank);
  }
}

void DbGraphO::_updateCumul(Id rank, const VectorDouble& v, VectorDouble& cumul)
{
  Id nech = static_cast<Id>(v.size());

  for (Id iech = 0; iech < nech; iech++)
  {
    if (v[iech] <= 0.) continue;
    cumul[iech] = MAX(cumul[iech], v[iech] + cumul[rank]);
  }
}

VectorInt DbGraphO::_getRanks(const VectorDouble& v)
{
  VectorInt retvec;
  for (Id iech = 0, nech = static_cast<Id>(v.size()); iech < nech; iech++)
  {
    if (v[iech] > 0.) retvec.push_back(iech);
  }
  return retvec;
}

VectorInt DbGraphO::getIndicesNextDown(Id node) const
{
  if (!_isValidNode(node)) return VectorInt();
  auto nech = getNNode();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v1[node] = 1.;
  _downArcs.prodVecMatInPlace(v1, v2);
  return _getRanks(v2);
}

VectorInt DbGraphO::getIndicesNextUp(Id node) const
{
  if (!_isValidNode(node)) return VectorInt();
  auto nech = getNNode();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v1[node] = 1.;
  _downArcs.prodMatVecInPlaceC(v1, v2);
  return _getRanks(v2);
}

bool DbGraphO::isEndDown(Id node) const
{
  if (!_isValidNode(node)) return false;
  VectorInt inds = getIndicesNextDown(node);
  return inds.empty();
}

bool DbGraphO::isEndUp(Id node) const
{
  if (!_isValidNode(node)) return false;
  VectorInt inds = getIndicesNextUp(node);
  return inds.empty();
}

bool DbGraphO::areConnected(Id node1, Id node2) const
{
  if (!_isValidNode(node1)) return false;
  if (!_isValidNode(node2)) return false;
  auto nech = getNNode();

  VectorInt order(nech, 0);
  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v2[node1] = 1.;

  while (VH::cumul(v2) > 0.)
  {
    v1 = v2;
    _downArcs.prodVecMatInPlace(v1, v2);
    if (v2[node2] > 0.) return true;
  }
  return false;
}

VectorInt DbGraphO::getEndsDown() const
{
  VectorInt vec;
  auto nech = getNNode();
  for (Id iech = 0; iech < nech; iech++)
    if (isEndDown(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getEndsUp() const
{
  VectorInt vec;
  auto nech = getNNode();
  for (Id iech = 0; iech < nech; iech++)
    if (isEndUp(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getOrphans() const
{
  VectorInt vec;
  auto nech = getNNode();
  for (Id iech = 0; iech < nech; iech++)
    if (isEndUp(iech) && isEndDown(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getOrderDown(Id node) const
{
  if (!_isValidNode(node)) return VectorInt();
  auto nech = getNNode();

  VectorInt order(nech, 0);
  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v2[node] = 1.;

  Id rank = 1;
  _updateOrder(rank, v2, order);
  while (VH::cumul(v2) > 0.)
  {
    rank++;
    v1 = v2;
    _downArcs.prodVecMatInPlace(v1, v2);
    _updateOrder(rank, v2, order);
  }

  return order;
}

VectorInt DbGraphO::_getNoneZeroIndices(const VectorDouble& v)
{
  Id nech = static_cast<Id>(v.size());

  VectorInt vall;
  for (Id iech = 0; iech < nech; iech++)
  {
    if (v[iech] <= 0.) continue;
    vall.push_back(iech);
  }
  return vall;
}

/**
 * @brief Local recursive function for finding cumul
 *
 * @param inds  List of indices of nodes connceted downwards
 * @param cumul Array contaiing the cumulative values per node
 * @param v1    Working array
 * @param v2    Working array
 */
void DbGraphO::_iterateCumul(const VectorInt& inds,
                             VectorDouble& cumul,
                             VectorDouble& v1,
                             VectorDouble& v2) const
{
  for (Id ind = 0, nind = static_cast<Id>(inds.size()); ind < nind; ind++)
  {
    Id rank = inds[ind];
    v1.fill(0.);
    v1[rank] = 1.;
    _downArcs.prodVecMatInPlace(v1, v2);
    _updateCumul(rank, v2, cumul);
    VectorInt indbis = _getNoneZeroIndices(v2);
    _iterateCumul(indbis, cumul, v1, v2);
  }
}

VectorDouble DbGraphO::getCumulDown(Id node) const
{
  if (!_isValidNode(node)) return VectorDouble();
  auto nech = getNNode();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  VectorDouble cumul(nech, 0);

  _iterateCumul({node}, cumul, v1, v2);
  return cumul;
}

void DbGraphO::setArcLine(const VectorInt& nodes, double value)
{
  Id number = static_cast<Id>(nodes.size());
  for (Id i = 1; i < number; i++)
  {
    Id i1 = nodes[i - 1];
    Id i2 = nodes[i];
    _downArcs.setValue(i1, i2, value);
  }
}

#ifdef HDF5
bool DbGraphO::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto dbG = SerializeHDF5::getGroup(grp, "DbGraphO");
  if (!dbG) return false;

  /* Read the grid characteristics */
  bool ret  = true;
  Id ndim  = 0;
  Id narcs = 0;

  ret = ret && SerializeHDF5::readValue(*dbG, "NDim", ndim);
  ret = ret && SerializeHDF5::readValue(*dbG, "Narcs", narcs);

  // Reading the set of arcs for the Oriented Graph organization
  auto dbgs = SerializeHDF5::getGroup(*dbG, "Arcs");
  if (!dbgs) return false;
  NF_Triplet nft;
  VectorDouble tab(3);
  for (Id i = 0; i < narcs; i++)
  {
    String locName = "Arc" + std::to_string(i);
    auto arcg      = SerializeHDF5::getGroup(*dbgs, locName);
    if (!arcg) return false;

    ret = ret && SerializeHDF5::readVec(*arcg, "Arc", tab);
    nft.add(static_cast<Id>(tab[0]), static_cast<Id>(tab[1]), tab[2]);
  }
  _downArcs.resetFromTriplet(nft);

  // Writing the set of addresses for Line organization

  ret = ret && Db::_deserializeH5(*dbG, verbose);

  return ret;
}

bool DbGraphO::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto dbG = grp.createGroup("DbGraphO");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(dbG, "NDim", getNDim());
  ret = ret && SerializeHDF5::writeValue(dbG, "Narcs", getNArc());

  // Writing the set of arcs for the Oriented Graph organization

  auto dbGs = dbG.createGroup("Arcs");

  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  VectorDouble tab(3);
  for (Id iarc = 0, narcs = getNArc(); iarc < narcs; iarc++)
  {
    String locName = "Arc" + std::to_string(iarc);
    auto arcG      = dbGs.createGroup(locName);

    tab[0] = nft.getRow(iarc);
    tab[1] = nft.getCol(iarc);
    tab[2] = nft.getValue(iarc);
    ret    = ret && SerializeHDF5::writeVec(arcG, "Arc", tab);
  }

  /* Writing the tail of the file */

  ret = ret && Db::_serializeH5(dbG, verbose);

  return ret;
}
#endif
} // namespace gstlrn