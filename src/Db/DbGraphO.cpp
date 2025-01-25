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
#include "Matrix/NF_Triplet.hpp"
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/DbGraphO.hpp"
#include "Db/DbStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Stats/Classical.hpp"

#include <math.h>

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

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
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

DbGraphO* DbGraphO::createFromSamples(int nech,
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

DbGraphO* DbGraphO::createFromMatrix(int nech,
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
  //message("arcs\n");
  //message("nrows=%d ncols=%d\n", MatArcs.getNRows(), MatArcs.getNCols());
  return dbgraphO;
}

int DbGraphO::_arcLinkage(NF_Triplet& NF_arcs, int nech)
{
  NF_arcs.force(nech, nech);
  _downArcs.resetFromTriplet(NF_arcs);
  return 0;
}

void DbGraphO::_checkForceDimension(int nech)
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
 * @return int Error returned code
 */
int DbGraphO::resetFromSamples(int nech,
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
 * @return int Error returned code
 */
int DbGraphO::resetFromMatrix(int nech,
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

bool DbGraphO::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  int narcs = 0;
  VectorString locators;
  VectorString names;
  VectorDouble values;
  VectorDouble allvalues;

    /* Initializations */

    bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);

  // Reading the set of arcs for the Oriented Graph organization

  NF_Triplet nft;
  ret = ret && _recordRead<int>(is, "Number of arcs", narcs);
  VectorDouble tab(3);
  for (int i = 0; i < narcs; i++)
  {
    ret = ret && _recordReadVec<double>(is, "", tab, 3);
    nft.add((int) tab[0], (int) tab[1], tab[2]);
  }
  _downArcs.resetFromTriplet(nft);

  // Writing the set of addresses for Line organization

  ret = ret && Db::_deserialize(is, verbose);

  return ret;
}

bool DbGraphO::_serialize(std::ostream& os, bool verbose) const
{
  bool ret       = true;

  /* Writing the header */

  ret            = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  // Writing the set of arcs for the Oriented Graph organization

  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  ret = ret && _recordWrite<int>(os, "Number of arcs", getArcNumber());

  VectorDouble tab(3);
  for (int i = 0, n = getArcNumber(); i < n; i++)
  {
    tab[0] = nft.getRow(i);
    tab[1] = nft.getCol(i);
    tab[2] = nft.getValue(i);
    ret    = ret && _recordWriteVec<double>(os, "", tab);
  }

  /* Writing the tail of the file */

  ret = ret && Db::_serialize(os, verbose);

  return ret;
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (Db format)
 * @param verbose         Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbGraphO* DbGraphO::createFromNF(const String& neutralFilename, bool verbose)
{
  DbGraphO* dbgraphO = new DbGraphO;
  std::ifstream is;
  bool success = false;
  if (dbgraphO->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = dbgraphO->deserialize(is, verbose);
  }
  if (! success)
  {
    delete dbgraphO;
    dbgraphO = nullptr;
  }
  return dbgraphO;
}

/**
 * @brief Check if the contents of private member of this class is compatible
 * with the number of samples stored in the Db
 * @return true if everything is OK; false if a problem occurs
 */
bool DbGraphO::isConsistent() const
{
  // Check on the count of addresses
  int nech = getNodeNumber();
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
  for (int irow = 0, nrows = _downArcs.getNRows(); irow < nrows; irow++)
    for (int icol = 0, ncols = _downArcs.getNCols(); icol < ncols; icol++)
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

int DbGraphO::getArcNumber() const
{
  return _downArcs.getNonZeros();
}

int DbGraphO::getNodeNumber() const
{
  return getNSample();
}

bool DbGraphO::_isValidArcRank(int iarc) const
{
  if (iarc < 0)
  {
    messerr("Argument 'iarc' (%d) should not be negative", iarc);
    return false;
  }
  int narcs = getArcNumber();
  if (iarc >= narcs)
  {
    messerr("Argument 'iarc' (%d) should be smaller than Number of arcs (%d)", iarc, narcs);
    return false;
  }
  return true;
}

bool DbGraphO::_isValidNode(int node) const
{
  if (node < 0)
  {
    messerr("Argument 'node' (%d) should not be negative", node);
    return false;
  }
  int nodeNumber = getNodeNumber();
  if (node >= nodeNumber)
  {
    messerr("Argument 'node' (%d) should be smaller than Number of Samples (%d)",
      node, nodeNumber);
    return false;
  }
  return true;
}

VectorDouble DbGraphO::getArc(int iarc, int idim) const
{
  if (!_isValidArcRank(iarc)) return VectorDouble();
  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  VectorDouble vec(2);
  vec[0] = getCoordinate(nft.getRow(iarc), idim);
  vec[1] = getCoordinate(nft.getCol(iarc), idim);
  return vec;
}

double DbGraphO::getArcValue(int iarc) const
{
  if (!_isValidArcRank(iarc)) return TEST;
  NF_Triplet nft = _downArcs.getMatrixToTriplet();
  return nft.getValue(iarc);
}

void DbGraphO::_updateOrder(int rank, const VectorDouble& v, VectorInt& order)
{
  int nech = (int)v.size();

  for (int iech = 0; iech < nech; iech++)
  {
    if (v[iech] > 0.) order[iech] = MAX(order[iech], rank);
  }
}

void DbGraphO::_updateCumul(int rank, const VectorDouble& v, VectorDouble& cumul)
{
  int nech = (int)v.size();

  for (int iech = 0; iech < nech; iech++)
  {
    if (v[iech] <= 0.) continue;
    cumul[iech] = MAX(cumul[iech], v[iech] + cumul[rank]);
  }
}

VectorInt DbGraphO::_getRanks(const VectorDouble& v)
{
  VectorInt retvec;
  for (int iech = 0, nech = (int)v.size(); iech < nech; iech++)
  {
    if (v[iech] > 0.) retvec.push_back(iech);
  }
  return retvec;
}

VectorInt DbGraphO::getIndicesNextDown(int node) const
{
  if (!_isValidNode(node)) return VectorInt();
  int nech = getNodeNumber();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v1[node] = 1.;
  _downArcs.prodVecMatInPlace(v1, v2);
  return _getRanks(v2);
}

VectorInt DbGraphO::getIndicesNextUp(int node) const
{
  if (!_isValidNode(node)) return VectorInt();
  int nech = getNodeNumber();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v1[node] = 1.;
  _downArcs.prodMatVecInPlace(v1, v2);
  return _getRanks(v2);
}

bool DbGraphO::isEndDown(int node) const
{
  if (!_isValidNode(node)) return false;
  VectorInt inds = getIndicesNextDown(node);
  return inds.empty();
}

bool DbGraphO::isEndUp(int node) const
{
  if (!_isValidNode(node)) return false;
  VectorInt inds = getIndicesNextUp(node);
  return inds.empty();
}

bool DbGraphO::areConnected(int node1, int node2) const
{
  if (!_isValidNode(node1)) return false;
  if (!_isValidNode(node2)) return false;
  int nech = getNodeNumber();

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
  int nech = getNodeNumber();
  for (int iech = 0; iech < nech; iech++)
    if (isEndDown(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getEndsUp() const
{
  VectorInt vec;
  int nech = getNodeNumber();
  for (int iech = 0; iech < nech; iech++)
    if (isEndUp(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getOrphans() const
{
  VectorInt vec;
  int nech = getNodeNumber();
  for (int iech = 0; iech < nech; iech++)
    if (isEndUp(iech) && isEndDown(iech)) vec.push_back(iech);
  return vec;
}

VectorInt DbGraphO::getOrderDown(int node) const
{
  if (!_isValidNode(node)) return VectorInt();
  int nech = getNodeNumber();

  VectorInt order(nech,0);
  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  v2[node] = 1.;

  int rank = 1;
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
  int nech = (int)v.size();

  VectorInt vall;
  for (int iech = 0; iech < nech; iech++)
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
  for (int ind = 0, nind = (int)inds.size(); ind < nind; ind++)
  {
    int rank = inds[ind];
    v1.fill(0.);
    v1[rank] = 1.;
    _downArcs.prodVecMatInPlace(v1, v2);
    _updateCumul(rank, v2, cumul);
    VectorInt indbis = _getNoneZeroIndices(v2);
    _iterateCumul(indbis, cumul, v1, v2);
  }
}

VectorDouble DbGraphO::getCumulDown(int node) const
{
  if (!_isValidNode(node)) return VectorDouble();
  int nech = getNodeNumber();

  VectorDouble v1(nech, 0.);
  VectorDouble v2(nech, 0.);
  VectorDouble cumul(nech, 0);

  _iterateCumul({node}, cumul, v1, v2);
  return cumul;
}

void DbGraphO::setArcLine(const VectorInt& nodes, double value)
{
  int number = (int)nodes.size();
  for (int i = 1; i < number; i++)
  {
    int i1 = nodes[i - 1];
    int i2 = nodes[i];
    _downArcs.setValue(i1, i2, value);
  }
}
