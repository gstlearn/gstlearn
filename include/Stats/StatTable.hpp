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

#include "Basic/Vector.hpp"
#include "Basic/ASerializable.hpp"

/**
 * Stores the multivariate statistics
 * - each line corresponds to a sample
 * - each column corresponds to a variable
 * The organization stands as a vector (variables) or samples.
 * This allows adding the statistics for all variables for a new sample
 */
class StatTable: public ASerializable
{
public:
  StatTable();
  StatTable(const StatTable &m);
  StatTable& operator= (const StatTable &m);
  virtual ~StatTable();

public:
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;

  bool isEmpty() const { return _stats.empty(); }
  int getRowNumber() const;
  int getColNumber() const;

private:
  VectorVectorDouble _stats;
};
