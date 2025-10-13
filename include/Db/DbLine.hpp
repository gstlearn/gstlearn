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

#include "Db/Db.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/ICloneable.hpp"

namespace gstlrn
{
/**
 * \brief
 * Class containing the Data Information organized as a set of Lines
 *
 * This class is derived from the Db class, with a specific decoration: samples
 * within the Db are organized sequenially along lines.
 *
 * Note that this particular Db does not allow the modification of the sample
 * number by addition or deletion.
 *
 * The Line decoration provides a vector giving the number of samples per line.
 * Within one line, samples are ordered sequentailly (the order must not be modified).
 */
class GSTLEARN_EXPORT DbLine: public Db
{
public:
  DbLine();
  DbLine(const DbLine& r);
  DbLine& operator=(const DbLine& r);
  virtual ~DbLine();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbLine)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isLine() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  Id resetFromSamples(Id nech,
                       const ELoadBy& order,
                       const VectorDouble& tab,
                       const VectorInt& lineCounts,
                       const VectorString& names        = VectorString(),
                       const VectorString& locatorNames = VectorString(),
                       bool flagAddSampleRank           = true);
  Id resetFromSamplesById(Id nech,
                           const ELoadBy& order,
                           const VectorDouble& tab,
                           const VectorInt& lineIds,
                           const VectorInt& ranksPerId,
                           const VectorString& names        = VectorString(),
                           const VectorString& locatorNames = VectorString(),
                           bool flagAddSampleRank           = true);

  static DbLine* createFromSamples(Id nech,
                                   const ELoadBy& order,
                                   const VectorDouble& tab,
                                   const VectorInt& lineCounts,
                                   const VectorString& names = VectorString(),
                                   const VectorString& locatorNames = VectorString(),
                                   bool flagAddSampleRank = true);
  static DbLine* createFromSamplesById(Id nech,
                                       const ELoadBy& order,
                                       const VectorDouble& tab,
                                       const VectorInt& lineIds,
                                       const VectorInt& ranksPerId,
                                       const VectorString& names = VectorString(),
                                       const VectorString& locatorNames = VectorString(),
                                       bool flagAddSampleRank           = true);
  static DbLine* createFromNF(const String& NFFilename, bool verbose = true);
  static DbLine* createFillRandom(Id ndim,
                                  Id nbline,
                                  Id nperline,
                                  double deltaX             = 5.0,
                                  const VectorDouble& delta = VectorDouble(),
                                  double unifDelta          = 0.3,
                                  Id seed                  = 13422);
  static DbLine* createVerticalFromGrid(const DbGrid& grid,
                                        const VectorString& names,
                                        const VectorInt& xranks,
                                        const VectorInt& yranks,
                                        Id byZ = 1);
  static DbLine* createMarkersFromGrid(const DbGrid& grid,
                                       const String& name,
                                       const VectorInt& xranks,
                                       const VectorInt& yranks,
                                       const VectorDouble& cuts);
  Db* createStatToHeader() const;

  Id getNLine() const;
  Id getNSamplePerLine(Id iline) const;
  double getLineLength(Id iline) const;
  Id getNTotal() const;
  Id getLineBySample(Id iech) const;
  VectorDouble _getHeaderCoordinate(Id idim) const;
  VectorDouble getCoordinatesPerLine(Id iline, Id idim) const;
  VectorDouble getLineLengths() const;
  Id getLineSampleRank(Id iline, Id isample) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "DbLine"; }

private:
  Id _lineLinkage(const VectorInt& lineCounts);
  Id _lineLinkageById(const VectorInt& linesId, const VectorInt& ranksPerId);
  bool _isLineNumberValid(Id iline) const;

private:
  // Information on addresses within the Db, per Line:
  // - first dimension: Number of Lines
  // - second dimension: Number of addresses (within Db) per Line
  VectorVectorInt _lineAdds;
};
}