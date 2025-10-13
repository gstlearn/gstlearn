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

#include "Basic/VectorNumT.hpp"
#include "LithoRule/Rule.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class DbGrid;

class GSTLEARN_EXPORT RuleShift: public Rule
{
public:
  RuleShift();
  RuleShift(const RuleShift& m);
  RuleShift& operator=(const RuleShift& m);
  virtual ~RuleShift();

  String displaySpecific() const override;

  Id resetFromNodes(const VectorInt& nodes, const VectorDouble& shift);
  Id resetFromNames(const VectorString& nodnames, const VectorDouble& shift);
  Id resetFromFaciesCount(Id nfacies, const VectorDouble& shift);
  Id resetFromNumericalCoding(const VectorInt& n_type,
                               const VectorInt& n_facs,
                               const VectorDouble& shift);

  static RuleShift* createFromNodes(const VectorInt& nodes,
                                    const VectorDouble& shift);
  static RuleShift* createFromNames(const VectorString& nodnames,
                                    const VectorDouble& shift);
  static RuleShift* createFromFaciesCount(Id nfacies,
                                          const VectorDouble& shift);
  static RuleShift* createFromNumericalCoding(const VectorInt& n_type,
                                              const VectorInt& n_facs,
                                              const VectorDouble& shift);

  Id particularities(Db* db,
                      const Db* dbprop,
                      Model* model,
                      Id flag_grid_check,
                      Id flag_stat) const override;
  Id gaus2facResult(PropDef* propdef,
                     Db* dbout,
                     Id* flag_used,
                     Id ipgs,
                     Id isimu,
                     Id nbsimu) const override;
  Id evaluateBounds(PropDef* propdef,
                     Db* dbin,
                     Db* dbout,
                     Id isimu,
                     Id igrf,
                     Id ipgs,
                     Id nbsimu) const override;

  bool checkModel(const Model* model, Id nvar = 0) const override;

  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope() const { return _slope; }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(Id idim) const { return _shift[idim]; }

protected:
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "RuleShift"; }

private:
  Id _st_shift_on_grid(Db* db, Id ndim, Id flag_grid_check) const;

private:
  double _shDsup;      /* Upper limit */
  double _shDown;      /* Downwards limit */
  double _slope;       /* Slope used for shadow option */
  VectorDouble _shift; /* Shadow or translation orientation */

  mutable double _incr;
  mutable VectorDouble _xyz;
  mutable VectorInt _ind1;
  mutable VectorInt _ind2;
};
} // namespace gstlrn