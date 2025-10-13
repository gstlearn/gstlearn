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
#include "LithoRule/Rule.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class PropDef;
class DbGrid;

class GSTLEARN_EXPORT RuleShadow: public Rule
{
public:
  RuleShadow();
  RuleShadow(double slope,
             double sh_dsup,
             double sh_down,
             const VectorDouble& shift);
  RuleShadow(const RuleShadow& m);
  RuleShadow& operator=(const RuleShadow& m);
  virtual ~RuleShadow();

  String displaySpecific() const override;

  Id particularities(Db *db,
                      const Db *dbprop,
                      Model *model,
                      Id flag_grid_check,
                      Id flag_stat)  const override;
  Id gaus2facData(PropDef *propdef,
                   Db *dbin,
                   Db *dbout,
                   Id *flag_used,
                   Id ipgs,
                   Id isimu,
                   Id nbsimu) override;
  Id gaus2facResult(PropDef *propdef,
                     Db *dbout,
                     Id *flag_used,
                     Id ipgs,
                     Id isimu,
                     Id nbsimu) const override;
  Id evaluateBounds(PropDef *propdef,
                     Db *dbin,
                     Db *dbout,
                     Id isimu,
                     Id igrf,
                     Id ipgs,
                     Id nbsimu) const override;


  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope()  const { return _slope;  }
  double getDMax()   const { return _dMax;   }
  double getTgte()   const { return _tgte;   }
  double getIncr()   const { return _incr;   }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(Id idim) const { return _shift[idim]; }

protected:
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "RuleShadow"; }

private:
  void _st_shadow_max(const Db *dbprop,
                      Id flag_stat,
                      double *sh_dsup_max,
                      double *sh_down_max) const;
  double _st_grid_eval(DbGrid *dbgrid,
                       Id isimu,
                       Id icase,
                       Id nbsimu,
                       VectorDouble& xyz0) const;
  void _normalizeShift();

private:
  double _shDsup;      /* Upper limit */
  double _shDown;      /* Downwards limit */
  double _slope;       /* Slope used for shadow option */
  VectorDouble _shift; /* Shadow or translation orientation */

  mutable double _dMax;
  mutable double _tgte;
  mutable double _incr;
  mutable VectorDouble _xyz;
  mutable VectorInt    _ind1;
  mutable VectorInt    _ind2;
};
}