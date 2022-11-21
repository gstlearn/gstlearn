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

#include "gstlearn_export.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "RuleStringFormat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"

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
  RuleShadow(const RuleShadow& r);
  RuleShadow& operator=(const RuleShadow& r);
  virtual ~RuleShadow();

  String displaySpecific() const override;

  int particularities(Db *db,
                      const Db *dbprop,
                      Model *model,
                      int flag_grid_check,
                      int flag_stat)  const override;
  int gaus2facData(PropDef *propdef,
                   Db *dbin,
                   Db *dbout,
                   int *flag_used,
                   int ipgs,
                   int isimu,
                   int nbsimu) override;
  int gaus2facResult(PropDef *propdef,
                     Db *dbout,
                     int *flag_used,
                     int ipgs,
                     int isimu,
                     int nbsimu) const override;
  int evaluateBounds(PropDef *propdef,
                     Db *dbin,
                     Db *dbout,
                     int isimu,
                     int igrf,
                     int ipgs,
                     int nbsimu) const override;


  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope()  const { return _slope;  }
  double getDMax()   const { return _dMax;   }
  double getTgte()   const { return _tgte;   }
  double getIncr()   const { return _incr;   }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(int idim) const { return _shift[idim]; }

protected:
  /// Interface for ASerializable
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  String _getNFName() const override { return "RuleShadow"; }

private:
  void _st_shadow_max(const Db *dbprop,
                      int flag_stat,
                      double *sh_dsup_max,
                      double *sh_down_max) const;
  double _st_grid_eval(DbGrid *dbgrid,
                       int isimu,
                       int icase,
                       int nbsimu,
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
