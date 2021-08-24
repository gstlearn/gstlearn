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

#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"
#include "geoslib_enum.h"

class PropDef;

class RuleShadow: public Rule
{
public:
  RuleShadow(double slope,
             double sh_dsup,
             double sh_down,
             const VectorDouble& shift);

  RuleShadow(const RuleShadow& r);
  RuleShadow& operator=(const RuleShadow& r);
  virtual ~RuleShadow();

  int deSerializeSpecific() override;
  void serializeSpecific() const override;
  String displaySpecific(int flagProp, int flagThresh) const override;

  int particularities(Db *db,
                      const Db *dbprop,
                      Model *model,
                      int flag_grid_check,
                      int flag_stat) override;
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
                     int nbsimu) override;
  int evaluateBounds(PropDef *propdef,
                     Db *dbin,
                     Db *dbout,
                     int isimu,
                     int igrf,
                     int ipgs,
                     int nbsimu) override;

  double st_grid_eval(Db *dbgrid,
                      int isimu,
                      int icase,
                      int nbsimu,
                      VectorDouble& xyz0);

  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope() const  { return _slope;  }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(int idim) const { return _shift[idim]; }

private:
  void _st_shadow_max(const Db *dbprop,
                      int flag_stat,
                      double *sh_dsup_max,
                      double *sh_down_max);

private:
  double _shDsup;      /* Upper limit */
  double _shDown;      /* Downwards limit */
  double _slope;       /* Slope used for shadow option */
  double _dMax;        /* Longest distance for shadow */
  double _tgte;        /* Tangent of the slope */
  double _incr;        /* Increments used for creating replicates */
  VectorDouble _shift; /* Shadow or translation orientation */
  mutable VectorDouble _xyz;
  VectorInt    _ind1;
  VectorInt    _ind2;
};
