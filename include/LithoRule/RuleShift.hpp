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

class RuleShift: public Rule
{
public:
  RuleShift(const VectorInt& nodes, const VectorDouble& shift);
  RuleShift(const VectorString& nodnames,const VectorDouble& shift);
  RuleShift(int nfacies, const VectorDouble& shift);
  RuleShift(const VectorInt& n_type,
            const VectorInt& n_facs,
            const VectorDouble& shift);
  RuleShift(const RuleShift& r);
  RuleShift& operator=(const RuleShift& r);
  virtual ~RuleShift();

  int deSerializeSpecific() override;
  void serializeSpecific() const override;
  String displaySpecific(int flagProp, int flagThresh) const override;

  int particularities(Db *db,
                      const Db *dbprop,
                      Model *model,
                      int flag_grid_check,
                      int flag_stat) override;
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


  bool checkModel(const Model* model, int nvar = 0) const;

  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope() const  { return _slope;  }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(int idim) const { return _shift[idim]; }

private:
  int _st_shift_on_grid(Db *db, int ndim, int flag_grid_check);

private:
  double _shDsup;      /* Upper limit */
  double _shDown;      /* Downwards limit */
  double _slope;       /* Slope used for shadow option */
  double _tgte;        /* Tangent of the slope */
  double _incr;        /* Increments used for creating replicates */
  VectorDouble _shift; /* Shadow or translation orientation */
  mutable VectorDouble _xyz;
  VectorInt    _ind1;
  VectorInt    _ind2;
};
