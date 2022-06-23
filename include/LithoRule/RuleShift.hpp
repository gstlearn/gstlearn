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

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"

class DbGrid;

class GSTLEARN_EXPORT RuleShift: public Rule
{
public:
  RuleShift();
  RuleShift(const RuleShift& r);
  RuleShift& operator=(const RuleShift& r);
  virtual ~RuleShift();

  String displaySpecific() const override;

  int resetFromNodes(const VectorInt& nodes, const VectorDouble& shift);
  int resetFromNames(const VectorString& nodnames,const VectorDouble& shift);
  int resetFromFaciesCount(int nfacies, const VectorDouble& shift);
  int resetFromNumericalCoding(const VectorInt& n_type,
                               const VectorInt& n_facs,
                               const VectorDouble& shift);

  static RuleShift* createFromNodes(const VectorInt& nodes,
                                    const VectorDouble& shift);
  static RuleShift* createFromNames(const VectorString& nodnames,
                                    const VectorDouble& shift);
  static RuleShift* createFromFaciesCount(int nfacies,
                                          const VectorDouble& shift);
  static RuleShift* createFromNumericalCoding(const VectorInt& n_type,
                                              const VectorInt& n_facs,
                                              const VectorDouble& shift);

  int particularities(Db *db,
                      const Db *dbprop,
                      Model *model,
                      int flag_grid_check,
                      int flag_stat) const override;
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


  bool checkModel(const Model* model, int nvar = 0) const override;

  double getShDown() const { return _shDown; }
  double getShDsup() const { return _shDsup; }
  double getSlope()  const { return _slope;  }
  const VectorDouble& getShift() const { return _shift; }
  double getShift(int idim) const { return _shift[idim]; }

protected:
  /// Interface for ASerializable
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  String _getNFName() const override { return "RuleShift"; }

private:
  int _st_shift_on_grid(Db *db, int ndim, int flag_grid_check) const;

private:
  double _shDsup;       /* Upper limit */
  double _shDown;       /* Downwards limit */
  double _slope;        /* Slope used for shadow option */
  VectorDouble _shift;  /* Shadow or translation orientation */

  mutable double       _incr;
  mutable VectorDouble _xyz;
  mutable VectorInt    _ind1;
  mutable VectorInt    _ind2;
};
